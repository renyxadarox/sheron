#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <strings.h>
#include <time.h>
#include <sys/timeb.h>

#ifdef BOINC
    #include "boinc_api.h"
#endif

#ifdef __linux__
	#include <sys/utsname.h>
#endif

#ifdef _WIN64
    const char* OS = "Windows 64-bit";
#elif _WIN32
    const char* OS = "Windows 32-bit";
#elif __APPLE__ || __MACH__
    const char* OS = "Mac OS X";
#elif __FreeBSD__
    const char* OS = "FreeBSD";
#else
    const char* OS = "Other";
#endif

#define APP_NAME "SHeron"
#define VERSION "1.00"
#define COPYRIGHT "Copyright 2021"
#define AUTHOR "Alexander Belogourov aka x3mEn"
#define TITLE "Heron Triangles with Square Sides search program"

// progress - режим із відображенням прогресу
int progress = 0;
// quiet - подавити вивід на stdout
int quiet = 0;
// output - відправити stdout у файл out_%
int output = 0;
// report - створити файл із статистикою задачі rep_%
int report = 0;
// skip - вважати такими, що виконані, завдання, для яких є out і немає chk
int skip = 0;
// verbose - режим із виводом результату в stderr
int verbose = 0;
// found - кількість знайдених трикутників
int found = 0;

struct timeb starttime, endtime;
uint32_t minc, minb, maxc, maxb, curc, curb;
char repfname[256] = "rep", outfname[256] = "out", chkfname[256] = "chk";
FILE * fout, * frep, * fchk;
uint64_t checksum = 0, counter;

uint32_t * primes = NULL;
uint32_t psize = 0;

#define MAX_FACTORES_CNT 63
typedef struct { uint64_t prime; uint8_t power ;} TFactor;
TFactor factors[MAX_FACTORES_CNT];
uint8_t fsize = 0;

#define max(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)

static __inline__ uint64_t string_to_u64(const char * s) {
  uint64_t i;
  char c ;
  int scanned = sscanf(s, "%" SCNu64 "%c", &i, &c);
  if (scanned == 1) return i;
  if (scanned > 1) {
    // TBD about extra data found
    return i;
    }
  // TBD failed to scan;
  return 0;
}

void u128_to_string(const __uint128_t n, char * s)
{
    uint64_t d4, d3, d2, d1, d0, q;
	const int size = 40; // floor(log10(2^128-1))
    char u[40];
    char * t = u;

	// n = d3*2^96 + d2*2^64 + d1*2^32 + d0
	// n = d3*79228162514264337593543950336 + d2*18446744073709551616 + d1*4294967296 + d0

	// n = d3*79_228162514_264337593_543950336 + d2*18_446744073_709551616 + d1*4_294967296 + d0

	// n = d3*79*10^27 + d3*228162514*10^18 + d3*264337593*10^9 + d3*543950336
	//                 + d2*       18*10^18 + d2*446744073*10^9 + d2*709551616
	//                                      + d1*        4*10^9 + d1*294967296
	//                                                          + d0*000000001

	// define constants

	const uint32_t k3 = 79;
	const uint32_t k2 = 228162514;
	const uint32_t k1 = 264337593;
	const uint32_t k0 = 543950336;

	const uint32_t l2 = 18;
	const uint32_t l1 = 446744073;
	const uint32_t l0 = 709551616;

	const uint32_t m1 = 4;
	const uint32_t m0 = 294967296;

	const uint32_t dec_unit = 1000000000;

    d0 = (uint32_t)n;
    d1 = (uint32_t)(n >> 32);
    d2 = (uint32_t)(n >> 64);
    d3 = n >> 96;

    d0 = (k0 * d3) + (l0 * d2) + (m0 * d1) + d0;
    q  = d0 / dec_unit;
    d0 = d0 % dec_unit;

    d1 = q + (k1 * d3) + (l1 * d2) + (m1 * d1);
    q  = d1 / dec_unit;
    d1 = d1 % dec_unit;

    d2 = q + (k2 * d3) + (l2 * d2);
    q  = d2 / dec_unit;
    d2 = d2 % dec_unit;

    d3 = q + (k3 * d3);
    q  = d3 / dec_unit;
    d3 = d3 % dec_unit;

    d4 = q;

    memset(t, 0, size);
	sprintf(t,"%u%.9u%.9u%.9u%.9u",(uint32_t)d4,(uint32_t)d3,(uint32_t)d2,(uint32_t)d1,(uint32_t)d0);

	// trim leading zeros
	while (*t && *t == '0') t++;
	if ( *t == 0x0 ) t--; // in case number = 0

    strcpy(s, t);
}

// Евклідів алгоритм обчислення НСД (Найбільший спільний дільник)
uint32_t gcd(uint32_t a, uint32_t b)
{
    if (!b) return a;
    uint32_t c;
    while (a)
    {
        c = a;
        a = b%a;
        b = c;
    }
    return b;
}

static __inline__ __uint128_t is_square()
{
    __uint128_t r = 1;
    uint8_t j;
    for (int8_t i = fsize-1; i >= 0; i--) {
        if (factors[i].power & 1) return 0;
        else for (j = 0; j < (factors[i].power >> 1); j++) r *= factors[i].prime;
    }
    return r;
}

void clear_factors()
{
    int i = 0;
    while (i < MAX_FACTORES_CNT && factors[i].prime) {
        factors[i] = (TFactor){0, 0};
        i++;
    }
    fsize = 0;
}

void free_primes(void)
{
    free(primes);
}

void init_primes(uint64_t max_n)
{
    uint32_t sq = ceil(sqrtl(max_n)), cb = ceil(sqrtl(sqrtf(max_n)));
    uint32_t sSize = ceil((float)sq / 128);
    uint32_t i, j;
    uint64_t * sieve;
    sieve = (uint64_t*) calloc (sSize, sizeof(uint64_t));
    if (sieve == NULL) {
#ifdef BOINC
        boinc_finish(EXIT_FAILURE);
#endif
        exit(EXIT_FAILURE);
    }
    sieve[0] = 1;
    for (i = 1; i < sSize; i++) sieve[i] = 0;
    for (i = 3; i <= cb; i += 2) {
        for (j = 3*i; j <= sq; j += 2*i) {
            sieve[j >> 7] |= ((uint64_t)1 << ((j >> 1)&63));
        }
    }
    psize = 1;
    for (i = 3; i <= sq; i += 2) {
        if (!(sieve[i >> 7]&((uint64_t)1 << ((i >> 1)&63)))) {
            psize++;
        }
    }
    primes = (uint32_t*) malloc (sizeof(uint32_t)*psize);
    if (primes == NULL) {
#ifdef BOINC
        boinc_finish(EXIT_FAILURE);
#endif
        exit(EXIT_FAILURE);
    }
    primes[0] = 2;
    psize = 1;
    for (i = 3; i <= sq; i += 2) {
        if (!(sieve[i >> 7]&((uint64_t)1 << ((i >> 1)&63)))) {
            primes[psize++] = i;
        }
    }
    free(sieve);
}

void add_factor(uint64_t f)
{
    uint32_t i = 0;
    while (factors[i].prime && factors[i].prime != f) i++;
    if (!(factors[i].prime)) {
        factors[i].prime = f;
        fsize++;
    }
    factors[i].power++;
}

void factorize(uint64_t n)
{
    uint32_t i = 0, sr = lrintl(sqrtl(n));
    int changed = 0;
    while (n > 1 && i < psize && primes[i] <= sr) {
        if (!(n % primes[i])) {
            n = n / primes[i];
            changed = 1;
            counter += 1;
            add_factor(primes[i]);
        } else {
            i++;
            if (changed) {
                sr = lrintl(sqrtl(n));
                changed = 0;
            }
        }
    }
    if (n > 1) add_factor(n);
}

int read_checkpoint(void)
{
    fchk = fopen(chkfname, "r");
    if(fchk == NULL)
        return (EXIT_FAILURE);
    char c;
    uint64_t dif;
    int scanned = fscanf(fchk, "%" SCNu32
                              ",%" SCNu32
                              ",%" SCNu64
                              ",%" SCNu64
                              ",%c"
                              , &curc
                              , &curb
                              , &checksum
                              , &dif
                              , &c);
    fclose(fchk);
    if (scanned != 4) {
#ifdef BOINC
        boinc_finish(EXIT_FAILURE);
#endif
        exit(EXIT_FAILURE);
    }
    if (!curc || !curb) return 1;
    curb += 1;

    starttime.time -= dif / 1000;
    long int millisec = (dif % 1000);
    if (starttime.millitm < millisec) {
        starttime.millitm += 1000 - millisec;
        starttime.time--;
    } else starttime.millitm -= millisec;

    return 0;
}

void save_checkpoint(uint32_t curc, uint32_t curb)
{
    fchk = fopen(chkfname, "w");
    if(fchk == NULL) {
#ifdef BOINC
        boinc_finish(EXIT_FAILURE);
#endif
        exit(EXIT_FAILURE);
    }

    struct timeb curtime;
    ftime(&curtime);
    uint64_t dif = (curtime.time - starttime.time) * 1000 + (curtime.millitm - starttime.millitm);
    fprintf(fchk,  "%" PRIu32
                  ",%" PRIu32
                  ",%" PRIu64
                  ",%" PRIu64
                ,curc
                ,curb
                ,checksum
                ,dif
           );
    fflush(fchk);
    fclose(fchk);
#if defined BOINC
	boinc_checkpoint_completed();
#endif
}

int init_task(void)
{
    if (minc > maxc || minc < 1) return 1;
    minb = min(max(minb, ceil(((float)minc + 1) / 2)), minc + 1);
    maxb = max(min(maxb, maxc + 1), ceil(((float)maxc + 1) / 2));
    if (minc == maxc) {
        if (minb >= maxb) return 1;
    }
    curc = minc;
    curb = minb;
    return 0;
}

#define PBSTR "========================================================================"
#define PBWIDTH 72
#define SCRWIDTH 80
void do_progress( double percentage )
{
    int val = (int) (percentage);
    int lpad = (int) (percentage  * (val==100?SCRWIDTH:PBWIDTH) / 100);
    int rpad = (val==100?SCRWIDTH:PBWIDTH) - lpad;
    fprintf(stderr, "\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
}

void print_usage(void)
{
#ifdef _WIN32
	char pref[3] = "";
#elif __linux__ || unix || __unix__ || __unix
	char pref[3] = "./";
#endif // __linux__
    fprintf(stderr, "Usage: %ssheron <low a> <low b> <high a> <high b> [switches]\n", pref);
    fprintf(stderr, "\t<low a>\t\tthe lower border of the biggest side\n");
    fprintf(stderr, "\t<low b>\t\tthe lower border of the middle side\n");
    fprintf(stderr, "\t<high a>\tthe higher border of the biggest side\n");
    fprintf(stderr, "\t<high b>\tthe higher border of the middle side\n");
    fprintf(stderr, "The following switches are accepted:\n");
    fprintf(stderr, "\t-q\t\tsuppress output to stdout\n");
    fprintf(stderr, "\t-p\t\tdisplay progress bar\n");
    fprintf(stderr, "\t-o\t\twrite results to output file\n");
    fprintf(stderr, "\t-r\t\twrite task stat to report file\n");
//    fprintf(stderr, "\t-s\t\tskip task if output file exists\n");
//    fprintf(stderr, "\t-v\t\tverbose mode\n");
}

int main(int argc, char** argv)
{
#if defined BOINC
	boinc_init();
#endif

#ifdef _WIN32
#elif __linux__ || unix || __unix__ || __unix
	char OS[256];
	struct utsname name;
	if(uname(&name)) exit(EXIT_FAILURE);
	sprintf(OS, "%s %s", name.sysname, name.release);
#endif // __linux__
    fprintf(stderr, "%s %s (%s)\n%s, %s\n%s\n\n", APP_NAME, VERSION, OS, COPYRIGHT, AUTHOR, TITLE);
    if (argc < 5) {
        print_usage();
#ifdef BOINC
        boinc_finish(EXIT_FAILURE);
#endif
        exit(EXIT_FAILURE);
    }

    minc = string_to_u64(argv[1]);
    minb = string_to_u64(argv[2]);
    maxc = string_to_u64(argv[3]);
    maxb = string_to_u64(argv[4]);

    if (!minc || !maxc) {
        print_usage();
#ifdef BOINC
        boinc_finish(EXIT_FAILURE);
#endif
        exit(EXIT_FAILURE);
    }

    for (int i = 5; i < argc; i++) {
        if (!strcmp(argv[i],"-q")) {quiet = 1; continue;}
        if (!strcmp(argv[i],"-p")) {progress = 1; continue;}
        if (!strcmp(argv[i],"-o")) {output = 1; continue;}
        if (!strcmp(argv[i],"-r")) {report = 1; continue;}
        if (!strcmp(argv[i],"-s")) {skip = 1; continue;}
        if (!strcmp(argv[i],"-v")) {verbose = 1; continue;}
        print_usage();
#ifdef BOINC
        boinc_finish(EXIT_FAILURE);
#endif
        exit(EXIT_FAILURE);
    }

    ftime(&starttime);

    time_t timer;
    char curdatetime[26];
    struct tm* tm_info;
    time(&timer);
    tm_info = localtime(&timer);
    strftime(curdatetime, 26, "%d.%m.%Y %H:%M:%S", tm_info);

#ifndef BOINC
    sprintf(repfname, "rep_%06" PRIu32 ".%06" PRIu32 "-%06" PRIu32 ".%06" PRIu32, minc, minb, maxc, maxb);
    sprintf(outfname, "out_%06" PRIu32 ".%06" PRIu32 "-%06" PRIu32 ".%06" PRIu32, minc, minb, maxc, maxb);
    sprintf(chkfname, "chk_%06" PRIu32 ".%06" PRIu32 "-%06" PRIu32 ".%06" PRIu32, minc, minb, maxc, maxb);
#endif

    int error_code, checkpoint_code;
    error_code = checkpoint_code = read_checkpoint();
    if (error_code) error_code = init_task();
    if (error_code) return error_code;

    fout = fopen(outfname, "r");
    if (skip && fout != NULL && checkpoint_code) {
        fclose(fout);
#ifdef BOINC
	boinc_finish(EXIT_SUCCESS);
#endif
        exit (EXIT_SUCCESS);
    }
    if (output) {
        if (!checkpoint_code && fout == NULL) {
#ifdef BOINC
            boinc_finish(EXIT_FAILURE);
#endif
            exit(EXIT_FAILURE);
        }
        if (checkpoint_code) {
            fout = fopen(outfname, "w");
        } else {
            fout = fopen(outfname, "a");
        }
        if (fout == NULL) {
#ifdef BOINC
            boinc_finish(EXIT_FAILURE);
#endif
            exit(EXIT_FAILURE);
        }
    }

    fprintf(stderr, "Command line    :");
    for (int i = 1; i < argc; i++)
        fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "%s from   : %" PRIu32" %" PRIu32"\n", checkpoint_code ? "Starting" : "Resuming", curc, curb);
    fprintf(stderr, "Local timestamp : %s\n", curdatetime);
#ifdef BOINC
    fprintf(stderr, "\n");
#endif

    init_primes((uint64_t)maxc * (uint64_t)maxc * 3);

    int cpcnt, ctpcnt = 0;
    float cstep = 0.001;
    int fpcnt, ftpcnt = 0;
    float fstep = 0.0001;

    if (progress)
        do_progress(ctpcnt);
#ifdef BOINC
    boinc_fraction_done(0);
#endif

    uint32_t a, b, c, gbc, bmin, total = (maxc - minc) * (maxc + 1) / 2 + maxb - minb, step;
    uint64_t aa, bb, cc;
    __uint128_t area;
    for (c = curc; c <= maxc; c++) {
        cc = (uint64_t)c * (uint64_t)c;
        bmin = ceil(((float)c+1)/2);
        for (b = max(max(bmin, curb), ((c == minc) ? minb : 0)); b < ((c == maxc) ? maxb : c+1); b++) {
            curb = 0;
            if (!(c&1) && !(b&1)) continue;
            bb = (uint64_t)b * (uint64_t)b;
            gbc = gcd(b, c);
            counter = 0;
            for (a = lrintl(sqrtl(cc - bb)) + 1; a < b+1; a++) {
                if ((a & 1) + (b & 1) + (c & 1) != 2 || gcd(a, gbc) != 1) continue;
                aa = (uint64_t)a * (uint64_t)a;
                clear_factors();
                factorize(aa + bb - cc);
                factorize(aa - bb + cc);
                factorize(-aa + bb + cc);
                factorize(aa + bb + cc);
                area = is_square() >> 2;
                if (area) {
                    char s128[40];
                    u128_to_string(area, s128);
                    if (!quiet) {
                        fprintf(stdout, "\na=%" PRIu32 "^2, b=%" PRIu32 "^2, c=%" PRIu32 "^2, area=%s\n", a, b, c, s128);
                        fflush(stdout);
                    }
                    if (output) {
                        fprintf(fout, "a=%" PRIu32 "^2, b=%" PRIu32 "^2, c=%" PRIu32 "^2, area=%s\n", a, b, c, s128);
                        fflush(fout);
                    }
                    found += 1;
                }
            }
            checksum += counter;
            step = (c - minc) * (maxc + 1) / 2 + b - bmin;
            fpcnt = (int)((double)step / total / fstep);
            if (ftpcnt != fpcnt) {
                ftpcnt = fpcnt;
#ifdef BOINC
                boinc_fraction_done((double)step / total);
#endif
            }
            cpcnt = (int)((double)step / total / cstep);
            if (ctpcnt != cpcnt) {
                ctpcnt = cpcnt;
                if (progress)
                    do_progress((double)step / total * 100);
                save_checkpoint(c, b);
                fflush(stdout);
            }
        }
    }

    do_progress((double)100);

    if (output) fclose(fout);
    remove(chkfname);

    ftime(&endtime);
    uint64_t dif = (endtime.time - starttime.time) * 1000 + (endtime.millitm - starttime.millitm);

#ifndef BOINC
    fprintf(stderr, "\n");
#endif
    fprintf(stderr, "Elapsed time    : %02d:%02d:%02d.%03d\n", (unsigned char)(dif/60/60/1000), (unsigned char)((dif/60/1000)%60), (unsigned char)((dif/1000)%60), (unsigned char)(dif%1000));
#ifndef BOINC
    fprintf(stderr, "Check sum       : %" PRIu64 "\n", checksum);
#endif

    if (report) {
        frep = fopen(repfname, "w");
        if(frep == NULL) {
            perror("Failed to open rep file");
#ifdef BOINC
			boinc_finish(EXIT_FAILURE);
#endif
            exit(EXIT_FAILURE);
        }
        fprintf(frep,  "%s,%s,%s,%s"
                      ",%" PRIu64
                      ",%" PRIu8
#ifndef BOINC
                      ",%s,%s,%02d:%02d:%02d.%03d"
#endif
                      "\n"
                    ,argv[1]
                    ,argv[2]
                    ,argv[3]
                    ,argv[4]
                    ,checksum
                    ,found
#ifndef BOINC
                    ,VERSION
                    ,OS
                    ,(unsigned char)(dif/60/60/1000), (unsigned char)((dif/60/1000)%60), (unsigned char)((dif/1000)%60), (unsigned char)(dif%1000)
#endif
               );
        fclose(frep);
    }
    free_primes();
#ifdef BOINC
    boinc_finish(EXIT_SUCCESS);
#endif

    return (EXIT_SUCCESS);
}
