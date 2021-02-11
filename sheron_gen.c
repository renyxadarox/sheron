#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

#define min(a,b) ((a) < (b) ? a : b)

void generate_task(uint32_t k, uint32_t l, uint32_t m, uint32_t n)
{
    FILE * fout;
    char fname[256];
    sprintf(fname, "sheron_%06" PRIu32 ".%06" PRIu32 "-%06" PRIu32 ".%06" PRIu32, k, l, m, n);
    fout = fopen(fname, "w");
    fprintf(fout, "%" PRIu32 " %" PRIu32 " %" PRIu32 " %" PRIu32, k, l, m, n);
    fclose(fout);
}

int main(int argc, char** argv)
{
    char *p;
    uint32_t c = strtoul(argv[1], &p, 10), d, s, minc, maxc, minb, maxb;
    uint32_t n = strtoul(argv[2], &p, 10);
    long double l;
    for (uint32_t i = 0; i < n;) {
        l = (0.00000151809654627994L*c*c-0.01733792620682095276L*c+120.38447044848726363853L)/3600L;
        minc = c;
        if (l < 1) {
            d = min((int)(1/l), 1000);
            maxc = c + d;
            minb = maxb = 0;
            generate_task(minc, minb, maxc, maxb);
            c += d;
            i++;
        } else {
            d = (int)l+1;
            maxc = c;
            minb = 0;
            maxb = s = ceil(sqrtl((long double)c * (long double)c / 2));
            for (uint32_t j = 0; j < d; j++) {
                if (j) minb = maxb;
                if (j == d-1) {maxc = c + 1; maxb = 0;}
                else maxb = s + (c - s) * sqrtl((long double)(j + 1)/d);
                generate_task(minc, minb, maxc, maxb);
                i++;
            }
            c++;
        }
    }
    fprintf(stdout, "%" PRIu32, c);
    return 0;
}
