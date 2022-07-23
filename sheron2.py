import argparse
from math import gcd, isqrt
from sympy.ntheory.primetest import is_square
from datetime import datetime
from sys import maxsize
from timeit import default_timer as timer

# The application searching for a Heronian Triangles with Square Sides

start_time = timer()

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-s', '--start', type=int, help='lower bound', default=2)
parser.add_argument('-f', '--finish', type=int, help='upper bound', default=maxsize)
args = parser.parse_args()
m = args.start
finish = args.finish
dif = 0
chk_file = 'sheron2.chk'
out_file = 'sheron2.out'

triangles = 0   # Heronian triangles
almost = 0      # almost Heronian Triangles with Square Sides
htss = 0        # Heronian Triangles with Square Sides (HTSS)


def current_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def print_current_state(dur):
    print(f'{current_time()} - {m=}, {triangles=}, {almost=}, {htss=}, duration={dur.strftime("%H:%M:%S.%f")[:-3]}')


try:
    with open(chk_file, 'r') as f:
        m, triangles, almost, htss, dif = map(int, f.readline().split(","))
    print(f'{current_time()} - Resuming from')
    start_time = timer() - dif / 1000
except:
    print(f'{current_time()} - Starting from')

start_loop = start_time

# Leonard Euler method of generating Heronian triangles
# a = mn(p²+q²)
# b = pq(m²+n²)
# c = (mq+np)(mp-nq)
# for integers m, n, p and q such that (m, n)=1, (p, q)=1

while m <= finish:
    print_current_state(datetime.utcfromtimestamp(timer() - start_loop))
    start_loop = timer()
    for n in range(1, m):
        if gcd(m, n) > 1:
            continue
        mn = m*n
        mmnn = m*m + n*n
        for p in range(2, m+1):
            mnpp = mn * p * p
            mmnnp = mmnn * p
            mp = m * p
            np = n * p
            for q in range(1, p):
                if gcd(p, q) > 1:
                    continue
                a = mmnnp * q
                b = mnpp + mn * q * q
                for i in (0, 1):
                    c = (m*q + np) * (mp - n*q) if i == 0 else (n*q + mp) * (np - m*q)
                    if c < 0:
                        continue
                    triangles += 1
                    # There is no isosceles Heron triangles with square sides
                    # Stănică, Pantelimon; Sarkar, Santanu; Sen Gupta, Sourav; Maitra, Subhamoy; Kar, Nirupam (2013).
                    # "Counting Heron triangles with constraints". Integers. 13: Paper No. A3, 17pp., p. 10
                    if a == b or b == c or c == a:
                        continue
                    if i:
                        m, n = n, m
                    proportion = gcd(a, b, c)
                    x, y, z = map(lambda w: w // proportion, sorted([a, b, c]))
                    # Check Heronian triangle (x, y, z) for square sides
                    if is_square(x) and is_square(y):
                        rootx = isqrt(x)
                        rooty = isqrt(y)
                        if is_square(z):
                            rootz = isqrt(z)
                            # Heron's formula the area of a triangle whose sides have lengths x, y, and z
                            # s is semi-perimeter
                            s = (x + y + z) // 2
                            area = isqrt(s * (s-x) * (s-y) * (s-z))
                            htss += 1
                            print(f"{current_time()} - HTSS:{htss} ({rootx}², {rooty}², {rootz}², {area=}) <- ({m=}, {n=}, {p=}, {q=})")
                            with open(out_file, 'a') as f:
                                f.write(f"{m=},{n=},{p=},{q=},{x=},{y=},{z=},{area=}\n")
                        else:
                            almost += 1
                            print(f"{current_time()} - almost:{almost} ({rootx}², {rooty}², {z}) <- ({m=}, {n=}, {p=}, {q=})")
                    if i:
                        m, n = n, m
    m += 1
    dif = int(round((timer() - start_time) * 1000))
    with open(chk_file, 'w') as f:
        f.write(",".join(map(str, (m, triangles, almost, htss, dif))))
