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
z = args.finish
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

while m <= z:
    print_current_state(datetime.utcfromtimestamp(timer() - start_loop))
    start_loop = timer()
    for n in range(1, m):
        if gcd(m, n) > 1:
            continue
        mn = m*n
        mmnn = m*m + n*n
        for p in range(2, m+1):
            for q in range(1, p):
                if gcd(p, q) > 1:
                    continue
                mnpp = mn * p * p
                mmnnp = mmnn * p
                a = mmnnp * q                     # a = p * q * (m*m + n*n)
                b = mnpp + mn*q*q                 # b = m * n * (p*p + q*q)
                c = (m*q + n*p) * (m*p - n*q)     # c = (m*q + n*p) * (m*p - n*q)
                triangles += 1
                # There is no isosceles Heron triangles with square sides
                # Stănică, Pantelimon; Sarkar, Santanu; Sen Gupta, Sourav; Maitra, Subhamoy; Kar, Nirupam (2013).
                # "Counting Heron triangles with constraints". Integers. 13: Paper No. A3, 17pp., p. 10
                if a == b or b == c or c == a:
                    continue
                proportion = gcd(a, b, c)
                a, b, c = map(lambda w: w // proportion, sorted([a, b, c]))
                # Check Heronian triangle (a, b, c) for square sides
                if is_square(a) and is_square(b):
                    almost += 1
                    roota = isqrt(a)
                    rootb = isqrt(b)
                    print(f"{current_time()} - almost:{almost} ({roota}², {rootb}², {c}) <- ({m=}, {n=}, {p=}, {q=})")
                    if is_square(c):
                        # Heron's formula the area of a triangle whose sides have lengths a, b, and c
                        # s is semi-perimeter
                        s = (a + b + c) // 2
                        area = isqrt(s * (s-a) * (s-b) * (s-c))
                        htss += 1
                        print(f"{current_time()} - HTSS:{htss} ({a=}, {b=}, {c=}, {area=}) <- ({m=}, {n=}, {p=}, {q=})")
                        with open(out_file, 'a') as f:
                            f.write(f"{m=},{n=},{p=},{q=},{a=},{b=},{c=},{area=}")
    m += 1
    dif = int(round((timer() - start_time) * 1000))
    with open(chk_file, 'w') as f:
        f.write(",".join(map(str, (m, triangles, almost, htss, dif))))
