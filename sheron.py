import argparse
import multiprocessing
from math import gcd, sqrt, ceil
from sympy.ntheory.primetest import is_square
from datetime import datetime
from timeit import default_timer as timer
from sys import maxsize
from time import sleep


def get_next():
    global n, squares, finish, exit_flag
    if not exit_flag:
        r = n
        n += 1
        squares.append(n*n)
        exit_flag = r == finish
    else:
        return None
    return r


def run(c, sq):
    run_time = timer()
    if c is None:
        return
    cc = sq[c]
    for b in range(ceil((c + 1) / 2), c + 1):
        if not c & 1 and not b & 1:
            continue
        bb = sq[b]
        gbc = gcd(b, c)
        for a in range(int(sqrt(cc - bb)) + 1, b + 1):
            if (a & 1) + (b & 1) + (c & 1) != 2 or gcd(a, gbc) != 1:
                continue
            aa = sq[a]
            area2 = (aa + bb + cc) * (aa + bb - cc) * (aa - bb + cc) * (- aa + bb + cc)
            if is_square(area2):
                area = int(sqrt(area2)) // 4
                print(f'\r{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} - (a={a}², b={b}², c={c}², area={area}) - '
                      f'Heron Triangle with Square Sides!')
                with open('sheron.out', 'a') as f:
                    f.write(f'\ra={a}^2, b={b}^2, c={c}^2, area={area}')
    run_duration = datetime.utcfromtimestamp(timer() - run_time)
    print(f'\r{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} - c={c}² - '
          f'Cycle: {run_duration.strftime("%H:%M:%S.%f")[:-3]}', end='')


if __name__ == '__main__':
    start_time = timer()

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--start', type=int, help='lower bound of interval', default=1)
    parser.add_argument('-f', '--finish', type=int, help='upper bound of interval', default=maxsize)
    parser.add_argument('-mt', '--threads', type=int, help='multi threads count to use', default=1)
    args = parser.parse_args()
    n = args.start
    finish = args.finish
    thread_count = args.threads

    try:
        with open('sheron.chk', 'r') as f:
            n, dif = map(int, f.readline().split(","))
            n += 1
        start_time = timer() - dif / 1000
    except:
        pass

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} - Resuming from c={n}² with {thread_count} thread(s)')

    squares = [i*i for i in range(n+1)]
    threads = [None for _ in range(thread_count)]
    threadLock = multiprocessing.Lock()
    exit_flag = False

    while not exit_flag:
        for i in range(len(threads)):
            if threads[i] and not threads[i].is_alive() and min([t.name for t in threads]) == threads[i].name:
                dif = int(round((timer() - start_time) * 1000))
                threadLock.acquire()
                with open('sheron.chk', 'w') as f:
                     f.write(",".join(map(str, (threads[i].name, dif))))
                threadLock.release()
            if not threads[i] or not threads[i].is_alive():
                m = get_next()
                thread = multiprocessing.Process(target=run, args=(m, squares))
                thread.name = str(m)
                thread.start()
                threads[i] = thread
        sleep(.2)
