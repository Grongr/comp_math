import matplotlib.pyplot as plt
import numpy as np  
import random
import math

MINH = 1e-5
MAXH = 1
points = [-2, -1, 0, 1, 2, 3]

def arrToStr(arr):
    s = ''
    for i in arr:
        s += str(i)+' '
    return s[:-1]

# def generateData(n, dots):
#     f = open('data.txt','w')
#     f.write(f'{n} {dots}\n')
#     arr = []
#     for i in range(1, n+1):
#         arr.append([j*(i*(MAXH - MINH)/n + MINH) for j in points])
#         f.writ
#         e(f'{arrToStr(arr[i-1])}\n')
#     print("DONE!")
#     f.close()

def aim(x):
    return ((math.log(x) + 1) * (x**x))

def genX(h):
    x = []
    for hi in h:
        x.append([])
        #print(len(points))
        for dot in points:
            x[-1].append(dot*hi + math.pi)
            #print(i*hi + math.pi, end=' ')
        #print()
    return x

def calculateD(x, c, h):
    ans = 0
    for i in range(len(x)):
        ans += x[i]**x[i]*c[i]
    return ans/h

def main():
    h    = [2**i for i in range(-20, 1)]
    x    = genX(h)
    coef = [0.05, -0.5, -0.333333, 1, -0.25, 0.0333333]
    err  = [abs(aim(math.pi) - calculateD(x[i], coef, h[i])) for i in range(len(h))]
    
    h = [np.log(h[i]) / np.log(2) for i in range(len(h))]
    err = [np.log(err[i]) / np.log(2) for i in range(len(err))]
    
    fig, ax = plt.subplots()
    plt.plot(h, err, label = 'err(h)', color='red')
    
    ax.legend()
    plt.show()
main()
