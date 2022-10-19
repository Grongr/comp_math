import matplotlib.pyplot as plt
import numpy as np  
import math

def genData(n, a = 0, b = 3):
    l = (b - a) / (n+1)
    print(n)
    for i in range(n):
        print(a + l*(i+1), end=' ')
    print()
    
def aim(x):
    return math.cos(x)
    
def printData():
    num = [2**i for i in range(1, 8)]
    for n in num:
        genData(n)
    
def main():
    
    fig, ax = plt.subplots()
    
    num = [2**i for i in range(1, 8)]
    err = [0.459698, 0.256359, 0.0971266, 0.0292395, 0.00796214, 0.00207277, 0.000528458]
    
    plt.plot(num, err)
    
    ax.legend()
    plt.show()

# printData()
main()