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
    
    # fig, ax = plt.subplots(2, 2)
    fig, ax = plt.subplots()
    
    num_3 = [1, 2, 3, 4, 5, 6, 7]
    
    err_3  = [1.64045, 0.405861, 0.0552836, 0.00429133, 0.000282585, 1.78916e-05, 1.12184e-06]
    err_5  = [0.528715, 0.0258387, 0.000827636, 1.57711e-05, 2.58476e-07, 4.08675e-09, 6.40442e-11]
    err_10 = [0.00235608, 1.8643e-06, 2.07866e-09, 3.02747e-12, 3.44169e-15, 2.22045e-16, 3.33067e-16]


    # ax[0].plot(num_3, err_3,  label="For  3 points")
    # ax[1].plot(num_3, err_3,  label="For  5 points")
    # ax[2].plot(num_3, err_10, label="For 10 points")
    
    # ax[3].plot(num_3, err_3,  label="For  3 points")
    # ax[3].plot(num_3, err_5,  label="For  5 points")
    # ax[3].plot(num_3, err_10, label="For 10 points")
    
    ax.plot(num_3, err_3,  label="For  3 points")
    ax.plot(num_3, err_5,  label="For  5 points")
    ax.plot(num_3, err_10, label="For 10 points")

    ax.set_ylabel("Errors")
    ax.set_xlabel("For ranges: 8, 4, 2, 1, 0.5, 0.25, 0.125")
    
    ax.legend()
    plt.show()

main()
