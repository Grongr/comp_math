import matplotlib.pyplot as plt
import numpy as np  
import math

def lin_ls(x, y, through_null=False):
    """
    Eng:
        Returns coefficients b, k and their deviations s_k, s_b
    from linear approximation (using least squares) y = k * x + b to
    given points (x, y).
        If the flag through_null is True, then approximation does not include
    coefficient b (as well as s_b), b is considered to be zero.
    Rus:
        Возвращает коэффициенты b, k и их погрешности s_b, s_k при
    линейной аппроксимации y = k * x + b методом МНК по введённым точкам (х, у).
        Если флаг through_null равен True, то аппроксимация проводится через
    точку (0,0) и b, s_b не включены в вывод функции.

    Parameters
    ----------
    x : numpy.ndarray
        Numpy array with x coordinates of points.

    y : numpy.ndarray
        Numpy array with y coordinates of points. Vectors x and y must be the same length.

    through_null: bool, optional
        Flag for making approximation through (0, 0) point.

    Returns
    -------
    out : tuple
        If through_null is True: (k, s_k),
        If through_null is False: (k, s_k, b, s_b)
    """
    if isinstance(x, np.ndarray) and isinstance(y, np.ndarray):
        if len(x) != len(y):
            raise ValueError("Incompatible x and y vectors. They must have the same length.")
        if through_null:
            k = np.mean(x * y) / np.mean(x * x)
            s_k = np.sqrt(1 / len(x)) * np.sqrt(np.mean(y * y) / np.mean(x * x) - k ** 2)
            return k, s_k
        else:
            xy = np.mean(x * y)
            x1y = np.mean(x) * np.mean(y)
            x2 = np.mean(x * x)
            x12 = np.mean(x) ** 2
            y2 = np.mean(y * y)
            y12 = np.mean(y) ** 2
            k = (xy - x1y) / (x2 - x12)
            b = np.mean(y) - k * np.mean(x)
            s_k = np.sqrt(1 / len(x)) * np.sqrt((y2 - y12) / (x2 - x12) - k ** 2)
            s_b = s_k * np.sqrt(x2 - x12)
            return k, s_k, b, s_b
    else:
        raise ValueError("Invalid x or/and y type. Must be numpy.ndarray.")

def genData(n, a = 0, b = 3):
    l = (b - a) / (n-1)
    print(n)
    for i in range(n):
        print(a + l*(i), end=' ')
    print()
    
def aim(x):
    return math.cos(x)
    
def printData():
    num = [2**i for i in range(1, 8)]
    for n in num:
        genData(n)
    
# n = int(input())
# genData(n)
    
fig, ax = plt.subplots()
    
num = [2**i for i in range(1, 8)]
##err = [0.459698, 0.256359, 0.0971266, 0.0292395, 0.00796214, 0.00207277, 0.000528458]
err = [0.229413, 0.057662, 0.0013495, 0.00037062, 0.00028901, 0.00001333, 3.97e-6]

num = [np.log(num[i]) / np.log(2) for i in range(len(num))]
err = [np.log(err[i]) / np.log(2) for i in range(len(err))]

k, sk, b, sb = lin_ls(np.array(num[2:]), np.array(err[2:]))
x = [1, 7]
y = [k*i + b + 0.2 for i in x]

plt.plot(x, y, color = 'red', label = f'y = {"%.3f" % k}*x + {"%.3f" % b} + 0.2')
plt.plot(num, err)

ax.legend()
ax.grid()
plt.show()

#printData()
