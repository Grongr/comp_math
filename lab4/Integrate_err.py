import matplotlib.pyplot as plt
import numpy as np  
import math

f = open("Integrate_err.txt", 'r')
data = f.readlines()

err_ = list(map(float, data[0].split()))
num_ = [i for i in range(1, len(err_)+1)]

num = []
err = []

for i in range(len(num_)):
    if (err_[i] != 0):
        num.append(num_[i])
        err.append(err_[i])

num = [np.log(num[i]) / np.log(2) for i in range(len(num))]
err = [np.log(err[i]) / np.log(2) for i in range(len(err))]

fig, ax = plt.subplots()

plt.plot(num, err)

#plt.xlim(0, 5)
#plt.ylim(-5, 0)
#ax.legend()
plt.show()