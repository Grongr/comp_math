import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


def plot(df):
    plt.figure(figsize=(5, 3))
    plt.plot(df['B'], df['A'], '-', color='lightblue')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid()
    plt.show()


files = ['out.csv']
for file in files:
    data = pd.read_csv(file, delimiter=',')
    plot(data)
