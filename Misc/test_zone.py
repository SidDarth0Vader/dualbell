import matplotlib.pyplot as mat
import numpy as np
import time


i = 0
j = 0
while i < 4:
    j = 0
    while j < 4:
        if i != 1 or j != 3:
            print(i, j)
        j+=1
    i+=1
    #print('umm')