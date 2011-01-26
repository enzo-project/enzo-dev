import numpy as na
import matplotlib.pyplot as plt
import os

# Assume enzo has been run with "enzo.exe TestOrbit > TestOrbit.out"

FileIn = open("TestOrbit.out", "r")
xpos = []
ypos = []
for line in FileIn:
   if 'id=1' in line:
      lst = line.split()
      xpos.append(lst[1])
      ypos.append(lst[2])
FileIn.close()
x = na.array(xpos)
y = na.array(ypos)

plt.plot(x[::2], y[::2], marker='o')
plt.axis("equal")
plt.savefig('TestOrbit.png')
