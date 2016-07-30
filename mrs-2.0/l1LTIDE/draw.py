import numpy as np
import matplotlib.pyplot as plt
import fileinput

line = plt.figure()


points = []
hyper_planes = []
read_points = True
for line in fileinput.input():
    print(line)
    if line == "\n":
        read_points = False
        continue
    if read_points == True:
        entries = line.split(" ")
        x = np.float(entries[1])
        y = np.float(entries[2])
        points.append((x,y))
    else:
        hyp = line.split(" ")
        x = np.float(hyp[1])
        y = np.float(hyp[2])
        c = np.float(hyp[3])
        hyper_planes.append((x,y,c))

plt.plot(*zip(*points), "o")

hyp_points = []
x_range = np.arange(-1, 1, 0.001)
for h in hyper_planes:
    x = h[0] / h[1]
    c = h[2] / h[1]
    y = x * x_range + c
    hyp_points.append(y)

print(len(hyper_planes))
for pnts in hyp_points:
    plt.plot(x_range, pnts)
#    print(pnts)
#for point in points:
#    print(point)

#np.random.seed()
#x = np.random.uniform(-1, 1, 500)
#y = np.random.uniform(-1, 1, 500)
axes = plt.gca()
axes.set_xlim([-1, 1])
axes.set_ylim([-1, 1])

plt.show()
