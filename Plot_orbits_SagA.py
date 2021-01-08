import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy import optimize

# Scripts to plot the orbits of the stars around Sagitarious A from their orbital parameters
# from Eisenhauer et al. 2005 ApJ, 628, 246

x0 = 0.01

a = [0.412,0.1226,0.329,0.286,0.219,0.225]
e = [0.358,0.8760,0.927,0.902,0.395,0.9389]
w = [129.8 * 0.0174533,62.6 * 0.0174533,159.2 * 0.0174533,311.8 * 0.0174533,250 * 0.0174533,344.7 * 0.0174533]
i = [120.5 * 0.0174533,131.9 * 0.0174533,60.6 * 0.0174533,32.8 * 0.0174533,11 * 0.0174533,97.3 * 0.0174533]
omega = [341.5 * 0.0174533,221.9 * 0.0174533,141.4 * 0.0174533,233.3 * 0.0174533,100 * 0.0174533,228.5 * 0.0174533]
p = [94.1,15.24,67.2,54.4,36,38]
t0 = [2002.6,2002.315,1987.71,1995.628,2006.1,2000.156]

readyx1,readyx2,readyx8,readyx12,readyx13,readyx14 = [],[],[],[],[],[]
readyx = [readyx1,readyx2,readyx8,readyx12,readyx13,readyx14]
readyy1,readyy2,readyy8,readyy12,readyy13,readyy14 = [],[],[],[],[],[]
readyy = [readyy1,readyy2,readyy8,readyy12,readyy13,readyy14]

def f(Et, e, M):
    return Et - e * np.sin(Et) - M

for q in range(len(t0)):
    rotx, roty = [],[]
    for t in arange(t0[q], t0[q] + p[q], 0.01):
        M = (360 - 360 * (t0[q] - t) / p[q]) % 360
        M = M * 2 * np.pi / 180

        zero = optimize.newton(f, x0, args=[e[q], M])
        nu = 2.0 * np.arctan2(np.sqrt(1.0 + e[q]) * np.sin((zero) / 2.0), np.sqrt(1.0 - e[q]) * np.cos((zero) / 2.0))
        r = a[q] * (1.0 - e[q] * np.cos(zero))
        ox = r * np.cos(nu)
        oy = r * np.sin(nu)
        rx = (ox * (np.cos(w[q]) * np.cos(omega[q]) - np.sin(w[q]) * np.cos(i[q]) * np.sin(omega[q])) - oy * (
                np.sin(w[q]) * np.cos(omega[q]) + np.cos(w[q]) * np.cos(i[q]) * np.sin(omega[q])))

        ry = (ox * (np.cos(w[q]) * np.sin(omega[q]) + np.sin(w[q]) * np.cos(i[q]) * np.cos(omega[q])) + oy * (
                np.cos(w[q]) * np.cos(i[q]) * np.cos(omega[q]) - np.sin(w[q]) * np.sin(omega[q])))

        readyx[q].append(rx * np.cos(np.pi / 2.0) - ry * np.sin(np.pi / 2.0))
        readyy[q].append(rx * np.sin(np.pi / 2.0) + ry * np.cos(np.pi / 2.0))

plt.figure(figsize=(8, 9))
plt.xlim(0.5, -0.3)
plt.ylim(-0.45, 0.55)
circle = plt.Circle((0, 0), 0.006, color='black')
plt.gcf().gca().add_artist(circle)
plt.title('Orbits of S stars around Sagitarius A')
plt.xlabel('Offest in RA (arcsec)')
plt.ylabel('Offest in DEC (arcsec)')
plt.grid(color='grey', linestyle='-', linewidth=0.5)
plt.plot(-np.array(readyx1), np.array(readyy1), '-', linewidth=1.5, color='darkblue', label='S1')
plt.plot(-np.array(readyx2), np.array(readyy2), '-', linewidth=1.5, color='red', label='S2')
plt.plot(-np.array(readyx8), np.array(readyy8), '-', linewidth=1.5, color='darkcyan', label='S8')
plt.plot(-np.array(readyx12), np.array(readyy12), '-', linewidth=1.5, color='black', label='S12')
plt.plot(-np.array(readyx13), np.array(readyy13), '-', linewidth=1.5, color='magenta', label='S13')
plt.plot(-np.array(readyx14), np.array(readyy14), '-', linewidth=1.5, color='darkgreen', label='S14')
plt.legend()
plt.savefig('orbits.png')
plt.show()