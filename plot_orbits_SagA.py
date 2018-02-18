import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import scipy.optimize
from scipy import optimize
from scipy import ndimage

x0 = 0.01

a1     = 0.412;     a2 = 0.1226;    a8 = 0.329;     a12 = 0.286;     a13 = 0.219;   a14 = 0.225
e1     = 0.358;     e2 = 0.8760;    e8 = 0.927;     e12 = 0.902;     e13 = 0.395;   e14 = 0.9389
w1     = 129.8*0.0174533;     w2 = 62.6*0.0174533;      w8 = 159.2*0.0174533;     w12 = 311.8*0.0174533;     w13 = 250*0.0174533;     w14 = 344.7*0.0174533
i1     = 120.5*0.0174533;     i2 = 131.9*0.0174533;     i8 = 60.6*0.0174533;      i12 = 32.8*0.0174533;      i13 = 11*0.0174533;      i14 = 97.3*0.0174533 
omega1 = 341.5*0.0174533; omega2 = 221.9*0.0174533; omega8 = 141.4*0.0174533; omega12 = 233.3*0.0174533; omega13 = 100*0.0174533; omega14 = 228.5*0.0174533;
p1     = 94.1;      p2 = 15.24;     p8 = 67.2;      p12 = 54.4;      p13 = 36;      p14 = 38 
t01    = 2002.6;   t02 = 2002.315; t08 = 1987.71;  t012 = 1995.628; t013 = 2006.1; t014 = 2000.156

rotx1=[]; rotx2=[]; rotx8=[]; rotx12=[]; rotx13=[]; rotx14=[]
roty1=[]; roty2=[]; roty8=[]; roty12=[]; roty13=[]; roty14=[]

for t in arange (t01,t01+p1,0.01):
	M1 = (360 - 360*(t01-t)/p1)%360
	M1 = M1*2*np.pi/180
	def f1(Et1):
		return Et1 - e1*np.sin(Et1) - M1
	shit1 = optimize.newton(f1, x0)
	nu1 = 2.0*np.arctan2(np.sqrt(1.0+e1)*np.sin((shit1)/2.0),np.sqrt(1.0-e1)*np.cos((shit1)/2.0))
	r1 = a1*(1.0-e1*np.cos(shit1))
	ox1 = r1*np.cos(nu1)
	oy1 = r1*np.sin(nu1)
	rx1=(ox1*(np.cos(w1)*np.cos(omega1)-np.sin(w1)*np.cos(i1)*np.sin(omega1))- oy1*(np.sin(w1)*np.cos(omega1)
	+ np.cos(w1)*np.cos(i1)*np.sin(omega1)))
	
	ry1=(ox1*(np.cos(w1)*np.sin(omega1)+np.sin(w1)*np.cos(i1)*np.cos(omega1))+ oy1*(np.cos(w1)*np.cos(i1)
	*np.cos(omega1) - np.sin(w1)*np.sin(omega1)))
	rotx1.append(rx1*np.cos(np.pi/2.0)-ry1*np.sin(np.pi/2.0))
	roty1.append(rx1*np.sin(np.pi/2.0)+ry1*np.cos(np.pi/2.0))
readyx1 = np.array(rotx1)
readyy1 = np.array(roty1)

for t in arange (t02,t02+p2,0.01):
	M2 = (360.0 - 360.0*(t02-t)/p2)%360.0
	M2 = M2*2.0*np.pi/180.0
	def f2(Et2):
		return Et2 - e2*np.sin(Et2) - M2
	shit2 = optimize.newton(f2, x0)
	nu2 = 2.0*np.arctan2(np.sqrt(1.0+e2)*np.sin((shit2)/2.0),np.sqrt(1.0-e2)*np.cos((shit2)/2.0))
	r2 = a2*(1.0-e2*np.cos(shit2))
	ox2 = r2*np.cos(nu2)
	oy2 = r2*np.sin(nu2)
	rx2=(ox2*(np.cos(w2)*np.cos(omega2)-np.sin(w2)*np.cos(i2)*np.sin(omega2))- oy2*(np.sin(w2)*np.cos(omega2)
	+ np.cos(w2)*np.cos(i2)*np.sin(omega2)))
	
	ry2=(ox2*(np.cos(w2)*np.sin(omega2)+np.sin(w2)*np.cos(i2)*np.cos(omega2))+ oy2*(np.cos(w2)*np.cos(i2)
	*np.cos(omega2) - np.sin(w2)*np.sin(omega2)))
	rotx2.append(rx2*np.cos(np.pi/2.0)-ry2*np.sin(np.pi/2.0))
	roty2.append(rx2*np.sin(np.pi/2.0)+ry2*np.cos(np.pi/2.0))
readyx2 = np.array(rotx2)
readyy2 = np.array(roty2)

for t in arange (t08,t08+p8,0.01):
	M8 = (360.0 - 360.0*(t08-t)/p8)%360.0
	M8 = M8*2.0*np.pi/180.0
	def f8(Et8):
		return Et8 - e8*np.sin(Et8) - M8
	shit8 = optimize.newton(f8, x0)
	nu8 = 2.0*np.arctan2(np.sqrt(1.0+e8)*np.sin((shit8)/2.0),np.sqrt(1.0-e8)*np.cos((shit8)/2.0))
	r8 = a8*(1.0-e8*np.cos(shit8))
	ox8 = r8*np.cos(nu8)
	oy8 = r8*np.sin(nu8)
	rx8=(ox8*(np.cos(w8)*np.cos(omega8)-np.sin(w8)*np.cos(i8)*np.sin(omega8))- oy8*(np.sin(w8)*np.cos(omega8)
	+ np.cos(w8)*np.cos(i8)*np.sin(omega8)))
	
	ry8=(ox8*(np.cos(w8)*np.sin(omega8)+np.sin(w8)*np.cos(i8)*np.cos(omega8))+ oy8*(np.cos(w8)*np.cos(i8)
	*np.cos(omega8) - np.sin(w8)*np.sin(omega8)))
	rotx8.append(rx8*np.cos(np.pi/2.0)-ry8*np.sin(np.pi/2.0))
	roty8.append(rx8*np.sin(np.pi/2.0)+ry8*np.cos(np.pi/2.0))
readyx8 = np.array(rotx8)
readyy8 = np.array(roty8)
	
for t in arange (t012,t012+p12,0.01):
	M12 = (360.0 - 360.0*(t012-t)/p12)%360.0
	M12 = M12*2.0*np.pi/180.0
	
	def f12(Et12):
		return Et12 - e12*np.sin(Et12) - M12
	shit12 = optimize.newton(f12, x0)
	nu12 = 2.0*np.arctan2(np.sqrt(1.0+e12)*np.sin((shit12)/2.0),np.sqrt(1.0-e12)*np.cos((shit12)/2.0))
	r12 = a12*(1.0-e12*np.cos(shit12))
	ox12 = r12*np.cos(nu12)
	oy12 = r12*np.sin(nu12)
	rx12=(ox12*(np.cos(w12)*np.cos(omega12)-np.sin(w12)*np.cos(i12)*np.sin(omega12))- oy12*(np.sin(w12)*np.cos(omega12)
	+ np.cos(w12)*np.cos(i12)*np.sin(omega12)))
	
	ry12=(ox12*(np.cos(w12)*np.sin(omega12)+np.sin(w12)*np.cos(i12)*np.cos(omega12))+ oy12*(np.cos(w12)*np.cos(i12)
	*np.cos(omega12) - np.sin(w12)*np.sin(omega12)))
	rotx12.append(rx12*np.cos(np.pi/2.0)-ry12*np.sin(np.pi/2.0))
	roty12.append(rx12*np.sin(np.pi/2.0)+ry12*np.cos(np.pi/2.0))
readyx12 = np.array(rotx12)
readyy12 = np.array(roty12)

for t in arange (t013,t013+p13,0.001):
	M13 = (360.0 - 360.0*(t013-t)/p13)%360.0
	M13 = M13*2.0*np.pi/180.0
	def f13(Et13):
		return Et13 - e13*np.sin(Et13) - M13
	shit13 = optimize.newton(f13, x0)
	nu13 = 2.0*np.arctan2(np.sqrt(1+e13)*np.sin((shit13)/2.0),np.sqrt(1.0-e13)*np.cos((shit13)/2.0))
	r13 = a13*(1.0-e13*np.cos(shit13))
	ox13 = r13*np.cos(nu13)
	oy13 = r13*np.sin(nu13)
	rx13=(ox13*(np.cos(w13)*np.cos(omega13)-np.sin(w13)*np.cos(i13)*np.sin(omega13))- oy13*(np.sin(w13)*np.cos(omega13)
	+ np.cos(w13)*np.cos(i13)*np.sin(omega13)))
	
	ry13=(ox13*(np.cos(w13)*np.sin(omega13)+np.sin(w13)*np.cos(i13)*np.cos(omega13))+ oy13*(np.cos(w13)*np.cos(i13)
	*np.cos(omega13) - np.sin(w13)*np.sin(omega13)))
	rotx13.append(rx13*np.cos(np.pi/2)-ry13*np.sin(np.pi/2.0))
	roty13.append(rx13*np.sin(np.pi/2)+ry13*np.cos(np.pi/2.0))
readyx13 = np.array(rotx13)
readyy13 = np.array(roty13)
	
for t in arange (t014,t014+p14,0.01):
	M14 = (360.0 - 360.0*(t014-t)/p14)%360.0
	M14 = M14*2.0*np.pi/180
	def f14(Et14):
		return Et14 - e14*np.sin(Et14) - M14
	shit14 = optimize.newton(f14, x0)
	nu14 = 2.0*np.arctan2(np.sqrt(1+e14)*np.sin((shit14)/2),np.sqrt(1-e14)*np.cos((shit14)/2))
	r14 = a14*(1.0-e14*np.cos(shit14))
	ox14 = r14*np.cos(nu14)
	oy14 = r14*np.sin(nu14)
	rx14=(ox14*(np.cos(w14)*np.cos(omega14)-np.sin(w14)*np.cos(i14)*np.sin(omega14))- oy14*(np.sin(w14)*np.cos(omega14)
	+ np.cos(w14)*np.cos(i14)*np.sin(omega14)))
	
	ry14=(ox14*(np.cos(w14)*np.sin(omega14)+np.sin(w14)*np.cos(i14)*np.cos(omega14))+ oy14*(np.cos(w14)*np.cos(i14)
	*np.cos(omega14) - np.sin(w14)*np.sin(omega14)))
	rotx14.append(rx14*np.cos(np.pi/2.0)-ry14*np.sin(np.pi/2.0))
	roty14.append(rx14*np.sin(np.pi/2.0)+ry14*np.cos(np.pi/2.0))
readyx14 = np.array(rotx14)
readyy14 = np.array(roty14)

plt.figure(figsize=(8,9))
plt.xlim(0.5,-0.3)
plt.ylim(-0.45,0.55)
circle=plt.Circle((0,0),0.006, color='black')
plt.gcf().gca().add_artist(circle)
plt.title('Orbits of S stars around Sagitarius A')
plt.xlabel('Offest in RA (arcsec)')
plt.ylabel('Offest in DEC (arcsec)')
plt.grid(color='grey', linestyle='-', linewidth=0.5)
plt.plot(-readyx1, readyy1, '-', linewidth=1.5,color='darkblue', label='S1')	
plt.plot(-readyx2,readyy2, '-', linewidth=1.5, color='red', label='S2')
plt.plot(-readyx8,readyy8, '-', linewidth=1.5, color='darkcyan', label='S8')
plt.plot(-readyx12,readyy12, '-', linewidth=1.5, color='black', label='S12')
plt.plot(-readyx13,readyy13, '-', linewidth=1.5, color='magenta', label='S13')
plt.plot(-readyx14,readyy14, '-', linewidth=1.5, color='darkgreen', label='S14')
plt.legend()
plt.savefig('orbits.png')
plt.show()
