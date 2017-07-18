import numpy as np
import matplotlib.pyplot as plt
with open('output_primal.bin', 'rb') as f:
	u1 = f.read()
Ng = 10
eta1 = np.zeros(Ng)
eta2 = np.zeros(Ng)
u1 = np.fromstring(u1,dtype='<f8')

eta1[:] = u1[0:Ng]
eta2[:] = u1[Ng:2*Ng]
print("eta1",eta1)
print("eta2",eta2)
x = np.linspace(0.,1.,Ng)
cospix = lambda x: [np.cos(np.pi*j*x) for j in range(1,Ng+1)]
sinpix = lambda x: [np.sin(np.pi*j*x) for j in range(1,Ng+1)]
px = [ -1.e0*np.dot(eta2,sinpix(x[j])) for j in range(Ng)]
ux = [ np.dot(eta1,cospix(x[j])) for j in range(Ng)]

plt.figure()
plt.plot(x,ux,linewidth=2.0)
plt.xlabel('x')
plt.ylabel('u(x)')
plt.show()
plt.figure()
plt.plot(x,px,linewidth=2.0)
plt.xlabel('x')
plt.ylabel('p(x)')
plt.show()



