import os
import sys
import shutil
import tempfile
import numpy as np
from subprocess import *
import matplotlib.pyplot as plt

my_path = os.path.dirname(os.path.abspath(__file__))
solver = os.path.join(my_path, 'solver_phaseportrait')
u0 = np.loadtxt(os.path.join(my_path, 'u0'))

def solve(u, s, nsteps):
    tmp_path = tempfile.mkdtemp()
    with open(os.path.join(tmp_path, 'input_phaseportrait.bin'), 'wb') as f:
        f.write(np.asarray(u, dtype='>d').tobytes())
    with open(os.path.join(tmp_path, 'param.bin'), 'wb') as f:
        f.write(np.asarray(s, dtype='>d').tobytes())
        f.write(np.asarray([0.05,0.865,0.3,0.04], dtype='>d').tobytes())
    call([solver, str(int(nsteps))], cwd=tmp_path)
    with open(os.path.join(tmp_path, 'output_phaseportrait.bin'), 'rb') as f:
        out = np.frombuffer(f.read(), dtype='>d')
    shutil.rmtree(tmp_path)
    return out

fig, ax = plt.subplots(1,2)
fig.set_facecolor('black')
plt.setp([a.get_xaxis().set_visible(False) for a in ax[:]])
plt.setp([a.get_yaxis().set_visible(False) for a in ax[:]])
plt.setp([a.set_facecolor('black') for a in ax[:]])
s = 0.05
Nsteps = 50000
Nplotsteps = Nsteps
k1 = 0
colors=["r","b"]

for i in range(2):
    eta = np.zeros(20)
    if(i!=0):
        eta[0] = 0.1*(i+20)
        eta[10] = 0.1*(i+20)
    else:
        eta = np.random.rand(20)
    p = np.zeros(Nplotsteps)
    u = np.zeros(Nplotsteps)
    heat = np.zeros(Nplotsteps)
    u1 = solve(eta,s,Nsteps)
    print(max(u1),min(u1),np.mean(u1))
    k = 0
    for n in range(Nsteps-Nplotsteps,Nsteps):
        p[k] = u1[3*n]
        u[k] = u1[(3*n+1)]
        heat[k] = u1[3*n+2]
        k = k + 1
    ax[k1].plot(p[5000:],u[5000:],linewidth=1.0,color=colors[i])
    k1 = k1+1
#plt.legend(loc='upper left')
fig.savefig('../../demo/plots/juniper_attractor.eps', facecolor=fig.get_facecolor(), edgecolor='none', format='eps', dpi=2000)
