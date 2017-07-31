import os
import sys
import shutil
import tempfile
import numpy as np
from scipy.stats import linregress
from subprocess import *
import matplotlib.pyplot as plt

my_path = os.path.dirname(os.path.abspath(__file__))
solver = os.path.join(my_path, 'solver_attractor')


def solve(u, s, nsteps):
    tmp_path = tempfile.mkdtemp()
    with open(os.path.join(tmp_path, 'input_attractor.bin'), 'wb') as f:
        f.write(np.asarray(u, dtype='>d').tobytes())
    with open(os.path.join(tmp_path, 'param.bin'), 'wb') as f:
        f.write(np.asarray(s, dtype='>d').tobytes())
        f.write(np.asarray([0.05,0.865,0.3,0.04], dtype='>d').tobytes())
    call([solver, str(int(nsteps))], cwd=tmp_path)
    with open(os.path.join(tmp_path, 'pf_history.bin'), 'rb') as f:
        out = np.frombuffer(f.read(), dtype='>d')
    shutil.rmtree(tmp_path)
    return out

fig, ax = plt.subplots()
Nsteps = 50000
Ns = 150
k = 0
p = np.zeros(Ns)
ss = np.linspace(0.005,0.035,Ns)
for s in ss:
    eta = np.zeros(31)
    eta[0] = 2.0
    eta[10] = 2.0
    u = np.zeros(Ns)
    pf_hist = solve(eta,s,Nsteps)
    print(max(pf_hist),min(pf_hist),np.mean(pf_hist))
    p[k] = np.mean(pf_hist[5000:])
    k = k + 1

slope, intercept, r_value, p_value, std_err = linregress(ss, p)
print("slope is: ",slope) 
ax.plot(ss,p,"o",linewidth=1.0)
plt.legend(loc='upper left')
fig.savefig('../../demo/plots/pf_vs_c1_qiqi.eps', facecolor=fig.get_facecolor(), edgecolor='none', format='eps', dpi=2000)
