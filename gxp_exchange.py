from rasmodel.scenarios.default import model

import numpy as np
from matplotlib import pyplot as plt
from pysb.integrate import Solver
from pysb import *

from tbidbaxlipo.util import fitting

# Zero out all initial conditions
for ic in model.initial_conditions:
    ic[1].value = 0

KRAS = model.monomers['KRAS']
GDP = model.monomers['GDP']
GTP = model.monomers['GTP']

Expression('KRAS_mGXP_', model.observables['KRAS_mGTP_open_'] +
                         model.observables['KRAS_mGDP_open_'])

# Add an initial condition for HRAS with GDP or GTP pre-bound
# (Concentration units in nM)
Initial(KRAS(gtp=1, gef=None, p_loop=None, s1s2='open', CAAX=None,
             mutant='WT') % GDP(p=1, label='n'),
             Parameter('KRAS_WT_GDP_0', 0.))
Initial(KRAS(gtp=1, gef=None, p_loop=None, s1s2='open', CAAX=None,
             mutant='G13D') % GDP(p=1, label='n'),
             Parameter('KRAS_G13D_GDP_0', 0.))
Initial(KRAS(gtp=1, gef=None, p_loop=None, s1s2='open', CAAX=None,
             mutant='WT') % GTP(p=1, label='n'),
             Parameter('KRAS_WT_GTP_0', 0.))
Initial(KRAS(gtp=1, gef=None, p_loop=None, s1s2='open', CAAX=None,
             mutant='G13D') % GTP(p=1, label='n'),
             Parameter('KRAS_G13D_GTP_0', 0.))
#model.parameters['k_recycle_gtp_from_gdp'].value = 0
plt.ion()

# First simulate the data from Figure 1A (GDP exchange)
# WT, GDP:
model.parameters['mGDP_0'].value = 1500.
model.parameters['KRAS_WT_GDP_0'].value = 750.
t = np.linspace(0, 1000, 1000) # 1000 seconds
sol = Solver(model, t)
sol.run()
plt.figure()
plt.plot(t, sol.yexpr['KRAS_mGXP_'], label='WT')

k = fitting.Parameter(1.)
A = fitting.Parameter(100.)
def expfunc(t):
    return A() * (1 - np.exp(-k()*t))
res = fitting.fit(expfunc, [A, k], sol.yexpr['KRAS_mGXP_'], t)
print A(), k()

# G13D, GDP:
model.parameters['KRAS_WT_GDP_0'].value = 0
model.parameters['KRAS_G13D_GDP_0'].value = 750.
sol.run()
plt.plot(t, sol.yexpr['KRAS_mGXP_'], label='G13D')
plt.legend(loc='lower right')

k = fitting.Parameter(1.)
A = fitting.Parameter(100.)
def expfunc(t):
    return A() * (1 - np.exp(-k()*t))
res = fitting.fit(expfunc, [A, k], sol.yexpr['KRAS_mGXP_'], t)
print A(), k()

# Now simulate the data from Figure 1B (GTP exchange)
# WT, GTP
model.parameters['mGDP_0'].value = 0.
model.parameters['mGTP_0'].value = 1500.
model.parameters['KRAS_WT_GDP_0'].value = 0.
model.parameters['KRAS_G13D_GDP_0'].value = 0.
model.parameters['KRAS_WT_GTP_0'].value = 750.
model.parameters['KRAS_G13D_GTP_0'].value = 0.
sol.run()
plt.plot(t, sol.yexpr['KRAS_mGXP_'], label='WT')

# G13D, GTP
model.parameters['KRAS_WT_GTP_0'].value = 0.
model.parameters['KRAS_G13D_GTP_0'].value = 750.
sol.run()
plt.plot(t, sol.yexpr['KRAS_mGXP_'], label='G13D')
plt.legend(loc='lower right')

k = fitting.Parameter(1.)
A = fitting.Parameter(100.)
def expfunc(t):
    return A() * (1 - np.exp(-k()*t))
res = fitting.fit(expfunc, [A, k], sol.yexpr['KRAS_mGXP_'], t)
print A(), k()
