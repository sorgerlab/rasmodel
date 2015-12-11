from rasmodel.scenarios.default import model

import numpy as np
from matplotlib import pyplot as plt
from pysb.integrate import Solver
from pysb import *

from tbidbaxlipo.util import fitting

# Zero out all initial conditions
for ic in model.initial_conditions:
    ic[1].value = 0

# In this first experiment, 0.5 uM (500 nM) of HRAS is used:
model.parameters['HRAS_0'].value = 500.

# We use this expression because we still get fluorescence even if GTP
# is hydrolyzed to GDP:
Expression('HRAS_mGXP_', model.observables['HRAS_mGTP_closed_'] +
                   model.observables['HRAS_mGDP_closed_'])
# We use the parameters calculated for experiments with mGTP at 5C
model.parameters['bind_HRASopen_GTP_kf'].value = 1e-2 # nM^-1 s^-1
model.parameters['bind_HRASopen_GTP_kr'].value = 1e-2 / (6.1e4 * 1e-9) # s^-1
model.parameters['equilibrate_HRASopenGTP_to_HRASclosedGTP_kf'].value = 4.5 #s^-1
model.parameters['equilibrate_HRASopenGTP_to_HRASclosedGTP_kr'].value = 0 #s^-1

plt.ion()

t = np.linspace(0, 10, 1000)
sol = Solver(model, t)

plt.figure()
k_list = []
# Perform titration
mgtp_concs = np.arange(1, 15) * 1000 # nM (1 - 15 uM)
for mgtp_conc in mgtp_concs:
    # Titration of labeled GTP:
    model.parameters['mGTP_0'].value = mgtp_conc
    sol.run()

    # Fit to an exponential function to extract the pseudo-first-order rates
    k = fitting.Parameter(1.)
    def expfunc(t):
        return 500 * (1 - np.exp(-k()*t))
    res = fitting.fit(expfunc, [k], sol.yexpr['HRAS_mGXP_'], t)

    # Plot data and fits
    plt.plot(t, sol.yexpr['HRAS_mGXP_'], color='b')
    plt.plot(t, expfunc(t), color='r')

    # Keep the fitted rate
    k_list.append(k())

# Plot the scaling of rates with mGTP concentration
plt.figure()
plt.plot(mgtp_concs, k_list, marker='o')
plt.ylim(bottom=0)

# Figure 3:
# A constant amount of labeled GDP
model.parameters['mGDP_0'].value = 2.5 * 1000 # nM
model.parameters['mGTP_0'].value = 0

model.parameters['bind_HRASopen_GDP_kf'].value = 1e-2 # nM^-1 s^-1
model.parameters['bind_HRASopen_GDP_kr'].value = 1e-2 / (5.7e4 * 1e-9) # s^-1
model.parameters['equilibrate_HRASopenGDP_to_HRASclosedGDP_kf'].value = 3.2 #s^-1
model.parameters['equilibrate_HRASopenGDP_to_HRASclosedGDP_kr'].value = 5e-7 #s^-1

k_list = []
plt.figure()
gdp_concs = np.arange(0, 22) * 1000 # nM
for gdp_conc in gdp_concs:
    # Titration of unlabeled GDP
    model.parameters['GDP_0'].value = gdp_conc
    sol.run()

    k = fitting.Parameter(1.)
    A = fitting.Parameter(100.)
    def expfunc(t):
        return A() * (1 - np.exp(-k()*t))
    res = fitting.fit(expfunc, [A, k], sol.yexpr['HRAS_mGXP_'], t)

    plt.plot(t, sol.yexpr['HRAS_mGXP_'], color='b')
    plt.plot(t, expfunc(t), color='r')

    k_list.append(k())

plt.figure()
plt.plot(gdp_concs, k_list, marker='o')
plt.ylim(bottom=0)
