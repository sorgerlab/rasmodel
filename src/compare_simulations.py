from REM.chen_2009 import model
from pysb.integrate import Solver
from pysb.bng import generate_equations
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib.recfunctions import merge_arrays
import sympy
import scipy
import scipy.io
import os.path


plt.figure()

## Simulate PySB model.
# Implement "fixed" species (they are never consumed).
generate_equations(model)
for i in (5, 6, 7):
    model.odes[i] = sympy.numbers.Zero()
tspan = np.linspace(0, 9000, 9001)
solver = Solver(model, tspan, atol=1e-6, rtol=1e-8)
solver.run()
pysb_sim = merge_arrays((solver.yobs, solver.yexpr), flatten=True)

## Load MATLAB simulation results.
matlab_sim = scipy.io.loadmat(os.path.join(os.path.dirname(__file__),
                                           'matlab_sim.mat'))

## Plot results and differences.
for i, (pobs, mobs, color) in enumerate([('pErbB1', 'Y_perbb1', 'b'),
                                         ('pERK', 'Y_perk', 'g'),
                                         ('pAKT', 'Y_pakt', 'r')]):
    # Plot PySB and MATLAB simulations in the left column.
    plt.subplot(3, 2, (i + 1) * 2 - 1)
    plt.plot(tspan, pysb_sim[pobs], c=color, label='PySB')
    plt.plot(matlab_sim['t'][::10], matlab_sim[mobs][::10], color + 'o',
             ms=4, mec='none', label='MATLAB')
    if i == 0:
        plt.title('Simulation trace')
        plt.legend(loc='lower right')
    plt.ylabel(pobs)
    # Plot difference between simulations in the right column.
    plt.subplot(3, 2, (i + 1) * 2)
    pysb_interp = scipy.interp(matlab_sim['t'], tspan, pysb_sim[pobs])
    diff_pct = (pysb_interp - matlab_sim[mobs]) / max(matlab_sim[mobs]) * 100
    plt.plot(matlab_sim['t'], diff_pct, 'm')
    plt.xlim(xmax=500)
    if i == 0:
        plt.title('Difference (%)')
plt.show()
