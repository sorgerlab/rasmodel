from rasmodel.chen_2009 import model
from pysb.integrate import Solver
from pysb.bng import generate_equations
import numpy as np
import matplotlib.pyplot as plt
import sympy

# Replicate matlab simulation

# Implement "fixed" species (they are never consumed).
generate_equations(model)
for i in (5, 6, 7):
    model.odes[i] = sympy.numbers.Zero()

tspan = np.linspace(0, 9000, 9001)
solver = Solver(model, tspan, atol=1e-6, rtol=1e-8)
solver.run()

plt.figure()
for i, (arr, obs, color) in enumerate([(solver.yexpr, 'pErbB1', 'b'),
                                       (solver.yobs, 'pERK', 'g'),
                                       (solver.yobs, 'pAKT', 'r')]):
    plt.subplot(3, 1, i + 1)
    plt.plot(tspan, arr[obs], c=color, label=obs)
plt.show()
