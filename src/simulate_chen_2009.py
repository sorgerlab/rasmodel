from REM.chen_2009 import model
from pysb.integrate import Solver
import numpy as np
import matplotlib.pyplot as plt

# Replicate matlab simulation

tspan = np.linspace(0, 9000, 901)
s = Solver(model, tspan, atol=1e-6, rtol=1e-8)
s.run()
