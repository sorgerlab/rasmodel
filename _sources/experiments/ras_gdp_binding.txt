.. component::
   :module: rasmodel.experiments.ras_gdp_binding

Ras binding to GTP and GDP
==========================

As a simple test to make sure everything is working, we run a simulation of a
subset of the default model showing the binding of HRAS to GTP and GDP coming to
steady state.

The following experiments are drawn from [PMID2200519]_.

First we import the default model and take care of some preliminaries::

    from rasmodel.scenarios.default import model

    import numpy as np
    from matplotlib import pyplot as plt
    from pysb.integrate import Solver
    from pysb import *

    for ic in model.initial_conditions:
        ic_param.value = 0

    t = np.logspace(-3, 3, 1000)

    sol = Solver(model, t)
    sol.run()


In Figure 2 of [PMID2200519]_, a titration of mGTP is performed, with the
pseudo-first-order rate constant for Ras-mGTP binding calculated at each
concentration by fitting to a single exponential.



