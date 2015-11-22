Ras binding to GTP and GDP
==========================

As a simple test to make sure everything is working, we run a simulation
showing the binding of HRAS to GTP and GDP coming to steady state.

First we take care of some preliminaries::

    import numpy as np
    from matplotlib import pyplot as plt
    from rasmodel import model
    from pysb.integrate import Solver
    from pysb import *

Now we reset all initial conditions, other than the ones we're interested in,
to 0::

    ics_to_keep = ['GTP_0', 'GDP_0', 'HRAS_0']
    for ic in model.initial_conditions:
        ic_param = ic[1]
        if ic_param.name not in ics_to_keep:
            ic_param.value = 0

Set the timepoints of our simulation, in seconds::

    t = np.logspace(-3, 3, 1000)

Initialize the solver object and run the simulation::

    sol = Solver(model, t)
    sol.run()

Plot the fraction of HRAS-GTP vs. HRAS-GDP over time::

    plt.figure()
    HRAS_0 = model.parameters['HRAS_0'].value
    for obs_name in ['HRAS_GTP_closed_', 'HRAS_GTP_open_',
                     'HRAS_GDP_closed_', 'HRAS_GDP_open_', 'HRAS_nf_']:
        plt.plot(t, sol.yobs[obs_name] / HRAS_0, label=obs_name)
    ax = plt.gca()
    ax.set_xscale('log')
    plt.legend(loc='right')
    plt.xlabel('Time (sec)')
    plt.ylabel('Fraction of HRAS')
    plt.savefig('output/simulation_1.png')

.. image:: /../output/simulation_1.png

Dissociation of Ras from GTP and GDP
------------------------------------

To show the intrinsic rate of dissociation of GTP from GDP, we start the simulation
with HRAS tightly bound to GTP (i.e., in the 'open' state) and monitor the fraction of
HRAS bound over time.

First we reset all initial conditions::

    for ic in model.initial_conditions:
        ic_param = ic[1]
        ic_param.value = 0

Now we create a new initial condition for HRAS-GTP::

    HRAS = model.monomers['HRAS']
    GTP = model.monomers['GTP']
    HRAS_mGTP_0 = Parameter('HRAS_mGTP_0', 4e-6)
    model.parameters['GTP_0'].value = 10e-6 # Unlabeled competitor
    model.initial(HRAS(gtp=1, s1s2='open', gef=None, p_loop=None,
                       CAAX=None, oncogenic='n') % GTP(p=1, label='y'),
                  HRAS_mGTP_0)

Now we run the simulation and plot::

    t = np.logspace(1, 6, 1000)
    sol = Solver(model, t)
    sol.run()

    plt.figure()
    plt.plot(t, sol.yobs['HRAS_mGTP_'] / HRAS_mGTP_0.value)
    plt.ylim(0, 1.05)
    ax = plt.gca()
    ax.set_xscale('log')

Repeat the above for GDP::

    for ic in model.initial_conditions:
        ic_param = ic[1]
        ic_param.value = 0

    GDP = model.monomers['GDP']
    HRAS_mGDP_0 = Parameter('HRAS_mGDP_0', 4e-6)
    model.parameters['GDP_0'].value = 10e-6 # Unlabeled competitor
    model.initial(HRAS(gtp=1, s1s2='open', gef=None, p_loop=None,
                       CAAX=None, oncogenic='n') % GDP(p=1, label='y'),
                  HRAS_mGDP_0)

    sol = Solver(model, t)
    sol.run()

Plot on the same plot::

    plt.plot(t, sol.yobs['HRAS_mGDP_'] / HRAS_mGDP_0.value)
    plt.ylim(0, 1.05)
    ax = plt.gca()
    ax.set_xscale('log')

    plt.savefig('output/simulation_2.png')

.. image:: /../output/simulation_2.png

