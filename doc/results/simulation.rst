Sample simulation
=================

As a simple test to make sure everything is working, we run a simulation
showing the binding of HRAS to GTP and GDP coming to steady state.

First we take care of some preliminaries::

    import numpy as np
    from matplotlib import pyplot as plt
    from ras_model import model
    from pysb.integrate import Solver

Now we reset all initial conditions, other than the ones we're interested in,
to 0::

    ics_to_keep = ['GTP_0', 'GDP_0', 'HRAS_0']
    for ic in model.initial_conditions:
        ic_param = ic[1]
        if ic_param.name not in ics_to_keep:
            ic_param.value = 0

Set the timepoints of our simulation, in seconds::

    t = np.logspace(-3, 6, 1000)

Initialize the solver object and run the simulation::

    sol = Solver(model, t)
    sol.run()

Plot the fraction of HRAS-GTP vs. HRAS-GDP over time::

    plt.figure()
    HRAS_0 = model.parameters['HRAS_0'].value
    for obs_name in ['HRAS_GTP_', 'HRAS_GDP_', 'HRAS_nf_']:
        plt.plot(t, sol.yobs[obs_name] / HRAS_0, label=obs_name)
    ax = plt.gca()
    ax.set_xscale('log')
    plt.legend(loc='right')
    plt.xlabel('Time (sec)')
    plt.ylabel('Fraction of HRAS')
    plt.savefig('output/simulation_1.png')

.. image:: /../output/simulation_1.png

