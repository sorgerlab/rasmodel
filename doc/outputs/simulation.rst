Sample simulation
=================

As a simple test to make sure everything is working, we run a simulation
showing the binding of HRAS to GTP and GDP coming to steady state.

First we take care of some preliminaries::

    import numpy as np
    from ras_model import model
    from pysb.integrate import Solver

Now we reset all initial conditions, other than the ones we're interested in,
to 0::

    ics_to_keep = ['GTP_0', 'GDP_0', 'KRAS_0']
    for ic in model.initial_conditions:
        ic_param = ic[1]
        if ic_param.name not in ics_to_keep:
            ic_param.value = 0

Set the timepoints of our simulation, in seconds::

    t = np.linspace(0, 1e2)

Initialize the solver object and run the simulation::

    sol = Solver(model, t)
    sol.run()

Plot the amount of HRAS-GTP over time::

    plt.figure()
    plt.plot(t, sol.yobs['HRAS_GTP'])
    plt.savefig('simulation_1.png')

.. image:: /../output/simulation_1.png

