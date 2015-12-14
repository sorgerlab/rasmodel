.. component::
   :module: rasmodel.experiments.kras_mapk

Observables
-----------

Our main readout is the amount of dobule-phosphorylated MAP2K1::

    Observable('RAS_GTP', RAS(gtp=1) % GTP(p=1))
    Observable('RAS_RASGAP', RAS(gap=1) % RASA1(rasgap=1))
    Observable('RAS_RAF', RAS(s1s2=1) % RAF(ras=1))
    Observable('RAFd', RAF(raf=1) % RAF(raf=1))
    Observable('MEKpp', MAP2K1(S218='p', S222='p'))

Simulation
----------

Set up the simulation conditions::

    from pysb.integrate import Solver
    import numpy

    ts = numpy.linspace(0, 1000, 100)
    solver = Solver(model, ts)
    solver.run()

Plotting
--------

::

    import matplotlib.pyplot as plt
    for obs in model.observables:
        plt.plot(ts, solver.yobs[obs.name], label=obs.name)
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration (nM)')
    plt.legend()
    plt.show()
