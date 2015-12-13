.. component::
   :module: rasmodel.scenarios.kras_mapk

MAPK signaling initiated by mutant KRAS
=======================================

Introduction
------------

Mutations at the G12 and Q61 positions are known to affect KRAS GTP
hydrolysis rates as well as affinity to downstream effectors like
RAF proteins. This model scenario aims to explain differences in the
dynamics of the MAP kinase cascade resulting from KRAS mutations.

Setup
-----

We construct the model out of existing components::

    from pysb import Model, Initial, Parameter, Observable
    from rasmodel.components import ras, raf, mapk

    Model()

    ras.ras_monomers()
    ras.rasgef_monomers()
    ras.rasgap_monomers()
    ras.nucleotide_monomers()

We only include the elements involving KRAS, not HRAS or NRAS::

    ras.kras_binds_nucleotides(model)
    ras.kras_hydrolyzes_gtp(model)
    ras.kras_rasgaps(model)

Ensure GTP/GDP levels are maintained::

    ras.recycle_gtp_from_gdp(model)

Next we include the module for RAF activation. We assume that all isoforms
of RAF are present as wild-type and we model them jointy as RAF::

    raf.raf_monomers()
    RAS = model.monomers['KRAS']
    RAS.rename('RAS')
    RAS.sites.append('raf')
    raf.raf_minimal(model)

Finally, we instantiate the MAP kinase cascade::

    mapk.mapk_minimal()
    #mapk.erk_dynamics()

Setting up parameters
---------------------

These parameters were automatically generated. Here we set them to plausible values::

    model.parameters['kf_rr_bind_1'].value = 1e-2 # nMol^-1 s^-1
    model.parameters['kr_rr_bind_1'].value = 0.56 # s^-1
    model.parameters['kf_rr_bind_2'].value = 1e-2 # nMol^-1 s^-1
    model.parameters['kr_rr_bind_2'].value = 1e-3 # s^-1
    model.parameters['kf_rm_bind_1'].value = 1e-2
    model.parameters['kr_rm_bind_1'].value = 1e-2
    model.parameters['kc_rm_phos_1'].value = 1e-2
    model.parameters['kf_rm_bind_2'].value = 1e-2
    model.parameters['kr_rm_bind_2'].value = 1e-2
    model.parameters['kc_rm_phos_2'].value = 1e-2

Initial conditions
------------------
We first set placeholders for initial conditions and then set the specific values. ::

    #init_nonzero = [RAS, RAF, MAP2K1, MAPK1, PPP2CA, DUSP6, GTP]
    init_nonzero = [RAS, RASA1, GTP, RAF, MAP2K1]
    # Iterate over every monomer
    for m in init_nonzero:
        states_dict = {}
        # Iterate over every site in the monomer
        for s in m.sites:
            # If it's in the site states dict, assign it the first of the
            # listed states
            if s in m.site_states:
                states_dict[s] = m.site_states[s][0]
            # Otherwise (e.g., the site is used only for binding) assign it
            # a state None, meaning unbound:
            else:
                states_dict[s] = None

        # Create the initial condition parameter based on the protein name
        initial_value = Parameter('{0}_0'.format(m.name), 0)

        # Create the initial condition
        Initial(m(**states_dict), initial_value)

The volume of a HeLA cell is around 2425 um^3 [PMID17461436]_. N_A_nm denotes the number of molecules in one nanomol::

    V = 2425e-15
    N_A_nm = 6.022e14
    init_scale = 1.0 / (N_A_nm * V)

    model.parameters['GTP_0'].value = 1e9 * init_scale
    model.parameters['RAS_0'].value = 50000 * init_scale
    model.parameters['RASA1_0'].value = 30000 * init_scale
    model.parameters['RAF_0'].value = 35000 * init_scale
    model.parameters['MAP2K1_0'].value = 80000 * init_scale
    ##model.parameters['MAPK1_0'].value = 100000 * init_scale
    ##model.parameters['PPP2CA_0'].value = 10000 * init_scale
    ##model.parameters['DUSP6_0'].value = 5000 * init_scale

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
