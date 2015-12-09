.. component::
   :module: rasmodel.scenarios.vemurafenib_resistance
	    
Vemurafenib resistance
======================

::
   
    from pysb import *

    Model()

    import rasmodel.components.egfr as egfr
    egfr.egfr_minimal()

    import rasmodel.components.ras as ras
    ras.ras_minimal()

    import rasmodel.components.raf as raf
    raf.raf_monomers()
    raf.ras_activates_raf()
    raf.vemurafenib_monomers()
    raf.vemurafenib_binds_raf()
    raf.mek_phosphorylation()

    import rasmodel.components.mapk as mapk
    mapk.erk_dynamics()

Default initial conditions for model set to 1e-8

::

   from pysb import Parameter, Initial

    # Iterate over every monomer
    for m in model.monomers:
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
        initial_value = Parameter('{0}_0'.format(m.name), 1.0e-8)

        # Create the initial condition
        Initial(m(**states_dict), initial_value)
  
