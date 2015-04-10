Default initial conditions
==========================

For the time being, we simply iterate over all of the Monomer objects in their
default state and assign them an initial condition of 1. This allows certain
analytical procedures that depend on initial conditions to be run.

::

    # Iterate over every monomer
    for m in model.monomers:
        states_dict = {}

        # Iterate over every site in the monomer
        for s in m.sites:
            # If it's in the site states dict, assign it the first of the
            # listed states
            if s in m.site_states:
                #import pdb; pdb.set_trace()
                states_dict[s] = m.site_states[s][0]
            # Otherwise (e.g., the site is used only for binding) assign it
            # a state None, meaning unbound:
            else:
                states_dict[s] = None

        # Create the initial condition parameter based on the protein name
        initial_value = Parameter('{0}_0'.format(m.name), 1.0)

        # Create the initial condition
        Initial(m(**states_dict), initial_value)

