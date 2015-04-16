Default initial conditions
==========================

Proteins
--------

Our default state for unknown initial protein concentrations and states is to
iterate over all of the Monomer objects in their default state and assign them
an initial condition of 1e-8 (10 nanomolar). This allows certain analytical
procedures that depend on initial conditions to be run.

Afterwards, we fill in more specific concentration information below.

::

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

Nucleotides
-----------

For the concentration of NTPs:

    [PMID7877593]_: The concentrations of bases, nucleosides, and nucleosides
    mono-, di- and tri-phosphate are compared for about 600 published values.
    The data are predominantly from mammalian cells and fluids. For the most
    important ribonucleotides, average concentrations + SD (uM) are: ATP, 3,152
    + 1,698; GTP, 468 + 224;...For deoxynucleosides-triphosphate (dNTP), the
    concentrations in dividing cells are: dATP, 24 + 22; dGTP, 5.2 + 4.5;...
    For the 4 dNTPs, tumor cells have concentrations of 6--11 fold over
    normal cells, and for the 4 NTPs, tumor cells also have concentrations
    1.2-5 fold over the normal cells. 

::

    GTP_0.value = 468e-6

For the time being, we assume a ten-fold lower value for GDP::

    GDP_0.value = GTP_0.value / 10.

For simulations involving a labeled nucleotide (as in the experiments of
Wittinghofer in [PMID2200519]_ and [PMID9585556]_), we explicitly declare these
here as having an initial condition of 0 so that they can be set for specific
simulations later::

    Parameter('mGTP_0', 0.)
    Initial(GTP(p=None, label='y'), mGTP_0)

    Parameter('mGDP_0', 0.)
    Initial(GDP(p=None, label='y'), mGDP_0)

References
----------

.. [PMID7877593] Traut TW. **Physiological concentrations of purines and pyrimidines.** Mol Cell Biochem. 1994 Nov 9;140(1):1-22. :pmid:`7877593` :download:`PDF </pdf/7877593.pdf>`

