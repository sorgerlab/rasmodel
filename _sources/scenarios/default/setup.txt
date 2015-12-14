Setup
=====

First we initialize the model::

    from pysb import Model

    Model()

Then we pull in the components for Ras::

    from rasmodel.components import ras

    ras.ras_monomers()
    ras.rasgef_monomers()
    ras.rasgap_monomers()
    ras.nucleotide_monomers()

We include binding and hydrolysis for HRAS, NRAS, and KRAS::

    ras.hras_binds_nucleotides(model)
    ras.hras_hydrolyzes_gtp(model)
    ras.kras_binds_nucleotides(model)
    ras.kras_hydrolyzes_gtp(model)
    ras.nras_binds_nucleotides(model)
    ras.nras_hydrolyzes_gtp(model)

Because some in vitro experiments also contain RasGAPs (e.g., to measure
hydrolysis rates), we include these as well::

    ras.kras_rasgaps(model)
