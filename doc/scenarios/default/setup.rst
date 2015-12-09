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

We only include the elements involving HRAS, not NRAS or KRAS::

    ras.hras_binds_nucleotides(model)
    ras.hras_hydrolyzes_gtp(model)

Ensure GTP/GDP levels are maintained::

    ras.recycle_gtp_from_gdp(model)
