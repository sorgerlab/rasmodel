Observables
===========

The list of observables allows us to identify the sets of molecular species
that we are interested in monitoring during simulation. As a simple convention,
we terminate each name with an underscore to identify it as an observable.

::

    from pysb import Observable

    Observable('HRAS_GTP_', HRAS(gtp=1) % GTP(p=1))
    Observable('HRAS_mGTP_', HRAS(gtp=1) % GTP(p=1, label='y'))
    Observable('HRAS_GTP_closed_', HRAS(gtp=1, s1s2='closed') % GTP(p=1))
    Observable('HRAS_GTP_open_', HRAS(gtp=1, s1s2='open') % GTP(p=1))
    Observable('HRAS_GDP_', HRAS(gtp=1) % GDP(p=1))
    Observable('HRAS_mGDP_', HRAS(gtp=1) % GDP(p=1, label='y'))
    Observable('HRAS_GDP_closed_', HRAS(gtp=1, s1s2='closed') % GDP(p=1))
    Observable('HRAS_GDP_open_', HRAS(gtp=1, s1s2='open') % GDP(p=1))
    Observable('HRAS_nf_', HRAS(gtp=None))

Observables for overall nucleotide abundances::

    Observable('GTP_', GTP())
    Observable('GDP_', GDP())
    Observable('Pi_', Pi())

