Observables
===========

The list of observables allows us to identify the sets of molecular species
that we are interested in monitoring during simulation. As a simple convention,
we terminate each name with an underscore to identify it as an observable.

::

    from pysb import Observable

    # Create the following observables for HRAS, KRAS, NRAS:
    for ras in [HRAS, KRAS, NRAS]:
        Observable('%s_GTP_' % ras.name,
                   ras(gtp=1) % GTP(p=1))
        Observable('%s_mGTP_' % ras.name,
                   ras(gtp=1) % GTP(p=1, label='y'))
        Observable('%s_GTP_closed_' % ras.name,
                   ras(gtp=1, s1s2='closed') % GTP(p=1))
        Observable('%s_GTP_open_' % ras.name,
                   ras(gtp=1, s1s2='open') % GTP(p=1))
        Observable('%s_mGTP_closed_' % ras.name,
                   ras(gtp=1, s1s2='closed') % GTP(p=1, label='y'))
        Observable('%s_mGTP_open_' % ras.name,
                   ras(gtp=1, s1s2='open') % GTP(p=1, label='y'))
        Observable('%s_GDP_' % ras.name,
                   ras(gtp=1) % GDP(p=1))
        Observable('%s_mGDP_' % ras.name,
                   ras(gtp=1) % GDP(p=1, label='y'))
        Observable('%s_GDP_closed_' % ras.name,
                   ras(gtp=1, s1s2='closed') % GDP(p=1))
        Observable('%s_GDP_open_' % ras.name,
                   ras(gtp=1, s1s2='open') % GDP(p=1))
        Observable('%s_mGDP_closed_' % ras.name,
                   ras(gtp=1, s1s2='closed') % GDP(p=1, label='y'))
        Observable('%s_mGDP_open_' % ras.name,
                   ras(gtp=1, s1s2='open') % GDP(p=1, label='y'))
        Observable('%s_nf_' % ras.name, ras(gtp=None))

Observables for overall nucleotide abundances::

    Observable('GTP_', GTP())
    Observable('GDP_', GDP())
    Observable('Pi_', Pi())

