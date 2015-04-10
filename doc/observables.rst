Observables
===========

The list of observables allows us to identify the sets of molecular species
that we are interested in monitoring during simulation. As a simple convention,
we terminate each name with an underscore to identify it as an observable. For
example, while 'HRAS' denotes the Monomer component for HRAS, 'HRAS_' might
denote the Observable component for HRAS.

::

    Observable('HRAS_GTP_', HRAS(gtp=1) % GTP(p=1))
    Observable('HRAS_GDP_', HRAS(gtp=1) % GDP(p=1))

