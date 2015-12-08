.. _raf_minimal:

RAF (minimal)
=============

This passage was processed and assembled automatically by INDRA:

  RAF binds to the RAS-GTP complex.
	The RAF-RAS complex binds another RAF-RAS complex.


INDRA-assembled model components
--------------------------------

::

    def raf_minimal():
        Monomer('RAF', ['ras', 'raf'])
        Monomer('GTP', ['ras'])
        Monomer('RAS', ['gtp', 'raf'])

        Parameter('kf_rr_bind_1', 1e-06)
        Parameter('kr_rr_bind_1', 1e-06)
        Parameter('kf_rr_bind_2', 1e-06)
        Parameter('kr_rr_bind_2', 1e-06)

        Rule('RAF_RAS_GTP_bind', RAF(ras=None) + RAS(gtp=2, raf=None) % GTP(ras=2) <> RAF(ras=1) % RAS(gtp=2, raf=1) % GTP(ras=2), kf_rr_bind_1, kr_rr_bind_1)
        Rule('RAF_RAS_RAF_RAS_bind', RAF(ras=2, raf=None) % RAS(raf=2) + RAF(ras=3, raf=None) % RAS(raf=3) <> RAF(ras=2, raf=1) % RAS(raf=2) % RAF(ras=3, raf=1) % RAS(raf=3), kf_rr_bind_2, kr_rr_bind_2)

        Annotation(RAS, 'http://identifiers.org/pfam/PF00071.18', 'is')
