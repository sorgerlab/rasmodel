.. _raf_minimal:

RAF (minimal)
=============

This passage was processed and assembled automatically by INDRA:

    RAF binds to the RAS-GTP complex.
    The RAF-RAS complex binds another RAF-RAS complex.


INDRA-assembled model components
--------------------------------

::

    from pysb import Monomer, Parameter, Rule, Annotation, ANY
    from pysb.util import alias_model_components

    def raf_monomers():
        Monomer('RAF', ['raf', 'ras', 'map2k1'])

    def raf_minimal(model):
        RAF = model.monomers['RAF']

        Parameter('kf_rr_bind_1', 1e-06)
        Parameter('kr_rr_bind_1', 1e-06)
        Parameter('kf_rr_bind_2', 1e-06)
        Parameter('kr_rr_bind_2', 1e-06)

        alias_model_components()

        Rule('RAF_RAS_GTP_bind', RAF(ras=None) + RAS(gtp=ANY, s1s2=None) >>
            RAF(ras=1) % RAS(gtp=ANY, s1s2=1), kf_rr_bind_1)

        Rule('RAF_RAS_GTP_dissociate', RAF(ras=1) % RAS(s1s2=1) >>
            RAF(ras=None) + RAS(s1s2=None), kr_rr_bind_1)

        Rule('RAF_RAS_RAF_RAS_bind', RAF(ras=ANY, raf=None) +\
            RAF(ras=ANY, raf=None) >> RAF(ras=ANY, raf=1) %\
            RAF(ras=ANY, raf=1), kf_rr_bind_2)

        Rule('RAF_RAS_RAF_RAS_dissociate', RAF(raf=1) % RAF(raf=1) >>
            RAF(raf=None) + RAF(raf=None), kr_rr_bind_2)

        Annotation(RAS, 'http://identifiers.org/pfam/PF00071.18', 'is')
