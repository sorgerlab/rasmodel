.. _egfr_minimal:

EGFR (minimal)
==============

This passage was processed and assembled automatically by INDRA:

  The receptor tyrosine kinase EGFR binds the growth factor ligand EGF.
  The EGFR-EGF complex binds another EGFR-EGF complex.
  EGFR, bound to EGFR, transphosphorylates itself on a tyrosine residue.
  The adaptor protein GRB2 can bind EGFR that is phosphorylated on tyrosine.
  EGFR-bound GRB2 binds SOS1.

INDRA-assembled model components
--------------------------------

::

    from pysb import Monomer, Parameter, Rule, Annotation, ANY
    from pysb.util import alias_model_components

    def egfr_minimal():
        Monomer(u'GRB2', [u'egfr', u'sos1'])
        Monomer(u'EGFR', [u'egf', u'egfr', 'Y', u'grb2'], {'Y': ['u', 'p']})
        Monomer(u'EGF', [u'egfr'])
        Monomer(u'SOS1', [u'grb2', 'ras'])

        Parameter(u'kf_ee_bind_1', 1e-06)
        Parameter(u'kr_ee_bind_1', 1e-06)
        Parameter(u'kf_ee_bind_2', 1e-06)
        Parameter(u'kr_ee_bind_2', 1e-06)
        Parameter(u'kf_ge_bind_1', 1e-06)
        Parameter(u'kr_ge_bind_1', 1e-06)
        Parameter(u'kf_gs_bind_1', 1e-06)
        Parameter(u'kr_gs_bind_1', 1e-06)
        Parameter(u'kf_ee_transphos_1', 0.001)

        alias_model_components()

        Rule(u'EGFR_EGF_bind', EGFR(egf=None) + EGF(egfr=None) >>
            EGFR(egf=1) % EGF(egfr=1), kf_ee_bind_1)

        Rule(u'EGFR_EGF_dissociate', EGFR(egf=1) % EGF(egfr=1) >>
            EGFR(egf=None) + EGF(egfr=None), kr_ee_bind_1)

        Rule(u'EGFR_EGF_EGFR_EGF_bind', EGFR(egf=ANY, egfr=None) +
            EGFR(egf=ANY, egfr=None) >> EGFR(egf=ANY, egfr=1) %\
            EGFR(egf=ANY, egfr=1), kf_ee_bind_2)

        Rule(u'EGFR_EGF_EGFR_EGF_dissociate', EGFR(egfr=1) % EGFR(egfr=1) >>
            EGFR(egfr=None) + EGFR(egfr=None), kr_ee_bind_2)

        Rule(u'GRB2_EGFR_bind', GRB2(egfr=None) + EGFR(grb2=None) >>
            GRB2(egfr=1) % EGFR(grb2=1), kf_ge_bind_1)

        Rule(u'GRB2_EGFR_dissociate', GRB2(egfr=1) % EGFR(grb2=1) >>
            GRB2(egfr=None) + EGFR(grb2=None), kr_ge_bind_1)

        Rule(u'GRB2_EGFR_SOS1_bind', GRB2(egfr=ANY, sos1=None) +\
            SOS1(grb2=None) >> GRB2(egfr=ANY, sos1=1) % SOS1(grb2=1),
            kf_gs_bind_1)

        Rule(u'GRB2_EGFR_SOS1_dissociate', GRB2(sos1=1) % SOS1(grb2=1) >>
            GRB2(sos1=None) + SOS1(grb2=None), kr_gs_bind_1)

        Rule(u'EGFR_transphospho_EGFR_Y',
            EGFR(egfr=ANY) % EGFR(Y='u') >> EGFR(egfr=ANY) % EGFR(Y='p'),
            kf_ee_transphos_1)


        Annotation(GRB2, 'http://identifiers.org/uniprot/P62993', 'is')
        Annotation(GRB2, 'http://identifiers.org/hgnc/HGNC:4566', 'is')
        Annotation(EGFR, 'http://identifiers.org/uniprot/P00533', 'is')
        Annotation(EGFR, 'http://identifiers.org/hgnc/HGNC:3236', 'is')
        Annotation(EGF, 'http://identifiers.org/uniprot/P01133', 'is')
        Annotation(EGF, 'http://identifiers.org/hgnc/HGNC:3229', 'is')
        Annotation(SOS1, 'http://identifiers.org/pfam/PF00617', 'is')
        Annotation(SOS1, 'http://identifiers.org/uniprot/Q62245', 'is')
        Annotation(SOS1, 'http://identifiers.org/hgnc/HGNC:11187', 'is')
