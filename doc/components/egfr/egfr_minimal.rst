.. _egfr_minimal:

EGFR (minimal)
==============

Model description
-----------------
This passage was processed and assembled automatically by INDRA:

    The receptor tyrosine kinase EGFR binds the growth factor ligand EGF. The EGFR-EGF complex binds another EGFR-EGF complex.
    EGFR, bound to EGFR, transphosphorylates itself on a tyrosine residue.
    The adaptor protein GRB2 can bind EGFR that is phosphorylated on tyrosine. EGFR-bound GRB2 binds SOS1.

INDRA-assembled model component
-------------------------------

::

    def egfr_minimal():
        Monomer(u'GRB2', [u'EGFR', u'SOS1'])
        Monomer(u'EGFR', [u'EGF', u'EGFR', u'GRB2', 'Y'], {'Y': ['u', 'p']})
        Monomer(u'EGF', [u'EGFR'])
        Monomer(u'SOS1', [u'GRB2'])

        Parameter(u'kf_ee_bind_1', 1e-06)
        Parameter(u'kr_ee_bind_1', 1e-06)
        Parameter(u'kf_ee_bind_2', 1e-06)
        Parameter(u'kr_ee_bind_2', 1e-06)
        Parameter(u'kf_eg_bind_1', 1e-06)
        Parameter(u'kr_eg_bind_1', 1e-06)
        Parameter(u'kf_gs_bind_1', 1e-06)
        Parameter(u'kr_gs_bind_1', 1e-06)
        Parameter(u'kf_ee_transphos_1', 0.001)
        Parameter(u'GRB2_0', 100.0)
        Parameter(u'EGFR_0', 100.0)
        Parameter(u'EGF_0', 100.0)
        Parameter(u'SOS1_0', 100.0)

        Rule(u'EGFR_EGF_bind', EGFR(EGF=None) + EGF(EGFR=None) <>
                EGFR(EGF=1) % EGF(EGFR=1), kf_ee_bind_1, kr_ee_bind_1)
        Rule(u'EGFR_EGF_EGFR_EGF_bind', EGFR(EGF=2, EGFR=None) % EGF(EGFR=2) +
                EGFR(EGF=3, EGFR=None) % EGF(EGFR=3) <> EGFR(EGF=2, EGFR=1) %\
                EGF(EGFR=2) % EGFR(EGF=3, EGFR=1) % EGF(EGFR=3), kf_ee_bind_2, kr_ee_bind_2)
        Rule(u'EGFR_transphospho_EGFR_Y', EGFR(EGFR=ANY) % EGFR(Y='u') >>
                EGFR(EGFR=ANY) % EGFR(Y='p'), kf_ee_transphos_1)
        Rule(u'GRB2_EGFR_bind', GRB2(EGFR=None) + EGFR(GRB2=None, Y='p') <>
              GRB2(EGFR=1) % EGFR(GRB2=1), kf_eg_bind_1, kr_eg_bind_1)
        Rule(u'GRB2_EGFR_SOS1_bind', GRB2(EGFR=2, SOS1=None) % EGFR(GRB2=2) +\
                SOS1(GRB2=None) <> GRB2(EGFR=2, SOS1=1) % EGFR(GRB2=2) % SOS1(GRB2=1), kf_gs_bind_1, kr_gs_bind_1)

        Annotation(GRB2, 'http://identifiers.org/uniprot/P62993', 'is')
        Annotation(GRB2, 'http://identifiers.org/hgnc/HGNC:4566', 'is')
        Annotation(EGFR, 'http://identifiers.org/uniprot/P00533', 'is')
        Annotation(EGFR, 'http://identifiers.org/hgnc/HGNC:3236', 'is')
        Annotation(EGF, 'http://identifiers.org/uniprot/P01133', 'is')
        Annotation(EGF, 'http://identifiers.org/hgnc/HGNC:3229', 'is')
        Annotation(SOS1, 'http://identifiers.org/uniprot/Q07889', 'is')
        Annotation(SOS1, 'http://identifiers.org/hgnc/HGNC:11187', 'is')

    #

.. raw:: html

    <script>
        window.setTimeout(function() {
        $('div.highlight-python pre > span.c:last-child').each(
            function () {
                if ($(this).text() == '#') {
                    $(this.nextSibling).detach();
                    $(this).detach();
                }
            }
        );
        }, 1000);
    </script>
