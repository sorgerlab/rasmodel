.. _mapk_minimal:

MAPK (INDRA)
============

Model description
-----------------

    RAF, bound to RAF, phosphorylates MEK1 at Ser218 and Ser222.
    MEK1, phosphorylated at Ser218 and Ser222, is activated.
    Active MEK1 phosphorylates ERK2 at Tyr187 and Thr185.
    PP2A dephosphorylates MEK1 at Ser218 and Ser222.
    DUSP6 dephosphorylates ERK2 at Tyr187 and Thr185.

INDRA-assembled model
---------------------

::

    def mapk_minimal():
        Monomer(u'DUSP6', [u'MAPK1'])
        Monomer(u'MAP2K1', ['S218', u'RAF', 'S222'],
            {'S218': ['u', 'p'], 'S222': ['u', 'p']})
        Monomer(u'MAPK1', ['Y187', u'MAP2K1', 'T185'],
            {'Y187': ['u', 'p'], 'T185': ['u', 'p']})
        Monomer(u'RAF', [u'RAF', 'phospho', u'MAP2K1'],
            {'phospho': ['u','p']})

        Parameter(u'kf_r_autophos_1', 0.001)
        Parameter(u'kf_rm_bind_1', 1e-06)
        Parameter(u'kr_rm_bind_1', 0.001)
        Parameter(u'kc_rm_phos_1', 0.001)
        Parameter(u'kf_rm_bind_2', 1e-06)
        Parameter(u'kr_rm_bind_2', 0.001)
        Parameter(u'kc_rm_phos_2', 0.001)
        Parameter(u'kf_mm_bind_1', 1e-06)
        Parameter(u'kr_mm_bind_1', 0.001)
        Parameter(u'kc_mm_phos_1', 0.001)
        Parameter(u'kf_mm_bind_2', 1e-06)
        Parameter(u'kr_mm_bind_2', 0.001)
        Parameter(u'kc_mm_phos_2', 0.001)
        Parameter(u'kf_dm_bind_1', 1e-06)
        Parameter(u'kr_dm_bind_1', 0.001)
        Parameter(u'kc_dm_dephos_1', 0.001)
        Parameter(u'kf_dm_bind_2', 1e-06)
        Parameter(u'kr_dm_bind_2', 0.001)
        Parameter(u'kc_dm_dephos_2', 0.001)
        Parameter(u'DUSP6_0', 100.0)
        Parameter(u'MAP2K1_0', 100.0)
        Parameter(u'MAPK1_0', 100.0)
        Parameter(u'RAF_0', 100.0)

        Rule(u'RAF_phospho_bind_MAP2K1_S218_1', RAF(RAF=ANY) +
            MAP2K1(S218='u', RAF=None) <> RAF(RAF=ANY, MAP2K1=1) %\
            MAP2K1(S218='u', RAF=1), kf_rm_bind_1, kr_rm_bind_1)
        Rule(u'RAF_phospho_bind_MAP2K1_S222_1', RAF(RAF=ANY) +
            MAP2K1(S222='u', RAF=None) <> RAF(RAF=ANY, MAP2K1=1) %\
            MAP2K1(S222='u', RAF=1), kf_rm_bind_1, kr_rm_bind_1)
        Rule(u'RAF_phospho_MAP2K1_S218_1', RAF(RAF=ANY, MAP2K1=1) %\
            MAP2K1(S218='u', RAF=1) >> RAF(RAF=ANY, MAP2K1=None) +\
            MAP2K1(S218='p', RAF=None), kc_rm_phos_1)
        Rule(u'RAF_phospho_MAP2K1_S222_1', RAF(RAF=ANY, MAP2K1=1) %\
            MAP2K1(RAF=1, S222='u') >> RAF(RAF=ANY, MAP2K1=None) +\
            MAP2K1(RAF=None, S222='p'), kc_rm_phos_2)
        Rule(u'MAP2K1_phospho_bind_MAPK1_Y187_1', MAP2K1(MAPK1=None) +\
            MAPK1(Y187='u', MAP2K1=None) <> MAP2K1(MAPK1=1) %\
            MAPK1(Y187='u', MAP2K1=1), kf_mm_bind_1, kr_mm_bind_1)
        Rule(u'MAP2K1_phospho_bind_MAPK1_T185_1', MAP2K1(MAPK1=None) +\
            MAPK1(T185='u', MAP2K1=None) <> MAP2K1(MAPK1=1) %\
            MAPK1(T185='u', MAP2K1=1), kf_mm_bind_1, kr_mm_bind_1)
        Rule(u'MAP2K1_phospho_MAPK1_Y187_1', MAP2K1(MAPK1=1) %\
            MAPK1(Y187='u', MAP2K1=1) >> MAP2K1(MAPK1=None) +\
            MAPK1(Y187='p', MAP2K1=None), kc_mm_phos_1)
        Rule(u'MAP2K1_phospho_MAPK1_T185_1', MAP2K1(MAPK1=1) %\
            MAPK1(MAP2K1=1, T185='u') >> MAP2K1(MAPK1=None) +\
            MAPK1(MAP2K1=None, T185='p'), kc_mm_phos_2)


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
