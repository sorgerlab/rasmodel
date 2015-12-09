.. _mapk_minimal:

MAPK (minimal)
==============

This passage was processed and assembled automatically by INDRA:

	RAF, bound to RAF, phosphorylates MEK1 at Ser218 and Ser222.
	MEK1, phosphorylated at Ser218 and Ser222, is activated.
	Active MEK1 phosphorylate
	ERK2 at Tyr187 and Thr185.
	PP2A-alpha dephosphorylates MEK1 at Ser218 and Ser222.
	DUSP6 dephosphorylates ERK2 at Tyr187 and Thr185.

INDRA-assembled model components
--------------------------------

::


   from pysb import Monomer, Parameter, Rule, Annotation
   from pysb.util import alias_model_components


   def mapk_minimal():
        
       Monomer(u'MAP2K1', ['S218', u'ppp2ca', 'raf', 'S222', u'mapk1'],
			  {'S218': ['u', 'p'], 'S222': ['u', 'p']})

       Parameter(u'kf_rm_bind_1', 1e-06)
       Parameter(u'kr_rm_bind_1', 0.001)
       Parameter(u'kc_rm_phos_1', 0.001)
       Parameter(u'kf_rm_bind_2', 1e-06)
       Parameter(u'kr_rm_bind_2', 0.001)
       Parameter(u'kc_rm_phos_2', 0.001)

       alias_model_components()

       Rule(u'RAF_phospho_bind_MAP2K1_S218_1',
	    RAF(raf1=ANY, map2k1=None) +
	    MAP2K1(S218='u', raf=None) <>
	    RAF(raf1=ANY, map2k1=1) % MAP2K1(S218='u', raf=1),
	    kf_rm_bind_1, kr_rm_bind_1)

       Rule(u'RAF_phospho_MAP2K1_S218_1',
	    RAF(raf1=ANY, map2k1=1) %
	    MAP2K1(S218='u', raf=1) >>
	    RAF(raf1=ANY, map2k1=None) + MAP2K1(S218='p', raf=None),
	    kc_rm_phos_1)

       Rule(u'RAF_phospho_bind_MAP2K1_S222_1',
	    RAF(raf1=ANY, map2k1=None) + MAP2K1(raf=None, S222='u') <>
	    RAF(raf1=ANY, map2k1=1) % MAP2K1(raf=1, S222='u'),
	    kf_rm_bind_2, kr_rm_bind_2)

       Rule(u'RAF_phospho_MAP2K1_S222_1',
	    RAF(raf1=ANY, map2k1=1) % MAP2K1(raf=1, S222='u') >>
	    RAF(raf1=ANY, map2k1=None) + MAP2K1(raf=None, S222='p'),
	    kc_rm_phos_2)

 
   def erk_dynamics():
   
       Monomer(u'PPP2CA', [u'map2k1'])
       Monomer(u'MAPK1', ['Y187', u'map2k1', 'T185', u'dusp6'],
			  {'Y187': ['u', 'p'], 'T185': ['u', 'p']})
       Monomer(u'DUSP6', [u'mapk1'])

       Parameter(u'kf_mm_bind_1', 1e-06)
       Parameter(u'kr_mm_bind_1', 0.001)
       Parameter(u'kc_mm_phos_1', 0.001)
       Parameter(u'kf_mm_bind_2', 1e-06)
       Parameter(u'kr_mm_bind_2', 0.001)
       Parameter(u'kc_mm_phos_2', 0.001)
       Parameter(u'kf_pm_bind_1', 1e-06)
       Parameter(u'kr_pm_bind_1', 0.001)
       Parameter(u'kc_pm_dephos_1', 0.001)
       Parameter(u'kf_pm_bind_2', 1e-06)
       Parameter(u'kr_pm_bind_2', 0.001)
       Parameter(u'kc_pm_dephos_2', 0.001)
       Parameter(u'kf_dm_bind_1', 1e-06)
       Parameter(u'kr_dm_bind_1', 0.001)
       Parameter(u'kc_dm_dephos_1', 0.001)
       Parameter(u'kf_dm_bind_2', 1e-06)
       Parameter(u'kr_dm_bind_2', 0.001)
       Parameter(u'kc_dm_dephos_2', 0.001)

       alias_model_components()


       Rule(u'MAP2K1_phospho_bind_MAPK1_Y187_1',
	    MAP2K1(S218='p', S222='p', mapk1=None) +
	    MAPK1(Y187='u', map2k1=None) <>
	    MAP2K1(S218='p', S222='p', mapk1=1) % MAPK1(Y187='u', map2k1=1),
	    kf_mm_bind_1, kr_mm_bind_1)

       Rule(u'MAP2K1_phospho_MAPK1_Y187_1',
	    MAP2K1(S218='p', S222='p', mapk1=1) % MAPK1(Y187='u', map2k1=1) >>
	    MAP2K1(S218='p', S222='p', mapk1=None) +
	    MAPK1(Y187='p', map2k1=None),
	    kc_mm_phos_1)

       Rule(u'MAP2K1_phospho_bind_MAPK1_T185_1',
	    MAP2K1(S218='p', S222='p', mapk1=None) +
	    MAPK1(map2k1=None, T185='u') <>
	    MAP2K1(S218='p', S222='p', mapk1=1) % MAPK1(map2k1=1, T185='u'),
	    kf_mm_bind_2, kr_mm_bind_2)

       Rule(u'MAP2K1_phospho_MAPK1_T185_1',
	    MAP2K1(S218='p', S222='p', mapk1=1) % MAPK1(map2k1=1, T185='u') >>
	    MAP2K1(S218='p', S222='p', mapk1=None) +
	    MAPK1(map2k1=None, T185='p'), kc_mm_phos_2)

       Rule(u'PPP2CA_dephos_bind_map2k1_S218_1',
	    PPP2CA(map2k1=None) + MAP2K1(S218='p', ppp2ca=None) <>
	    PPP2CA(map2k1=1) % MAP2K1(S218='p', ppp2ca=1),
	    kf_pm_bind_1, kr_pm_bind_1)

       Rule(u'PPP2CA_dephos_map2k1_S218_1',
	    PPP2CA(map2k1=1) % MAP2K1(S218='p', ppp2ca=1) >>
	    PPP2CA(map2k1=None) + MAP2K1(S218='u', ppp2ca=None),
	    kc_pm_dephos_1)

       Rule(u'PPP2CA_dephos_bind_map2k1_S222_1',
	    PPP2CA(map2k1=None) + MAP2K1(ppp2ca=None, S222='p') <>
	    PPP2CA(map2k1=1) % MAP2K1(ppp2ca=1, S222='p'),
	    kf_pm_bind_2, kr_pm_bind_2)

       Rule(u'PPP2CA_dephos_map2k1_S222_1',
	    PPP2CA(map2k1=1) % MAP2K1(ppp2ca=1, S222='p') >>
	    PPP2CA(map2k1=None) + MAP2K1(ppp2ca=None, S222='u'),
	    kc_pm_dephos_2)

       Rule(u'DUSP6_dephos_bind_MAPK1_Y187_1',
	    DUSP6(mapk1=None) + MAPK1(Y187='p', dusp6=None) <>
	    DUSP6(mapk1=1) % MAPK1(Y187='p', dusp6=1),
	    kf_dm_bind_1, kr_dm_bind_1)

       Rule(u'DUSP6_dephos_MAPK1_Y187_1',
	    DUSP6(mapk1=1) % MAPK1(Y187='p', dusp6=1) >>
	    DUSP6(mapk1=None) + MAPK1(Y187='u', dusp6=None), kc_dm_dephos_1)

       Rule(u'DUSP6_dephos_bind_MAPK1_T185_1',
	    DUSP6(mapk1=None) + MAPK1(T185='p', dusp6=None) <>
	    DUSP6(mapk1=1) % MAPK1(T185='p', dusp6=1),
	    kf_dm_bind_2, kr_dm_bind_2)

       Rule(u'DUSP6_dephos_MAPK1_T185_1',
	    DUSP6(mapk1=1) % MAPK1(T185='p', dusp6=1) >>
	    DUSP6(mapk1=None) + MAPK1(T185='u', dusp6=None), kc_dm_dephos_2)


       Annotation(PPP2CA, 'http://identifiers.org/pfam/PF00149', 'is')
       Annotation(PPP2CA, 'http://identifiers.org/uniprot/P63330', 'is')
       Annotation(PPP2CA, 'http://identifiers.org/hgnc/HGNC:9299', 'is')
       Annotation(MAP2K1, 'http://identifiers.org/uniprot/Q02750', 'is')
       Annotation(MAP2K1, 'http://identifiers.org/hgnc/HGNC:6840', 'is')
       Annotation(MAPK1, 'http://identifiers.org/pfam/PF00069', 'is')
       Annotation(MAPK1, 'http://identifiers.org/uniprot/P63085', 'is')
       Annotation(MAPK1, 'http://identifiers.org/hgnc/HGNC:6871', 'is')
       Annotation(DUSP6, 'http://identifiers.org/uniprot/Q16828', 'is')
       Annotation(DUSP6, 'http://identifiers.org/hgnc/HGNC:3072', 'is')
