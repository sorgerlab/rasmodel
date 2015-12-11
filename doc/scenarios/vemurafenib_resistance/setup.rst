Vemurafenib resistance
======================

::
   
    from pysb import *

    Model()

    import rasmodel.components.egfr as egfr
    egfr.egfr_minimal()

    import rasmodel.components.ras as ras
    ras.ras_minimal()

    import rasmodel.components.raf as raf
    raf.raf_monomers()
    raf.ras_activates_raf()
    raf.vemurafenib_monomers()
    raf.vemurafenib_binds_raf()
    raf.mek_phosphorylation()

    import rasmodel.components.mapk as mapk
    mapk.erk_dynamics()

Model extensions
----------------
Negative feedback by MAPK1-dependent phosphorylation of SOS is an important mechanism that underlies
(a) RAS independence of BRAFV600E mutants and
(b) the 'imperfect adaptation' resulting in Vemurafenib resistance

The following extensions were made to the model assembled from components
1. A phosphorylation state is added to Monomer SOS
::
   a = model.monomers['SOS1']
   a.sites.append('mapk1')
   a.sites.append('state')
   a.site_states['state'] = ['up', 'p']

   b = model.monomers['MAPK1']
   b.sites.append('sos1')

2. SOS can bind Grb2 and activate RAS only when it is unphosphorylated
::
   
   r = model.rules['GRB2_EGFR_SOS1_bind']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   SOS1(grb2=None, state='up')
   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   SOS1(grb2=1, state='up')
   
   r = model.rules['SOS1_GRB2_RAS_bind']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[0] =\
   SOS1(grb2=ANY, ras=None, state='up')
   r.product_pattern.complex_patterns[0].monomer_patterns[0] =\
   SOS1(grb2=ANY, ras=1, state='up')

   r = model.rules['GRB2_EGFR_SOS1_dissociate']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[1] =\
   SOS1(grb2=1, state='up')
   r.product_pattern.complex_patterns[1].monomer_patterns[0] =\
   SOS1(grb2=None, state='up')

   raf.erk_feedback()
    

    
Model reduction
---------------
The following constraints are put on interactions for reducing combinatorial complexity in the model:-


1. EGFR dimers dissociate only when bound to EGF
::

   r = model.rules['EGFR_EGF_EGFR_EGF_dissociate']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[0] =\
   EGFR(egfr=1, egf=ANY)
   r.reactant_pattern.complex_patterns[0].monomer_patterns[1] =\
   EGFR(egfr=1, egf=ANY)

   r.product_pattern.complex_patterns[0].monomer_patterns[0] =\
   EGFR(egfr=None, egf=ANY)
   r.product_pattern.complex_patterns[1].monomer_patterns[0] =\
   EGFR(egfr=None, egf=ANY)

2. EGFR binds Grb2 only when when EGFR is in the active dimer form.
::
   r = model.rules['GRB2_EGFR_bind']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   EGFR(Y='p', grb2=None, egfr=ANY)

   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   EGFR(Y='p', grb2=1, egfr=ANY)

3. Grb2 and EGFR dissociate only wwhen EGFR is in the active dimer form.
::
   r = model.rules['GRB2_EGFR_dissociate']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[1] =\
   EGFR(grb2=1, Y='p', egfr=ANY)
   r.product_pattern.complex_patterns[1].monomer_patterns[0] =\
   EGFR(grb2=None, Y='p', egfr=ANY)

4. Grb2 and SOS1 dissociate when Grb2 is bound to bound to EGFR
::
   r = model.rules['GRB2_EGFR_SOS1_dissociate']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[0] =\
   GRB2(sos1=1, egfr=ANY)
   r.product_pattern.complex_patterns[0].monomer_patterns[0] =\
   GRB2(sos1=None, egfr=ANY)

5. SOS1 binds RAS when it is not bound to its substrate RAF
::

   r = model.rules['SOS1_GRB2_RAS_bind']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   RAS(sos1=None, raf=None)
   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   RAS(sos1=1, raf=None)

   
6. SOS1 and RAS dissociate only when SOS1 is bound to Grb2
::
   r = model.rules['SOS1_GRB2_RAS_dissociate']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[0] =\
   SOS1(ras=1, grb2=ANY)
   r.product_pattern.complex_patterns[0].monomer_patterns[0] =\
   SOS1(ras=None, grb2=ANY)

   r.reactant_pattern.complex_patterns[0].monomer_patterns[1] =\
   RAS(sos1=1, raf=None)
   r.product_pattern.complex_patterns[1].monomer_patterns[0] =\
   RAS(sos1=None, raf=None)
   
7. PP2A dephosphorylates MAP2K1 only when MAP2K1 is not bound to RAF or MAPK1
::
   r = model.rules['PPP2CA_dephos_bind_map2k1_S218_1']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   MAP2K1(S218='p', ppp2ca=None, raf=None, mapk1=None)
   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   MAP2K1(S218='p', ppp2ca=1, raf=None, mapk1=None)

   r = model.rules['PPP2CA_dephos_bind_map2k1_S222_1']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   MAP2K1(S222='p', ppp2ca=None, raf=None, mapk1=None)
   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   MAP2K1(S222='p', ppp2ca=1, raf=None, mapk1=None)

8. DUSP6 dephosphorylates MAPK1 only when MAPK is not bound to  MAP2K1
::
   r = model.rules['DUSP6_dephos_bind_MAPK1_Y187_1']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   MAPK1(Y187='p', dusp6=None, sos1=None, map2k1=None)
   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   MAPK1(Y187='p', dusp6=1, sos1=None, map2k1=None)

   r = model.rules['DUSP6_dephos_bind_MAPK1_T185_1']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   MAPK1(T185='p', dusp6=None, sos1=None, map2k1=None)
   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   MAPK1(T185='p', dusp6=1, sos1=None, map2k1=None)

9. MAPK21 phosphorylates MAPK1 only when PP2A or DUSP6 are not bound to MAPK1
::
   r = model.rules['MAP2K1_phospho_bind_MAPK1_Y187_1']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   MAPK1(Y187='u', map2k1=None, dusp6=None)
   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   MAPK1(Y187='u', map2k1=1, dusp6=None)

   r.reactant_pattern.complex_patterns[0].monomer_patterns[0] =\
   MAP2K1(S218='p', S222='p', mapk1=None, ppp2ca=None)
   r.product_pattern.complex_patterns[0].monomer_patterns[0] =\
   MAP2K1(S218='p', S222='p', mapk1=1, ppp2ca=None)
   
   r = model.rules['MAP2K1_phospho_bind_MAPK1_T185_1']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   MAPK1(T185='u', map2k1=None, dusp6=None)
   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   MAPK1(T185='u', map2k1=1, dusp6=None)

   r.reactant_pattern.complex_patterns[0].monomer_patterns[0] =\
   MAP2K1(S218='p', S222='p', mapk1=None, ppp2ca=None)
   r.product_pattern.complex_patterns[0].monomer_patterns[0] =\
   MAP2K1(S218='p', S222='p', mapk1=1, ppp2ca=None)
