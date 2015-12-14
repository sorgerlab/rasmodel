Setup for Vemurafenib resistance scenario
=========================================

::

   from pysb import *

   Model()

   import rasmodel.components.egfr as egfr
   egfr.egfr_minimal()

A phosphorylation state is added to SOS1
::
   a = model.monomers['SOS1']
   a.sites.append('mapk1')
   a.sites.append('state')
   a.sites.append('phos')
   a.site_states['state'] = ['up', 'p']
   

Minimal model of RAS activation
::
   import rasmodel.components.ras as ras
   ras.ras_minimal()

Model reduction
---------------

SOS can bind Grb2 only when it is unphosphorylated
::
   
   r = model.rules['GRB2_EGFR_SOS1_bind']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   SOS1(grb2=None, ras=None, state='up')
   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   SOS1(grb2=1, ras=None, state='up')

   r = model.rules['GRB2_EGFR_SOS1_dissociate']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[1] = SOS1(grb2=1, ras=None, state='up')
   r.product_pattern.complex_patterns[1].monomer_patterns[0] = SOS1(grb2=None, ras=None, state='up') 
   
SOS1 binds RAS when it is not bound to its substrate RAF or bound to GTP
::

   r = model.rules['SOS1_GRB2_RAS_bind']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   RAS(sos1=None, gtp=None, raf=None)
   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   RAS(sos1=1, gtp=None, raf=None)

   r.reactant_pattern.complex_patterns[0].monomer_patterns[0] =\
   SOS1(grb2=ANY, ras=None, state='up')
   r.product_pattern.complex_patterns[0].monomer_patterns[0] =\
   SOS1(grb2=ANY, ras=1, state='up')

   
SOS1 and RAS dissociate only when SOS1 is bound to Grb2
::
   r = model.rules['SOS1_GRB2_RAS_dissociate']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[0] =\
   SOS1(ras=1, grb2=ANY, state='up')
   r.product_pattern.complex_patterns[0].monomer_patterns[0] =\
   SOS1(ras=None, grb2=ANY, state='up')

   r.reactant_pattern.complex_patterns[0].monomer_patterns[1] =\
   RAS(sos1=1, gtp=None, raf=None)
   r.product_pattern.complex_patterns[1].monomer_patterns[0] =\
   RAS(sos1=None, gtp=None, raf=None)

Dissociation of RAS:GTP and SOS1
::

   r = model.rules
   r.add(Rule('RASGTP_SOS1_dissociate', RAS(gtp=ANY, sos1=1, raf=None) % SOS1(ras=1, grb2=ANY) >>\
   RAS(gtp=ANY, sos1=None, raf=None) + SOS1(ras=None, grb2=ANY), Parameter('ke4', 50)))
   
   import rasmodel.components.raf as raf
   raf.RAF_monomers()

   raf.RAF_dynamics()

   raf.SOS1_dephosphorylation()

   import rasmodel.components.mapk as mapk
   mapk.mapk_monomers()

   mapk.RAF_activates_MEK()

   mapk.MEK_phosphorylates_ERK()

   mapk.PP2A_phosphatase()

   mapk.DUSP_phosphatase()

   mapk.ERK_feedback()
   mapk.MEK_inhibitor()

   mapk.observables()


   
Model reduction
---------------
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


Initial conditions
------------------
::

   Initial(EGF(egfr=None), Parameter('EGF_0', 1e3))
   Initial(EGFR(egf=None, egfr=None, grb2=None, Y='u'), Parameter('EGFR_0', 1e5))
   Initial(GRB2(egfr=None, sos1=None), Parameter('GRB2_0', 1e5))
   Initial(SOS1(grb2=None, ras=None, mapk1=None, phos=None, state='up'), Parameter('SOS1_0', 1e3))
   Initial(GTP(ras=None), Parameter('GTP_0', 1e9))
   Initial(RAS(sos1=None, gtp=None, raf=None), Parameter('RAS_0', 1e5))

Parameters
----------
::

   kf_ee_bind_1.value =  1  # EGF bind EGFR
   kr_ee_bind_1.value =  0.1  # EGF: EGFR dissociate
   kf_ee_bind_2.value =  1  # EGFR binds EGFR
   kr_ee_bind_2.value =  0.1  # EGFR:EGFR dissociate
   kf_ge_bind_1.value = 1   # Grb2 binds EGFR
   kr_ge_bind_1.value = 0.1   # GRb2:EGFR dissociate
   kf_gs_bind_1.value =  1  #  Grb2 binds SOS1
   kr_gs_bind_1.value =  0.1  # Grb2:SOS1 complex dissociates
   kf_ee_transphos_1.value =  1 # transphosphorylation of EGFR in dimer

   kf_sr_bind_1.value = 1
   kr_sr_bind_1.value = 1e-3
   kf_rg_bind_1.value = 1
   kr_rg_bind_1.value = 0.5 
