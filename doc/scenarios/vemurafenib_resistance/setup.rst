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
