MAPK module with MEK inhibition
===============================

Import dependencies

::
   
   from pysb import Monomer, Rule, Parameter, Annotation, ANY, Initial, Observable
   from pysb.macros import bind, _macro_rule, catalyze_state
   from pysb.util import alias_model_components

MEK phosphorylation
::

   def mapk_monomers():
       Monomer('MEK', ['s', 'erk', 'inh', 'state'], {'state': ['up', 'p']})
       Monomer('PP2A', ['mek'])
       Monomer('ERK', ['s', 'sos1', 'state'], {'state': ['up', 'p']})
       Monomer('DUSP', ['erk'])

       Parameter('ERK_0', 1e5)  # 1e3
       Parameter('DUSP_0', 1e3)
       Parameter('MEK_0', 1e5)
       Parameter('PP2A_0', 1e3)

       alias_model_components()

       Initial(ERK(s=None, sos1=None,  state='up'), ERK_0)
       Initial(DUSP(erk=None), DUSP_0)
       Initial(MEK(s=None, erk=None, inh=None, state='up'), MEK_0)
       Initial(PP2A(mek=None), MEK_0)

   def RAF_activates_MEK():

       Parameter('k_bmf', 1)
       Parameter('k_bmr', 0.1)
       Parameter('k_bme', 3)

       alias_model_components()

       catalyze_state(RAF(vem=None), 'erk', MEK(), 's',
		      'state', 'up', 'p', (k_bmf, k_bmr, k_bme))

   def MEK_phosphorylates_ERK():

       Parameter('k_mef', 1)
       Parameter('k_mer', 0.1)
       Parameter('k_mee', 10)

       alias_model_components()

       catalyze_state(MEK(s=None, state='p', inh=None), 'erk', ERK(), 's',
		      'state', 'up', 'p', (k_mef, k_mer, k_mee))

   def PP2A_phosphatase():

       Parameter('k_pp2f', 1)
       Parameter('k_pp2r', 0.001)
       Parameter('k_pp2e', 10)

       alias_model_components()

       catalyze_state(PP2A(), 'mek', MEK(erk=None), 's', 'state', 'p', 'up',
		      (k_pp2f, k_pp2r, k_pp2e))		      


   def DUSP_phosphatase():

       Parameter('k_dspf', 1)
       Parameter('k_dspr', 0.001)
       Parameter('k_dspe', 10)

       alias_model_components()

       catalyze_state(DUSP(), 'erk', ERK(sos1=None), 's', 'state', 'p', 'up',
			 (k_dspf, k_dspr, k_dspe))


   def ERK_feedback():

       Parameter('k_epsf', 1e-4)
       Parameter('k_epsr', 0.1)
       Parameter('k_epse', 1)

       alias_model_components()

       catalyze_state(ERK(state='p', s=None), 'sos1', SOS1(ras=None, phos=None),
		      'mapk1', 'state', 'up', 'p', (k_epsf, k_epsr, k_epse))

   def MEK_inhibitor():
       Monomer('Trametinib', ['mek'])

       Parameter('Trametinib_0', 0)

       Parameter('kf_btm', 1)
       Parameter('kr_btm', 1e-5)

       alias_model_components()

       Initial(Trametinib(mek=None), Trametinib_0)

       Rule('bind_Trametinib', Trametinib(mek=None) + MEK(inh=None) <>
	    Trametinib(mek=1) % MEK(inh=1), kf_btm, kr_btm)
