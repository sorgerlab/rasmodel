""" MEK is phosphorylated (2-step catalysis) by BRAF (wt &  v600e) and CRAF 
MEK phosphorylates ERK """

from pysb import *
from pysb.macros import catalyze_state
from pysb.util import alias_model_components


def monomers():
    Monomer('MEK', ['s', 'erk', 'inh', 'state'], {'state': ['up', 'p']})
    Monomer('PP2A', ['mek'])
    Monomer('ERK', ['s', 'sos', 'state'], {'state': ['up', 'p']})
    Monomer('DUSP', ['erk'])

    Parameter('ERK_0', 1e5)  # 1e3
    Parameter('DUSP_0', 1e3)
    Parameter('MEK_0', 1e5)
    Parameter('PP2A_0', 1e3)

    alias_model_components()

    Initial(ERK(s=None, sos=None,  state='up'), ERK_0)
    Initial(DUSP(erk=None), DUSP_0)
    Initial(MEK(s=None, erk=None, inh=None, state='up'), MEK_0)
    Initial(PP2A(mek=None), MEK_0)
    

def by_BRAF_wt():

    Parameter('k_bwf', 1)
    Parameter('k_bwr', 0.1)
    Parameter('k_bwe', 1)

    alias_model_components()

    Rule('BRAF_binds_ERK',
         BRAF(d=1, vem=None, erk=None) % BRAF(d=1) + MEK(s=None, state='up') <>
         MEK(s=2, state='up') % BRAF(d=1, vem=None, erk=2) % BRAF(d=1),
         k_bwf, k_bwr)

    Rule('BRAF_activates_ERK',
         BRAF(d=ANY, vem=None, erk=1) % MEK(s=1, state='up') >>
         BRAF(d=ANY, vem=None, erk=None) + MEK(s=None, state='p'), k_bwe) 


def by_BRAF_mut():

    Parameter('k_bmf', 1)
    Parameter('k_bmr', 0.1)
    Parameter('k_bme', 3)

    alias_model_components()

    catalyze_state(BRAF(vem=None), 'erk', MEK(), 's',
                   'state', 'up', 'p', (k_bmf, k_bmr, k_bme))


def by_CRAF():

    Parameter('k_cwf', 1)
    Parameter('k_cwr', 0.1)
    Parameter('k_cwe', 1)

    alias_model_components()

    Rule('CRAF_binds_ERK',
         CRAF(d=1, vem=None, erk=None) % CRAF(d=1) + MEK(s=None, state='up') <>
         MEK(s=2, state='up') % CRAF(d=1, vem=None, erk=2) % CRAF(d=1),
         k_cwf, k_cwr)

    Rule('CRAF_activates_ERK',
         CRAF(d=ANY, vem=None, erk=1) % MEK(s=1, state='up') >>
         CRAF(d=ANY, vem=None, erk=None) + MEK(s=None, state='p'), k_cwe)


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


def DUSP_phospatase():

    Parameter('k_dspf', 1)
    Parameter('k_dspr', 0.001)
    Parameter('k_dspe', 10)

    alias_model_components()

    catalyze_state(DUSP(), 'erk', ERK(sos=None), 's', 'state', 'p', 'up',
                   (k_dspf, k_dspr, k_dspe))


def ERK_feedback():

    Parameter('k_epsf', 1e-4)
    Parameter('k_epsr', 0.1)
    Parameter('k_epse', 1)

    alias_model_components()

    catalyze_state(ERK(state='p', s=None), 'sos', SOS(ras=None, phos=None),
                   'erk', 'state', 'up', 'p', (k_epsf, k_epsr, k_epse))


def MEK_inhibitor():
    Monomer('Trametinib', ['mek'])

    Parameter('Trametinib_0', 0)

    Parameter('kf_btm', 1)
    Parameter('kr_btm', 1e-5)

    alias_model_components()

    Initial(Trametinib(mek=None), Trametinib_0)

    Rule('bind_Trametinib', Trametinib(mek=None) + MEK(inh=None) <>
         Trametinib(mek=1) % MEK(inh=1), kf_btm, kr_btm)
    

def declare_observables():
    alias_model_components()

    Observable('ERK_P', ERK(state='p'))    
    Observable('MEK_P', MEK(state='p'))
