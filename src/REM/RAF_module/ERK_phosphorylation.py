""" ERK is phosphorylated (2-step catalysis) by BRAF (wt &  v600e) and CRAF """

from pysb import *
from pysb.macros import catalyze_state
from pysb.util import alias_model_components


def monomers():
    Monomer('ERK', ['s', 'state'], {'state': ['up', 'p']})
    Monomer('DUSP', ['erk'])

    Parameter('ERK_0', 1e3)
    Parameter('DUSP_0', 1e2)

    alias_model_components()

    Initial(ERK(s=None, state='up'), ERK_0)
    Initial(DUSP(erk=None), DUSP_0)
    

def by_BRAF_wt():

    Parameter('k_bwf', 1)
    Parameter('k_bwr', 0.1)
    Parameter('k_bwe', 1)

    alias_model_components()

    Rule('BRAF_binds_ERK',
         BRAF(d=1, vem=None, erk=None) % BRAF(d=1) + ERK(s=None, state='up') <>
         ERK(s=2, state='up') % BRAF(d=1, vem=None, erk=2) % BRAF(d=1),
         k_bwf, k_bwr)

    Rule('BRAF_activates_ERK',
         BRAF(d=ANY, vem=None, erk=1) % ERK(s=1, state='up') >>
         BRAF(d=ANY, vem=None, erk=None) + ERK(s=None, state='p'), k_bwe) 


def by_BRAF_mut():

    Parameter('k_bmf', 1)
    Parameter('k_bmr', 0.1)
    Parameter('k_bme', 3)

    alias_model_components()

    catalyze_state(BRAF(vem=None), 'erk', ERK(), 's',
                   'state', 'up', 'p', (k_bmf, k_bmr, k_bme))


def by_CRAF():

    Parameter('k_cwf', 1)
    Parameter('k_cwr', 0.1)
    Parameter('k_cwe', 1)

    alias_model_components()

    Rule('CRAF_binds_ERK',
         CRAF(d=1, vem=None, erk=None) % CRAF(d=1) + ERK(s=None, state='up') <>
         ERK(s=2, state='up') % CRAF(d=1, vem=None, erk=2) % CRAF(d=1),
         k_cwf, k_cwr)

    Rule('CRAF_activates_ERK',
         CRAF(d=ANY, vem=None, erk=1) % ERK(s=1, state='up') >>
         CRAF(d=ANY, vem=None, erk=None) + ERK(s=None, state='p'), k_cwe)


def DUSP_phospatase():

    Parameter('k_dspf', 1)
    Parameter('k_dspr', 0.1)
    Parameter('k_dspe', 1)

    alias_model_components()

    catalyze_state(DUSP(), 'erk', ERK(), 's', 'state', 'p', 'up',
                   (k_dspf, k_dspr, k_dspe))


def declare_observables():
    alias_model_components()

    Observable('ERK_P', ERK(state='p'))    
