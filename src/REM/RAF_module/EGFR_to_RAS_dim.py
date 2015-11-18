""" Straw model of EGFR activation """

from pysb import *
from pysb.util import alias_model_components
from pysb.macros import catalyze_state


def monomers():
    # Declare monomer species in the model
    # ------------------------------------
    Monomer('EGF', ['egfr'])
    Monomer('EGFR', ['egf', 'd', 'grb2'])
    Monomer('Grb2', ['egfr', 'sos'])
    Monomer('SOS', ['grb2', 'ras', 'erk', 'phos', 'state'],
            {'state': ['up', 'p']})
    Monomer('KRAS', ['raf', 'sos', 'state'], {'state': ['gdp', 'gtp']})

    # IC values
    # ---------
    Parameter('EGF_0', 1e3)
    Parameter('EGFR_0', 1e5)
    Parameter('Grb2_0', 1e5)
    Parameter('SOS_0', 1e3)
    Parameter('KRAS_0', 2e5)

    alias_model_components()

    # Initial conditions
    # ------------------
    Initial(EGF(egfr=None), EGF_0)
    Initial(EGFR(egf=None, d=None, grb2=None), EGFR_0)
    Initial(Grb2(egfr=None, sos=None), Grb2_0)
    Initial(SOS(grb2=None, ras=None, erk=None, phos=None, state='up'), SOS_0)
    Initial(KRAS(sos=None, raf=None, state='gdp'), KRAS_0)


def KRAS_activation():

    # Rate constants
    # -------------
    Parameter('kf1', 1)
    Parameter('kr1', 0.1)
    Parameter('kf2', 1)
    Parameter('kr2', 0.1)
    Parameter('kf3', 1)
    Parameter('kr3', 0.1)
    Parameter('kf4', 1)
    Parameter('kr4', 0.001)
    Parameter('ke4', 50)
    Parameter('kf5', 1)
    Parameter('kr5', 1e-9)
    Parameter('k_dimf', 1)
    Parameter('k_dimr', 0.1)

    alias_model_components()

    # Rules
    # -----

    # EGF binds EGFR
    Rule('EGF_binds_EGFR',
         EGF(egfr=None) + EGFR(egf=None) <>
         EGF(egfr=1) % EGFR(egf=1), kf1, kr1)

    # The EGF-EGFR complex binds another EGF-EGFR complex
    Rule('EGFR_dimerization',
         EGFR(egf=ANY, d=None) + EGFR(egf=ANY, d=None) <>
         EGFR(egf=ANY, d=1) % EGFR(egf=ANY, d=1), k_dimf, k_dimr)

    # EGFR, bound to EGFR, binds Grb2
    Rule('EGFR_binds_Grb2',
         EGFR(d=ANY, grb2=None) + Grb2(egfr=None) <>
         EGFR(d=ANY, grb2=1) % Grb2(egfr=1), kf2, kr2)

    # Grb2, bound to EGFR, binds SOS
    Rule('Grb2_binds_SOS',
         Grb2(egfr=ANY, sos=None) + SOS(grb2=None, state='up') <>
         Grb2(egfr=ANY, sos=1) % SOS(grb2=1, state='up'), kf3, kr3)

    # SOS activates KRAS
    catalyze_state(SOS(grb2=ANY, state='up'), 'ras',
                   KRAS(raf=None), 'sos', 'state', 'gdp', 'gtp',
                   (kf4, kr4, ke4))

    # KRAS deactivates itself
    # Making this step reversible increased combinatorial complexity manifold
    Rule('KRAS_inactivation',
         KRAS(state='gtp') >>
         KRAS(state='gdp'),
         kf5)


def SOS_dephosphorylation():

    Monomer('SOS_phos', ['sos'])

    Parameter('SOS_phos_0', 100)
    Parameter('k_spf', 1)
    Parameter('k_spr', 0.1)
    Parameter('k_spe', 1e-4)

    alias_model_components()

    Initial(SOS_phos(sos=None), SOS_phos_0)

    catalyze_state(SOS_phos(), 'sos', SOS(ras=None, erk=None), 'phos',
                   'state', 'p', 'up', (k_spf, k_spr, k_spe))


def declare_observables():

    # Observables
    # -----------
    Observable('active_KRAS', KRAS(state='gtp'))
    Observable('active_SOS', SOS(state='up'))
