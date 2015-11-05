""" 
Overview
========

PySB implimentation of the Chen 2009 model, originally written in MATLAB's simbiology toolbox, and published as:

William Chen, Birgit Schoeberl et al.(2009) Molecular Systems Biology 5:239; doi:10.1038/msb.2008.74

"""

from pysb import*
from pysb.util import alias_model_components
import pysb.core
import pysb.macros
from pysb.macros import bind_table
from pysb.macros import catalyze
from pysb.macros import catalyze_state


# Monomers in the model
# =====================
def receptor_monomers():
    """ Declare the 4 ErbB receptors in the model.
        
        A description of sites on the receptors is given below
        ======================================================
        'lig' is the site to bind ligands EGF or HRG
        'd' is the site at which monomers bind to each other
        'atp' is the site for ATP binding
        'gap' is the site to bind GAP
        'gs' is the site to bind Shc or Grb2
        'rtk' is the site to bind RTK phosphatase
        'cpp' is the site to bind endosomol transport protein cPP
        'state' denotes the phoshorylation status of ErbB, with 'up' denoting unphosphorylated, and 'p' denting phoshorylated'
        'comp' denotes the compartment in which the receptor is present, plamsa membrane ('pm') or endosome ('endo')
        """
    
    receptor_names = ['ErbB1', 'ErbB2', 'ErbB3', 'ErbB4']
    
    for erb in receptor_names:
        Monomer(erb, ['lig', 'd', 'atp', 'state', 'gap', 'gs', 'rtk', 'cpp', 'comp'],
                {'state':['up', 'p','full_act', 'pm', 'endo'], 'comp': ['pm', 'endo']})
    
    # ERK is defined here since it is used in both MAPK and AKT pathways
    Monomer('ERK', ['mek', 'pase3', 'gab1','sos', 'state'], {'state':['up','p', 'pp']})

    alias_model_components()

def ligand_monomers():
    """ Declare ligands in the model. 
        
        A description of sites on the ligands is given below
        ====================================================
        'rec' denotes the site at which ligands binds receptors
        'comp' denotes the compartment in which the receptor is present, plamsa membrane ('pm') or endosome ('endo')
    """

    for ligand in 'EGF', 'HRG':
        Monomer(ligand, ['rec', 'comp'], {'comp': ['pm', 'endo']})
    
    alias_model_components()
    
def adapter_monomers():
    """ Declare components that directly bind ErbB receptors, such as ATP, GAP, Grb2, Shc, RTK phosphatase, and cPP.
        
        A description of sites on the monomers is given below
        ========================================================
        'erb' is a site on ATP, Grb2, Shc, RTK, and cPP used to bind receptors
        'rec' is a site on GAP to bind receptors. ### Consider change to 'erb' for consistency
        'gab1' is a site on ATP to bind Gab1
        'shc' is a site of Grb2 to bind Shc
        'sos' is a site on Grb2 to bind SOS
        'grb' is a site on SOS to binf Grb2
        'ras' and 'erk' are sites on SOS to bind RAS and ERK respectively.
        'state' denotes the phoshorylation status of the species, with 'up' denoting unphosphorylated, and 'p' denoting phoshorylated'
        'comp' denotes the compartment in which the species is present, plamsa membrane ('pm') or endosome ('endo')
    """
    
    Monomer('ATP', ['erb','gab1'])
    Monomer('RTK', ['erb']) # RTK phosphatase
    Monomer('GAP', ['rec']) # Currenlty unknown intermediate
    Monomer('Grb2', ['shc', 'erb', 'sos', 'gab1'])
    Monomer('SOS', ['grb', 'ras', 'erk','state'],{'state':['up','p']})
    Monomer('SHC', ['erb','grb', 'state'], {'state':['up','p']})
    Monomer('cPP', ['erb', 'comp'], {'comp':['pm', 'endo']})
    Monomer('RAS', ['sos', 'pi3k', 'raf', 'state'], {'state':['gdp', 'gtp', 'active_gtp']})
    
    alias_model_components()

########################################################
def ErbB1_priming():
    " v828 in sbml model"
    " Unlike other Erb monomers, ErbB1 can participate in reactions only after binding with ATP "
    
    # Initial amounts
    Parameter('ErbB1_0',	1.08e6)     # c531
    Parameter('ATP_0', 1.2e9)           # c105
    # Rate constants
    Parameter('k122_priming', 1.8704e-8)  # k122
    Parameter('kd122_priming', 1)       # kd122
    
    alias_model_components()
    
    # Initial conditions
    # ===================
    Initial(ErbB1(lig=None, d=None, atp=None, gap=None, gs=None, rtk=None, cpp=None, state='up', comp ='pm'), ErbB1_0)
    Initial(ATP(erb=None, gab1=None), ATP_0)

    # Rules
    # =============
    # ErbB1 + ATP <-> ErbB1:ATP
    Rule('ErbB1_atp', ErbB1(lig=None,d=None,atp=None, comp='pm') + ATP(erb=None,gab1=None) <>
         ErbB1(lig=None,d=None,atp=1, comp='pm') % ATP(erb=1,gab1=None), k122_priming, kd122_priming)

########################################################
def ligand_binding():
    " EGF binds ErbB1 (v1), HRG bind ErbB3 and ErbB4 (v784-785)"
    " Rate constants differ based on receptor and compartment "
    
    # Initial amounts
    Parameter('ErbB2_0', 4.62e5)        # c141
    Parameter('ErbB3_0', 6230)          # 140
    Parameter('ErbB4_0', 794)           # c143
    Parameter('EGF_0', 5e-9)            # c1
    Parameter('HRG_0', 0)               # c514
    Parameter('dummy_init', 0)
    # Rate constants
    Parameter('k1', 1e+7)
    Parameter('kd1', 0.0033)
    Parameter('k119', 1e+7)
    Parameter('kd119', 0.0103115)
    Parameter('k10b',0.05426)
    Parameter('kd10',0.011)
    
    alias_model_components()
    
    # Initial conditions
    # ==================
    for erb in receptors[1:]:
        Initial(erb(lig=None, d=None, atp=None, gap=None, gs=None, rtk=None, cpp=None, state='up', comp ='pm'), eval(erb.name + '_0'))
  
    for ligand in EGF, HRG:
        Initial(ligand(rec=None, comp ='pm'), eval(ligand.name + '_0'))

    " Initialize endo|HRG and endo| HRG:ErbB4 since there are no reactions leading to these species "
    Initial(HRG(rec=None, comp ='endo'), dummy_init)
    Initial(HRG(rec=1, comp='endo')% ErbB4(lig=1, d=None,atp=None, gap=None, gs=None, rtk=None, cpp=None, state='up', comp ='endo'), dummy_init)

    global comprtmnts
    comprtmnts = ['pm', 'endo']
    
    # ErbB1:ATP + EGF -> EGF:ErbB1:ATP (plasma membrane)
    Rule('ErbB1_EGF_pm', ErbB1(lig=None,d=None,atp=ANY, state='up', comp='pm') + EGF(rec=None, comp='pm') <>
     ErbB1(lig=1,d=None,atp=ANY, state='up', comp='pm') % EGF(rec=1, comp='pm'), k1, kd1)
    
    # ErbB1:ATP + EGF -> EGF:ErbB1:ATP (plasma membrane)
    Rule('ErbB1_EGF_endo', ErbB1(lig=None,d=None,atp=ANY, state='up', comp='endo') + EGF(rec=None, comp='endo') <>
         ErbB1(lig=1,d=None,atp=ANY, state='up', comp='endo') % EGF(rec=1, comp='endo'), k10b, kd10)
        
    # ErbB3 + HRG -> ErbB3:HRG (plasma membrane)
    Rule('ErbB3_HRG_pm', ErbB3(lig=None,d=None,atp=None,state='up', comp='pm') + HRG(rec=None, comp='pm') <>
              ErbB3(lig=1,d=None,atp=None, state='up', comp='pm') % HRG(rec=1, comp='pm'), k119, kd119)
         
    # ErbB3 + HRG -> ErbB3:HRG (endosome)
    Rule('ErbB3_HRG_endo', ErbB3(lig=None,d=None,atp=None,state='up', comp='endo') + HRG(rec=None, comp='endo') <>
              ErbB3(lig=1,d=None,atp=None, state='up', comp='endo') % HRG(rec=1, comp='endo'), k10b, kd10)
         
    # ErbB4 + HRG -> ErbB4:HRG (this reaction takes place only in pm compartment as per Chen 2009)
    Rule('ErbB4_HRG_pm', ErbB4(lig=None,d=None,atp=None, state='up', comp='pm') + HRG(rec=None, comp='pm') <>
          ErbB4(lig=1,d=None,atp=None, state='up', comp='pm') % HRG(rec=1, comp='pm'), k119, kd119)


########################################################
def transport():
    " Free ErbB monomers (v164,v176-178), phosphorylated dimers (v165-175,179-193) and cPP(v211) can \
    shuttle between plasma membrane and endosome "
    
    Parameter('k15', 1.667e-8)
    Parameter('kd15', 0)
    Parameter('k7', 5e-05)
    Parameter('kd7', 1.38e-4)
    Parameter('k6', 0.013)
    Parameter('kd6', 5e-05)
    Parameter('k6b', 0)
    Parameter('kd6b', 0)

    alias_model_components()
    
    global receptors
    receptors = [ErbB1, ErbB2, ErbB3, ErbB4]
    
    # v164 ErbB1:ATP <-> ErbB1:ATP
    Rule('Transport_ErbB1_ATP', ErbB1(d=None, lig=None, atp=ANY, state='up', comp='pm') <>
         ErbB1(d=None, lig=None, atp=ANY, state='up', comp='endo'), k6, kd6)
        
    " V176-178 and v181-184 have rate constants set to zero in both directions "
    # ErbB[2-4] monomers transported bw plasma membrane and endosome
    for erb in receptors[1:]:
        Rule('transport_' + erb.name, erb(lig=None, d=None, atp=None, state='up', comp='pm') <>
              erb(lig=None, d=None, atp=None, state='up', comp='endo'), k6b, kd6b)

    " v165 - 175 --> Transport of ErbB1 homodimers. Rules have to be carefully specified to prevent transport of dimers bound to Gab1, Erk, or SOS~P "
    " Specifying EGF in the rules is important to avoid creating species where EGF and bound ErbB1 are in different compartments "
    " Discuss refactoring "
    # 2(EGF:ErbB1)
    Rule('ErbB11_pm_to_endo', EGF(rec=2, comp='pm') % EGF(rec=3, comp='pm') %
         ErbB1(lig=2, d=1, rtk=None, cpp=None, gap=None, state='p', comp='pm')
         % ErbB1(lig=3, d=1, rtk=None, state='p', cpp=None, gap=None, comp='pm') <>
         EGF(rec=2, comp='endo') % EGF(rec=3, comp='endo') %
         ErbB1(lig=2, d=1, rtk=None, cpp=None, gap=None, state='p', comp='endo') %
         ErbB1(lig=3, d=1, rtk=None, state='p', cpp=None, gap=None, comp='endo'), k6, kd6)
    # 2(EGF:ErbB1):GAP
    Rule('ErbB11_gap_pm_to_endo', EGF(rec=2, comp='pm') % EGF(rec=3, comp='pm') %
         ErbB1(lig=2, d=1, atp=None, rtk=None, cpp=None, gap=ANY, gs=None, state='p', comp='pm')%
         ErbB1(lig=3, d=1, rtk=None, state='p', cpp=None, atp=None, gap=None, comp='pm') <>
          EGF(rec=2, comp='endo') % EGF(rec=3, comp='endo') %
         ErbB1(lig=2, d=1, rtk=None, cpp=None, atp=None, gap=ANY, gs=None,state='p', comp='endo') %
          ErbB1(lig=3, d=1, rtk=None, state='p', atp=None,cpp=None, gap=None, comp='endo') , k6, kd6)
    # 2(EGF:ErbB1):GAP:SHC['p' or 'up']
    Rule('ErbB11_gapshc_pm_to_endo', EGF(rec=2, comp='pm') % EGF(rec=3, comp='pm') %
         ErbB1(lig=2, d=1, atp=None, rtk=None, cpp=None, gap=ANY, gs=ANY, state='p', comp='pm')%
         ErbB1(lig=3, d=1, rtk=None, state='p', cpp=None, atp=None, gap=None, comp='pm') % SHC(erb=ANY,grb=None) <>
          EGF(rec=2, comp='endo') % EGF(rec=3, comp='endo') %
         ErbB1(lig=2, d=1, rtk=None, cpp=None, atp=None, gap=ANY, gs=ANY,state='p', comp='endo') %
          ErbB1(lig=3, d=1, rtk=None, state='p', atp=None,cpp=None, gap=None, comp='endo') % SHC(erb=ANY, grb=None) , k6, kd6)
    # 2(EGF:ErbB1):GAP:Grb2
    Rule('ErbB11_gapgrb_pm_to_endo', EGF(rec=2, comp='pm') % EGF(rec=3, comp='pm') %
         ErbB1(lig=2, d=1, atp=None, rtk=None, cpp=None, gap=ANY, gs=ANY, state='p', comp='pm')%
         ErbB1(lig=3, d=1, rtk=None, state='p', cpp=None, atp=None, gap=None, comp='pm') % Grb2(sos=None, gab1=None) <>
          EGF(rec=2, comp='endo') % EGF(rec=3, comp='endo') %
         ErbB1(lig=2, d=1, rtk=None, cpp=None, atp=None, gap=ANY, gs=ANY,state='p', comp='endo') %
          ErbB1(lig=3, d=1, rtk=None, state='p', atp=None,cpp=None, gap=None, comp='endo') % Grb2(sos=None, gab1=None) , k6, kd6)
    # 2(EGF:ErbB1):GAP:Grb2:SOS['up']
    Rule('ErbB11_gapsos_pm_to_endo', EGF(rec=2, comp='pm') % EGF(rec=3, comp='pm') %
         ErbB1(lig=2, d=1, atp=None, rtk=None, cpp=None, gap=ANY, gs=ANY, state='p', comp='pm')%
         ErbB1(lig=3, d=1, rtk=None, state='p', cpp=None, atp=None, gap=None, comp='pm') % SOS(erk=None, state='up') <>
          EGF(rec=2, comp='endo') % EGF(rec=3, comp='endo') %
         ErbB1(lig=2, d=1, rtk=None, cpp=None, atp=None, gap=ANY, gs=ANY,state='p', comp='endo') %
          ErbB1(lig=3, d=1, rtk=None, state='p', atp=None,cpp=None, gap=None, comp='endo') % SOS(erk=None, state='up') , k6, kd6)

    " Transport of ErbB1 heterodimers that are not bound to any other monomers (v185-187) "
    # ErbB1~P:ErbB[2-4]~P
    for erb in receptors[1:]:
             Rule('ErbB1_pm_to_endo'+erb.name, ErbB1(d=1, rtk=None, cpp=None, gap=None, state='p', comp='pm') %
                  erb(d=1, rtk=None, state='p', cpp=None,gap=None, comp='pm') <>
                  ErbB1(d=1, rtk=None, cpp=None, state='p', gap=None, comp='endo') %
                  erb(d=1, rtk=None, state='p', cpp=None, gap=None, comp='endo'), k7, kd7)

    " Transport of ErbB2 heterodimers that are not bound to any other monomers (v188-190) "
    # ErbB1~P:ErbB[2-4]~P
    for erb in receptors[1:]:
        Rule('ErbB2_pm_to_endo'+erb.name, ErbB2(d=1, rtk=None,atp=None, state='p', gap=None, comp='pm') %
             erb(d=1, atp=None, rtk=None, state='p', gap=None, comp='pm') <>
             ErbB2(d=1, atp=None, rtk=None,state='p', gap=None, comp='endo') %
             erb(d=1, atp=None, rtk=None, state='p', gap=None, comp='endo'), k7, kd7)

    " Transport of ErbB2 dimers bound to GAP "
    # v191-193 ErbB2:ErbB[2-4]:GAP
    Rule('transport_ErbB2_gap', ErbB2(d=1, rtk=None, atp=None, state='p', gap=ANY, gs=None, comp='pm') %
     ErbB2(d=1, rtk=None, atp=None, state='p', gap=None, comp='pm')
     <> ErbB2(d=1, rtk=None, atp=None, state='p', gap=ANY, gs=None, comp='endo') %
     ErbB2(d=1, rtk=None, atp=None, state='p', gap=None, comp='endo') , k7, kd7)

    " Transport of cPP "
    # v211  cPP(endo) <-> cPP(pm) , c9 <-> c12 "
    Rule('transport_cPP', cPP(erb=None, comp='endo') <> cPP(erb=None, comp='pm'), k15, kd15)

    " Transport of ErbB2 dimers bound to Shc. Seperate rate constants based on type of ErbB and Shc phosphorylation status "
    # v180 ErbB2:ErbB4:Shc
    Rule('transport_ErbB24_shc', ErbB2(d=1, rtk=None, atp=None, state='p', gap=ANY, comp='pm') %
         ErbB4(d=1, rtk=None, atp=None, state='p', gap=None, comp='pm') % SHC(grb=None, state='up') <>
         ErbB2(d=1, rtk=None, atp=None, state='p', gap=ANY, comp='endo') %
         ErbB4(d=1, rtk=None, atp=None, state='p', gap=None, comp='endo') % SHC(grb=None, state='up'), k6, kd6)
    # v184 ErbB2:ErbB4:Shc~P
    Rule('transport_ErbB24_shcp', ErbB2(d=1, rtk=None, atp=None, state='p', gap=ANY, comp='pm') %
         ErbB4(d=1, rtk=None, atp=None, state='p', gap=None, comp='pm') % SHC(grb=None, state='p') <>
         ErbB2(d=1, rtk=None, atp=None, state='p', gap=ANY, comp='endo') %
         ErbB4(d=1, rtk=None, atp=None, state='p', gap=None, comp='endo') % SHC(grb=None, state='p'), k6b, kd6b)
         
     # v179, v181 ErbB2:ErbB3:Shc['p' or 'up']
    Rule('transport_ErbB23_shc', ErbB2(d=1, rtk=None, atp=None, state='p', gap=ANY, comp='pm') %
         ErbB3(d=1, rtk=None, atp=None, state='p', gap=None, comp='pm') % SHC(grb=None) <>
         ErbB2(d=1, rtk=None, atp=None, state='p', gap=ANY, comp='endo') %
         ErbB3(d=1, rtk=None, atp=None, state='p', gap=None, comp='endo') % SHC(grb=None), k6b, kd6b)

    # ErbB2:ErbB2:Shc['p' or 'up']
    Rule('transport_ErbB22_shc', ErbB2(d=1, rtk=None, atp=None, state='p', gap=ANY, comp='pm') %
         ErbB2(d=1, rtk=None, atp=None, state='p', gap=None, comp='pm') % SHC(grb=None) <>
         ErbB2(d=1, rtk=None, atp=None, state='p', gap=ANY, comp='endo') %
         ErbB2(d=1, rtk=None, atp=None, state='p', gap=None, comp='endo') % SHC(grb=None), k7, kd7)

########################################################
def receptor_dimerization():
    
    " Heterodimers form when a ligand bound monomer binds a naked monomer, as shown below "
    " Lig:monomer:ATP + monomer <-> Lig:dimer *** Note that ATP is missing from the right hand side product "

    " ErbB1 homodimers form when two ligand bound, ATP bound monomers dimerize. ErbB2 homodimers cannot be formed by this route \
    since they do not bind ligand, instead by secondary dimerization described in a separate module "
    
    Parameter('k2', 7.44622e-6)
    Parameter('kd2', 0.16)
    Parameter('k2b', 3.73632e-8)
    Parameter('kd2b', 0.016)
    Parameter('k120', 1.48131e-8)
    Parameter('kd120', 0.1)
    Parameter('k120b', 5.92538e-11)
    Parameter('kd120b', 0.1)

    alias_model_components()
    
    # Note that we have to "uncorrect" for BNG's correction for symmetric
    # reactions to match the original model.
    Expression('k2_symmetric', k2 * 2)

    alias_model_components()

    # Rules
    # =====
    for place in comprtmnts:
        # EGF:ErbB1:ATP + ErbB[2-4] <-> ErbB1:Erb[2-4]
        for erb in receptors[1:]:
            Rule(place+'_bind_egf_E1_atp_' + erb.name, ATP(erb=2,gab1=None) % ErbB1(lig=ANY,atp=2,d=None, state='up',comp=place) +
                 erb(lig=None,d=None,atp=None, state='up', comp=place) <>
                 ErbB1(lig=ANY, atp=None, d=1, state='up', comp=place) % erb(lig=None,d=1, atp=None, state='up', comp=place),
                 k2b, kd2b)
    
        # EGf:ErbB1:ATP + EGF:ErbB1:ATP <-> 2(EGF:ErbB1:ATP)
        Rule(place+'_bind_egf_E1_atp_egf_ErbB1', ATP(erb=2,gab1=None) % ErbB1(lig=ANY,atp=2,d=None,state='up', comp=place) +
             ATP(erb=3,gab1=None) % ErbB1(lig=ANY, d=None,atp=3,state='up', comp=place) <>
             ATP(erb=2,gab1=None) % ErbB1(lig=ANY, atp=2, d=1,state='up', comp=place) % ErbB1(lig=ANY, d=1, atp=3,state='up', comp=place) %
             ATP(erb=3,gab1=None), k2_symmetric, kd2)
            
        # HRG:ErbB[3-4] + ErbB2 <-> ErbB2:ErbB[3-4]
        bind_table([[   ErbB2(atp=None, state='up', comp=place)],
                     [ErbB3(lig=ANY,atp=None, comp=place), (k120, kd120)],
                     [ErbB4(lig=ANY,atp=None, comp=place), (k120, kd120)]],'d', 'd')
         
        # HRG:Erb[3-4] + ErbB1:ATP <-> Erb[3-4]:ErbB1
        " Refactor to single rule. If I remeber correctly, the rules could not be put in a table beacuase of need to specify bound ATP "
        Rule(place+'_E1_atp__H_E3', ErbB3(lig=ANY, atp=None, d=None, comp=place) +
             ATP(erb=2,gab1=None) % ErbB1(lig=None,atp=2, d=None, comp=place)
              <> ErbB3(lig=ANY,atp=None, d=1, comp=place) % ErbB1(lig=None, atp=None, d=1, comp=place), k120b, kd120)
             
        Rule(place+'_E1_atp__H_E4', ErbB4(lig=ANY, atp=None, d=None, comp=place) +
             ATP(erb=2,gab1=None) % ErbB1(lig=None,atp=2, d=None,comp=place) <>
             ErbB4(lig=ANY,atp=None, d=1, comp=place) % ErbB1(lig=None, atp=None, d=1, comp=place), k120b, kd120)

########################################################
def secondary_dimerization():
    " v673-v681, When phosphorylated dimers dissociate, dissociation products can pair up in new combinations. "
    " ErbB2 homodimers are formed as a result of one such secondary dimerization process "
    
    Parameter('k102', 5e-07)
    Parameter('kd102', 5.61009)
    Parameter('k103', 8.36983e-9)
    Parameter('kd103', 0.016)
    Parameter('k96', 1.67e-6)
    Parameter('kd96', 0.1)
    
    alias_model_components()
    
    # Note that we have to "uncorrect" for BNG's correction for symmetric
    # reactions to match the original model.
    Expression('k102_symmetric', k102 * 2)
    Expression('k96_symmetric', k96 * 2)

    alias_model_components()

    " secondary dimerization does not take place in endo compartment. DO NOT loop over compartments "
    # 2(EGF:ErbB1)~P <-> EGF:ErbB1~P + EGF:ErbB1~P
    # Note that we have to "uncorrect" for BNG's correction for symmetric
    # reactions to match the original model.
    Rule('Unbind_EGF_ErbB1_P_', ErbB1(lig=ANY, state='p',atp=None,d=1,gap=None,rtk=None, comp='pm')
             % ErbB1(lig=ANY, state='p',atp=None,d=1,gap=None, rtk=None, comp='pm') <>
             ErbB1(lig=ANY, state='p',atp=None,d=None,gap=None, rtk=None, comp='pm') +
             ErbB1(lig=ANY, state='p',atp=None,d=None, gap=None, rtk=None, comp='pm'),
             kd102, k102_symmetric)
    
    # EGF:ErbB1~P + ErbB[2-4]~P <-> ErbB1:ErbB[2-4]~P
    for erb in receptors[1:]:
            Rule('pbind_ErbB1_'+ erb.name, EGF(rec=2, comp='pm') % ErbB1(lig=2, state='p', d=None, atp=None,gap=None, rtk=None, comp='pm') +
                 erb(lig=None, state='p', d=None, atp=None,gap=None, rtk=None, comp='pm') <>
                 ErbB1(lig=None, state='p', d=1,atp=None,gap=None, rtk=None, comp='pm') %
                 erb(lig=None, state='p', d=1, atp=None,gap=None, rtk=None, comp='pm'), k102, kd102)
    
    # ErbB2~P + Erb[1-4]~P <-> (ErbB2:ErbB[1-4]~P
    # "Uncorrect" for ErbB2:2 reaction here too.
    common = {'state': 'p', 'gap': None, 'rtk': None, 'comp': 'pm'}
    bind_table([[                 ErbB2(**common),       ErbB3(lig=None, **common), ErbB4(lig=None, **common)],
                [ErbB2(**common), (k96_symmetric, kd96), (k103, kd103),             (k103, kd103)]],
               'd', 'd')

########################################################
def lateral_signaling():
    
    Parameter('k103_ls', 8.36983e-9)    # k103
    Parameter('kd103_ls', 0.016)    # kd103
    Parameter('k122_ls', 1.8704e-8)     # k122
    Parameter('kd122_ls', 1)          # kd122
    Parameter('kd123_ls', 0.177828)   # kd123
    
    alias_model_components()
    
    # ErbB2~P + ErbB2 <-> ErbB2~P:ErbB2
    Rule('bind_ErbB2_ErbB2_up_', ErbB2(d=None, state='p', gap=None, rtk=None, comp='pm', atp=None) +
             ErbB2(d=None, state='up',gap=None, rtk=None, comp='pm', atp=None) <>
             ErbB2(d=1, state='p',gap=None, rtk=None, comp='pm', atp=None) %
            ErbB2(d=1,state='up',gap=None, rtk=None, comp='pm', atp=None),k103_ls,kd103_ls)

    # ErbB2~P:ErbB2 -> 2(ErbB2)~P
    Rule('ErbB2_tp_ErbB2_up', ErbB2(d=1, state='p', comp='pm', atp=None, gap=None, rtk=None, cpp=None) %
         ErbB2(d=1, state='up', atp=None, gap=None, comp='pm', rtk=None, cpp=None) +ATP(erb=None, gab1=None) <>
     ErbB2(d=1, state='p', atp=None, comp='pm', gap=None, rtk=None, cpp=None) %
         ErbB2(d=1, state='up', atp=2, gap=None, comp='pm', rtk=None, cpp=None)% ATP(erb=2, gab1=None), k122_ls, kd122_ls)
    # 2(ErbB2)~P --> ErbB2~P + ErbB2~P
    Rule('ErbB2_tp_cat_', ErbB2(d=1, state='p', atp=None, comp='pm') % ErbB2(d=1, state='up', atp=2, comp='pm') % ATP(erb=2, gab1=None) >>
              ErbB2(d=1, state='p', atp=None, comp='pm') % ErbB2(d=1, state='p', atp=None, comp='pm') +  ATP(erb=None, gab1=None), kd123_ls)

#########################
def ErbB2_magic_rxns():
    """ Classified as magic reaction since ErbB2 dimerizes with ErbB3 or ErbB4 in the absence of any ligand.
        The unholy dimers proceed to bind (consume) EGF even though it is promised only to ErbB1, while simultneosuly being 
        phoshorylated in the absence of ATP
    """
    Parameter('k103_magic', 8.36983e-9)
    Parameter('kd103_magic', 0.016)
    Parameter('k1c', 800)
    Parameter('kd1c', 1)
    Parameter('k1d', 518)
    Parameter('kd1d', 0.1)

    alias_model_components()

    #ErbB2  + ErbB[3-4] <-> ErbB2:ErbB[3-4]
    for erb in receptors[2:]:
        Rule('bind_ErbB2_' +erb.name, ErbB2(d=None, atp=None, state='up', comp='pm', gap=None, cpp=None) +
             erb(lig=None, d=None, atp=None, state='up', comp='pm', gap=None, cpp=None) <>
             ErbB2(d=1, atp=None, state='up', comp='pm', gap=None, cpp=None) %
             erb(lig=None, d=1, atp=None, state='up', comp='pm', gap=None, cpp=None), k103_magic, kd103_magic)
    # EGF + ErbB2:ErbB3 <-> ErbB2~P:ErbB3~P
    Rule('phosphorylate_ErbB2_ErbB3', EGF(rec=None, comp='pm') + ErbB2(d=1, atp=None, state='up', comp='pm', gap=None, cpp=None)
         % ErbB3(lig=None, d=1, atp=None, state='up', comp='pm', gap=None, cpp=None) <>
         ErbB2(d=1, atp=None, state='p', comp='pm', gap=None, cpp=None) %
         ErbB3(lig=None, d=1, atp=None, state='p', comp='pm', gap=None, cpp=None), k1c, kd1c)
    # EGF + ErbB2:ErbB4 <-> ErbB2~P:ErbB4~P
    Rule('phosphorylate_ErbB2_ErbB4', EGF(rec=None, comp='pm') + ErbB2(d=1, atp=None, state='up', comp='pm', gap=None, cpp=None) %
         ErbB4(lig=None, d=1, atp=None, state='up', comp='pm', gap=None, cpp=None) <>
         ErbB2(d=1, atp=None, state='p', comp='pm', gap=None, cpp=None) %
         ErbB4(lig=None, d=1, atp=None, state='p', comp='pm', gap=None, cpp=None), k1d, kd1d)

########################################################
def trans_phosphorylation():
    "v21-v29, v37-v41, v806-814, v822-827 "
    " Lig:dimer + ATP <-> Lig:dimer:ATP -> dimer~P + ATP *** Note that ligand dissapears in the ES -> E + P step "
    " check v29 since product end sup in the endosome"
    
    Parameter('k122', 1.8704e-8)    # k122
    Parameter('kd122', 1)         # kd122
    Parameter('kd123', 0.177828)  # kd123
    Parameter('k123', 0)          # k123
    
    alias_model_components()
    
    " Unlike other classes of dimers, note that Phoshporylation of ErbB[3-4]:ErbB1 dimers take place only in the plasma membrane "
    for s in receptors[2:]:
        # HRG:ErbB[3-4]:ErbB1 + ATP <-> HRG:Erb[3-4]:ErbB1:ATP
        Rule(s.name+'_tp_bind_', s(lig=ANY,d=1, state='up', atp=None, comp='pm') %
             ErbB1(lig=None,d=1, state='up', atp=None) + ATP(erb=None,gab1=None) <>
             s(lig=ANY,d=1, state='up', atp=2, comp='pm') % ErbB1(lig=None,d=1, state='up', atp=None)
             % ATP(erb=2, gab1=None), k122, kd122)
        # HRG:Erb[3-4]:ErbB1:ATP --> HRG:Erb[3-4]~P:ErbB1~P + ATP
        Rule(s.name+'_tp_cat_', HRG(rec=3, comp='pm') % s(lig=3,d=1, state='up', atp=2, comp='pm', gap=None) %
                  ErbB1(lig=None,d=1, state='up', atp=None, gap=None) % ATP(erb=2,gab1=None) >>
                  s(lig=None,d=1, state='p', atp=None, comp='pm', gap=None) % ErbB1(lig=None,d=1, state='p', atp=None, gap=None) +
                  ATP(erb=None,gab1=None), kd123)
    
    " Phoshporylation of dimers below takes place in both pm and endo compartments "
    for place in comprtmnts:
        for s in receptors[2:]:
            # HRG:Erb[3-4]:ErbB2 + ATP <-> HRG:Erb[3-4]:ErbB2:ATP
            Rule(place + '_' + s.name+'_tp_bind_', s(lig=ANY,d=1, state='up', atp=None, comp=place) %
                 ErbB2(lig=None,d=1, state='up', atp=None) + ATP(erb=None,gab1=None) <>
                 s(lig=ANY,d=1, state='up', atp=2, comp=place) % ErbB2(lig=None,d=1, state='up', atp=None) % ATP(erb=2, gab1=None), k122, kd122)
            # HRG:Erb[3-4]:ErbB2:ATP -> Erb[3-4]~P:ErbB2~P + ATP
            Rule(place + '_' + s.name+'_tp_cat_', HRG(rec=3, comp=place) % s(lig=3,d=1, state='up', atp=2, comp=place, gap=None) %
                      ErbB2(lig=None,d=1, state='up', atp=None, gap=None) % ATP(erb=2,gab1=None) >>
                      s(lig=None,d=1, state='p', atp=None, comp=place, gap=None) % ErbB2(lig=None,d=1, state='p', atp=None, gap=None) +
                      ATP(erb=None,gab1=None), kd123)
    
        for r in receptors[1:]:
            # EGF:ErbB1::Erb[2-4] + ATP -> ATP:EGF:ErbB1:Erb[2-4]
            Rule(place+'_ErbB1_tp_bind_'+ r.name, ErbB1(lig=ANY,d=1, state='up', atp=None, comp=place) %
                 r(lig=None,d=1, state='up', atp=None) + ATP(erb=None,gab1=None) <>
                 ErbB1(lig=ANY,d=1, state='up', atp=2, comp=place) % r(lig=None,d=1, state='up', atp=None)% ATP(erb=2,gab1=None), k122, kd122)
            # ATP:EGF:ErbB1:Erb[2-4] --> ErbB1~P:Erb[2-4]~P + ATP
            Rule(place +'_ErbB1_tp_cat_'+ r.name, EGF(rec=3)%ErbB1(lig=3,d=1, state='up', atp=2, comp=place) %
                      r(lig=None,d=1, state='up', atp=None)% ATP(erb=2,gab1=None) >>
                      ErbB1(lig=None,d=1, state='p', atp=None, comp=place) % r(lig=None,d=1, state='p', atp=None) +
                 ATP(erb=None,gab1=None), kd123)
        
        " Note that a thrid ATP is added to the ErbB1 homodimers, resulting in a species called 'full active' ErbB1 dimers "
        " In the PySB model, this is addressed by having a new state called full active. Hence the change of state in the phoshorylation\
        state is: 'up' <-> 'full_active' -> 'p'  "
        # 2(EGF:ErbB1:ATP) + ATP <-> 2(EGF:ErbB1:ATP)_full_active 2(EGF:ErbB1)~P
        Rule('ErbB1_tp_bind_ErbB1_'+ place, MatchOnce(ErbB1(lig=ANY,d=1, state='up', atp=3, comp=place)% ATP(erb=3, gab1=None) %
             ErbB1(lig=ANY,d=1, state='up', atp=4)% ATP(erb=4, gab1=None)) + ATP(erb=None, gab1=None) <>
             MatchOnce(ATP(erb=2,gab1=None) % ATP(erb=3, gab1=None) % ErbB1(lig=ANY,d=1, state='full_act', atp=3, comp=place)
             % ErbB1(lig=ANY,d=1, state='full_act', atp=2)) , k122, kd122)
        # 2(EGF:ErbB1:ATP)_full_active -> 2(EGF:ErbB1)~P + ATP
        Rule('ErbB1_tp_cat_ErbB1_'+ place, MatchOnce(ErbB1(lig=ANY,d=1, state='full_act', atp=2, comp=place) %
                  ErbB1(lig=ANY,d=1, state='full_act', atp=3)% ATP(erb=2,gab1=None)% ATP(erb=3,gab1=None)) >>
                  ErbB1(lig=ANY,d=1, state='p', atp=None, comp=place) % ErbB1(lig=ANY,d=1, state='p', atp=None)
                  + ATP(erb=None,gab1=None) , kd123)

########################################################
def GAP_binding():
    "v194-v207"
    " In the Chen 2009 model, each dimer binds to a single unit of downstream substrates like GAP, Grb2, SHC. To reproduce this \
    reaction, GAP (and Grb2, SHc) bind only the ErbB1 receptor in all its dimers or ErbB2 in its dimers. In ErbB1:ErbB2 dimers, GAP binds ErbB1 "
    " ErbB1:Erb[3-4] dimers in both compartments bind GAP with rate constants k8b and k8bd "
    " ErbB1:ErbB2 dimer in pm comparmtment binds with rate constants k8b and k8bd  "
    
    # Initial amount
    # ==============
    Parameter('GAP_0', 534751)
    # Rate constants
    # ==============
    Parameter('k8', 5.91474e-7)
    Parameter('kd8', 0.2)
    Parameter('k8b', 9.34641e-6)
    Parameter('kd8b', 0.02)
    
    alias_model_components()
    
    # Initial conditions
    # ==================
    Initial(GAP(rec=None), GAP_0)
    
    for place in comprtmnts:
        # ErbB1:ErbB[1-4] + GAP <-> GAP:ErbB1:ErbB[1-4]
        for erb in receptors:
            rates = (k8b, kd8b)
            if erb is ErbB1 or (erb is ErbB2 and place == 'pm'):
                rates = (k8, kd8)
            Rule(place+'_gap_bind_ErbB1_'+erb.name, MatchOnce(ErbB1(state='p', d=1, gap=None,gs=None, rtk=None, comp=place) %
                 erb(state='p',d=1, gap=None,gs=None, rtk=None)) + GAP(rec=None) <>
                 MatchOnce(ErbB1(state='p', d=1, gap=2,gs=None, rtk=None, comp=place) % erb(state='p',d=1, gap=None,gs=None, rtk=None) %
                 GAP(rec=2)), *rates)
        
        # ErbB2:Erb[2-4] + GAP <-> GAP:ErbB2:Erb[2-4]
        for erb in receptors[1:]:
            Rule(place+'_gap_bind_ErbB2_'+erb.name, MatchOnce(ErbB2(state='p', d=1, gap=None,gs=None, rtk=None, comp=place) %
                 erb(state='p',d=1,gap=None,gs=None, rtk=None)) + GAP(rec=None) <>
                 ErbB2(state='p', d=1, gap=2,gs=None, rtk=None, comp=place) % erb(state='p',d=1,gap=None,gs=None, rtk=None) %
                 GAP(rec=2), k8, kd8)

########################################################
def Grb2_binding():
    " v212- v240 Grb2 competes with SHC to bind GAP bound receptors, and dissocaites under all conditions except when SOS is bound to RAS "
    
    # Initial value
    # =============
    Parameter('Grb2_0', 1264.91)    # c22 //// Grb2:SOS has non-zero initial condition (c30)
    Parameter('SOS_0', 0)           # c24 // zero ic
    Parameter('Grb2SOS_0', 8.8914e7)  # c30
    # Rate constants
    # ==============
    Parameter('k16', 1.67e-05)     # k16
    Parameter('kd24', 0.55)         # kd24 = 0.55 or kd63 = 0.275; kd63 appears in 6 of the 29 backward reactions, no apparent pattern observed :(
    Parameter('kd63', 0.275)
    Parameter('k17',  1.67e-05)         # k17
    Parameter('kd17', 0.06)              # kd17
    Parameter('k34', 7.5e-6)
    Parameter('kd34', 0.03)
    Parameter('k35', 7.5e-6)
    Parameter('kd35', 0.0015)
    Parameter('k25', 1.67e-5)
    Parameter('kd25', 0.0214)
    Parameter('k101', 8.33e-07)         # 1e-8
    Parameter('kd101', 0.03)
    Parameter('k40', 5e-5)
    Parameter('kd40', 0.064)
    
    alias_model_components()
    
    # Initial conditions
    # ===================
    Initial(Grb2(shc=None, erb=None, sos=None, gab1=None), Grb2_0)
    Initial(SOS(grb=None, ras=None, erk=None, state='up'), SOS_0)
    Initial(Grb2(shc=None, erb=None, sos=1, gab1=None) % SOS(grb=1, ras=None, erk=None,state='up'), Grb2SOS_0)
    
    # Rules
    # =====
    " SOS binds Grb2 that is not bound to receptors "
    # Grb2:Shc~P binds unphoshphoryalted SOS
    Rule('bind_free_Grb2SchP_SOS', SHC(state='p', erb=None, grb=2) % Grb2(sos=None, erb=None, gab1=None, shc=2) + SOS(grb=None,ras=None,erk=None, state='up') <>
         SHC(state='p', erb=None, grb=2) % Grb2(sos=1, erb=None,gab1=None, shc=2) % SOS(grb=1, ras=None, erk=None, state='up'), k40, kd40)
        
    # Grb2 binds unphoshphoryalted SOS
    Rule('bind_free_Grb2_SOS', Grb2(sos=None, erb=None, gab1=None, shc=None) + SOS(grb=None,ras=None,erk=None, state='up') <>
          Grb2(sos=1, erb=None,gab1=None, shc=None) % SOS(grb=1, ras=None, erk=None, state='up'), k35, kd35)
         
         
    for place in comprtmnts:
        " Grb2  binds to receptors "
        for erb in receptors:
             # Grb2 binds ErbB1:Erb[1-4] receptor dimers
             rate_reverse = kd24
             if erb is ErbB1:
                 rate_reverse = kd63
             Rule(place+'_Grb2_bind_ErbB1_'+erb.name, ErbB1(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) +
                  Grb2(erb=None,shc=None,gab1=None, sos=None) <>
                  Grb2(erb=1,shc=None,gab1=None, sos=None) % ErbB1(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None), k16, rate_reverse)
             # Grb2:SOS binds ErbB1:Erb[1-4] receptor dimers
             Rule(place+'_Grb2sos_bind_ErbB1_'+erb.name, ErbB1(gap=ANY,gs=None,d=2,  cpp=None, comp=place) %
                  erb(d=2,gs=None) + SOS(grb=3, ras=None, erk=None,state='up') % Grb2(erb=None,shc=None,gab1=None, sos=3) <>
                  SOS(grb=3, ras=None,erk=None, state='up') % Grb2(erb=1,shc=None,gab1=None, sos=3) %
                  ErbB1(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),k34, kd34)
                     
        for erb in receptors[1:]:
             # Grb2 binds ErbB2:ErbB[2-4] receptor dimers
             rate_reverse = kd24
             if erb is ErbB2 or (erb is ErbB3 and place == 'pm') or (erb is ErbB4 and place == 'endo'):
                 rate_reverse = kd63
             Rule(place+'_Grb2_bind_ErbB2_'+ erb.name, ErbB2(gap=ANY,gs=None,d=2,  cpp=None, comp=place) % erb(d=2,gs=None) +
                  Grb2(erb=None,shc=None,gab1=None, sos =None) <>
                  Grb2(erb=1,shc=None,gab1=None, sos=None) % ErbB2(gap=ANY,gs=1,d=2,  cpp=None, comp=place) % erb(d=2,gs=None), k16, rate_reverse)
             # Grb2:SOS binds ErbB2:ErbB[2-4] receptor dimers
             Rule(place+'_Grb2sos_bind_ErbB2_'+ erb.name, ErbB2(gap=ANY,gs=None,d=2,  cpp=None, comp=place) % erb(d=2,gs=None) +
                  SOS(grb=3, ras=None,erk=None, state='up') %Grb2(erb=None,shc=None,gab1=None, sos=3) <>
                  SOS(grb=3, ras=None, erk=None, state='up') % Grb2(erb=1,shc=None,gab1=None, sos=3) %
                  ErbB2(gap=ANY,gs=1,d=2,  cpp=None, comp=place) % erb(d=2,gs=None),k34, kd34)
                 
        " SOS binds Grb2 that is bound to receptors "
        # SOS binds GRb2 in dimers
        for erb in receptors[:2]:
             # SOS binds GRb2 in dimers
             Rule(place+'_bind_Grb2_SOS_noSHC'+ erb.name, erb(cpp=None, comp=place, gap=ANY) %
                  Grb2(shc=None, erb=ANY, sos=None,gab1=None) + SOS(grb=None,erk=None, ras=None, state='up') <>
                  erb(cpp=None, gap=ANY, comp=place) % Grb2(shc=None, erb=ANY, sos=1,gab1=None) %
                  SOS(grb=1, ras=None, erk=None, state='up'), k17, kd17)
                  
             #for erb in receptors[:2]:
             # SHC~P:Grb2 binds SOS
             Rule(place+'_bind_Grb2_SOS_'+ erb.name, erb(cpp=None, comp=place, gap=ANY) % Grb2(shc=ANY, erb=None, sos=None,gab1=None) +
                  SOS(grb=None,erk=None, ras=None, state='up') <> erb(cpp=None, gap=ANY, comp=place) %
                  Grb2(shc=ANY, erb=None, sos=1,gab1=None) % SOS(grb=1, ras=None, erk=None, state='up'), k25, kd25)
                 
        " SOS~P binds Grb2 only in ErbB1 homodimers "
        Rule(place+'_bind_Grb2_SOS', ErbB1(d=1, cpp=None, comp=place, gap=ANY) % ErbB1(d=1, cpp=None, comp=place, gap=None) %
              Grb2(sos=None, gab1=None)+ SOS(grb=None,erk=None, ras=None, state='p') <>
              ErbB1(d=1, cpp=None, comp=place, gap=ANY)%ErbB1(d=1, cpp=None, comp=place, gap=None) %
              Grb2(sos=2, gab1=None) % SOS(grb=2,erk=None, ras=None, state='p'), k101, kd101)


########################################################
def SHC_binding():
    " v367-v380"
    " SHC competes with Grb2 to bind GAP primed receptors "
    " SHC can dissociate from receptors under any condition other than when cPP is bound to the receptor, or RAS is bound to SOS, \
    giving rise to dissocation products like SHC, SHC~P. SHC~P:Grb2, SHC~P:Grb2:SOS"
    " The binding affinities differ based on whether Shc is unphoshorylated (k22, kd22), phoshporyated (k37, kd37), \
    bound to Grb2(k37/kd37) or bound to Grb2:SOS(k32/kd32)"
    
    # Initial amount
    Parameter('SHC_0', 11e+5)   # c31
    # Rate constants
    Parameter('k22', 1.39338e-7)   # k22
    Parameter('kd22', 0.1)      # kd22
    Parameter('kd22b', 0.1)     # kd22
    Parameter('k37', 1.5e-6)    # For Shc~p:Grb2 and Shc~P
    Parameter('kd37', 0.3)
    Parameter('k32', 4e-7)      # For Shc~p:Grb2:SOS
    Parameter('kd32', 0.1)
    Parameter('k23', 6)
    Parameter('kd23', 0.06)
    Parameter('k36', 0.005)
    Parameter('kd36', 0)
    
    alias_model_components()
    
    # Initial condition
    # =================
    Initial(SHC(erb=None, grb=None, state='up'), SHC_0)
    
    # Rules
    # =====
    # dephoshorylation of free Shc~P
    Rule('shc_phos', SHC(erb=None, grb=None, state='p') <> SHC(erb=None, grb=None, state='up'), k36, kd36)
    
    for place in comprtmnts:
        # SHC binds ErbB1:Erb[1-4] receptor dimers
        for erb in receptors:
            rate_reverse = kd22b
            if erb is ErbB1:
                rate_reverse = kd22
            Rule(place+'_SHC_bind_ErbB1_'+erb.name, ErbB1(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) +
                 SHC(erb=None, grb=None, state='up') <> SHC(erb=1, grb=None, state='up') %
                 ErbB1(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None), k22, rate_reverse)
        
        # SHC binds ErbB2:Erb([2-4] receptor dimers
        for erb in receptors[1:]:
            rate_reverse = kd22
            if erb is ErbB2 or erb is ErbB3:
                rate_reverse = kd22b
            Rule(place+'_SHC_bind_ErbB2_'+erb.name, ErbB2(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) +
                 SHC(erb=None, grb=None, state='up') <> SHC(erb=1, grb=None, state='up') %
                 ErbB2(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None), k22, rate_reverse)
        
        # SHC~P binds ErbB1:Erb[1-4] receptor dimers
        for erb in receptors:
            Rule(place+'_SHCP_bind_ErbB1_'+erb.name, ErbB1(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) +
                 SHC(erb=None, grb=None, state='p') <> SHC(erb=1, grb=None, state='p') %
                 ErbB1(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),k37, kd37)
        
        # SHC~P binds ErbB2:Erb[2-4] receptor dimers
        for erb in receptors[1:]:
            Rule(place+'_SHCP_bind_ErbB2_'+erb.name, ErbB2(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) +
                 SHC(erb=None, grb=None, state='p') <> SHC(erb=1, grb=None, state='p') %
                 ErbB2(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),k37, kd37)
        
        # SHC:Grb binds ErbB1:Erb([1-4] receptor dimers
        for erb in receptors:
            Rule(place+'_SHC_grb_bind_ErbB1_'+erb.name, ErbB1(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) +
                 Grb2(shc=3, sos=None) % SHC(erb=None, grb=3) <> Grb2(shc=3, sos=None) % SHC(erb=1, grb=3) %
                 ErbB1(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),k37, kd37)
        
        # SHC:Grb binds ErbB2:Erb[2-4] receptor dimers
        for erb in receptors[1:]:
            Rule(place+'_SHC_grb_bind_ErbB2_'+erb.name, ErbB2(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) +
                 Grb2(shc=3, sos=None) %SHC(erb=None, grb=3) <> Grb2(shc=3, sos=None) % SHC(erb=1, grb=3) %
                 ErbB2(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),k37, kd37)
        
        # SHC:Grb:SOS binds ErbB1:Erb[1-4] receptor dimers
        for erb in receptors:
            Rule(place+'_SHC_grbsos_bind_ErbB1_'+erb.name, ErbB1(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) +
                 SOS(ras=None,erk=None, state='up') % SHC(erb=None, grb=ANY) <> SOS(ras=None, erk=None, state='up') %
                 SHC(erb=1, grb=ANY) % ErbB1(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),k32, kd32)
        
        # SHC:Grb:SOS binds ErbB2:Erb[2-4] receptor dimers
        for erb in receptors[1:]:
            Rule(place+'_SHC_grbsps_bind_ErbB2_'+erb.name, ErbB2(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) +
                 SOS(ras=None, erk=None, state='up') % SHC(erb=None, grb=ANY) <> SOS(ras=None, erk=None, state='up') %
                 SHC(erb=1, grb=ANY) % ErbB2(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),k32, kd32)
        
    # Phosphorylation of recepotr bound SHC, SHC <-> SHC~P
    Rule('shc_bound_phos', SHC(erb=ANY, grb=None, state='up') <> SHC(erb=ANY, grb=None, state='p'), k23, kd23)

########################################################
def secondary_Grb2_binding():
    " Grb2 binds SHC~P, and dissocaties under all conditions except when RAS is bound to SOS "
    
    Parameter('k16_scndry',  1.67e-05)     # k16
    Parameter('kd24_scndry',  0.55)     # kd24 = 0.55 or kd63 = 0.275; kd63 appears in 6 of the 29 backward reactions, no apparent pattern observed :(
    Parameter('k41',  5e-05)         # k17
    Parameter('kd41', 0.0429)
    Parameter('k33',  3.5e-05)         # k17
    Parameter('kd33', 0.2)
    
    alias_model_components()
    
    # free Grb2 binds free SHC~P
    Rule('shc_grb_binding', SHC(erb=None, grb=None, state='p') + Grb2(erb=None, sos=None, shc=None,gab1=None) <>
         SHC(erb=None, grb=1, state='p') % Grb2(erb=None, sos=None, shc=1,gab1=None), k16_scndry, kd24_scndry)
    
    # Shc~P + Grb2:SOS <-> Shc~P:Grn2:SOS
    Rule('SHCp_binds_grb2sos', SHC(grb=None, erb=None, state='p') + Grb2(erb=None, shc=None, sos=ANY) % SOS(ras=None, state='up') <>
         SHC(grb=1, erb=None, state='p') % Grb2(erb=None, shc=1, sos=ANY) % SOS(ras=None, state='up'), k33, kd33)
        
    for place in comprtmnts:
             
         for erb in receptors[:2]:
             # Grb2 unbinds from SHC
             Rule(place+'_shc_grb_unbinding_'+erb.name, erb(cpp=None,gs=2,comp=place) % SHC(grb=1,erb=2, state='p') %
                  Grb2(erb=None,sos=None, shc=1, gab1=None) <> erb(cpp=None,gs=2,comp=place) % SHC(grb=None,erb=2, state='p') +
                  Grb2(erb=None,sos=None, shc=None,gab1=None), kd24_scndry, k16_scndry)
             # Grb2:SOS unbinds from SHC    
             Rule(place+'_shc_grbsos_unbinding_'+erb.name, erb(cpp=None,gs=3,comp=place) % SHC(grb=1, erb=3, state='p') %
                  Grb2(erb=None,sos=2, shc=1, gab1=None) % SOS(grb=2, ras=None, erk=None, state='up') <>
                  erb(cpp=None,gs=3,comp=place) % SHC(grb=None, erb=3, state='p') + Grb2(erb=None,sos=2, shc=None,gab1=None) %
                  SOS(grb=2, erk=None, ras=None, state='up'), kd41, k41)

########################################################
def RAS_binds_sos():
    
    "4 typical reactions:(1) dimer + RAS_GDP <-> dimer:RAS_GDP (2) dimer + RAS_GTP <-> dimer:RAS_GDP \
    (3) dimer + RAS_GDP <-> dimer:RAS_GTP (4) dimer + RAS_active_GTP <-> dimer:RAS_GTP "
    " The following reactions are missing from the above pattern "
    "endosomal dimers + RASGTP/RAS_actvie GTP <> "
    
    Parameter('RAS_gdp_0', 58095.2) # c26
    
    Parameter('k18', 2.5e-5)    # k18
    Parameter('kd18', 1.3)     # kd18
    
    Parameter('k19', 1.667e-7)    # k19
    Parameter('kd19', 0.5)      # kd19
    
    Parameter('k21', 3.67e-7)    # k21
    Parameter('kd21', 0.23)     # kd21
    
    Parameter('k20', 1.1068e-5)  # k20
    Parameter('kd20', 0.4)     # kd20
    
    alias_model_components()
    
    # Initial conditions
    # ==================
    Initial(RAS(sos=None, pi3k=None, raf=None,state='gdp'), RAS_gdp_0)
    
    # RAS binds SOS
    for place in comprtmnts:
        for erb in receptors[:2]:
            Rule(place+'_Ras_gdp_gdp_'+erb.name, erb(cpp=None,gap=ANY,comp=place) % GAP() %
                 SOS(grb=ANY,erk=None, ras=None, state='up') + RAS(sos=None, pi3k=None, raf=None, state='gdp') <>
                 erb(cpp=None, gap=ANY, comp=place) % GAP() % SOS(grb=ANY,erk=None, ras=1, state='up') %
                 RAS(sos=1, pi3k=None, raf=None, state='gdp'), k18, kd18)
                
                
            Rule(place+'_Ras_gdp_gtp_'+erb.name, erb(cpp=None, gap=ANY, comp=place) % GAP() %
                      SOS(grb=ANY,erk=None, ras=None, state='up') + RAS(sos=None, pi3k=None, raf=None, state='gdp') <>
                      erb(cpp=None, gap=ANY, comp=place) % GAP() % SOS(grb=ANY,erk=None, ras=1, state='up') %
                      RAS(sos=1, pi3k=None, raf=None, state='gtp'), k21, kd21)

    "  dimers + RAS_gtp/active_gtp reactions occur only in plasma membrane "
    for erb in receptors[:2]:
        Rule('Ras_agtp_gtp'+erb.name, erb(cpp=None, gap=ANY,comp='pm') % GAP() % SOS(grb=ANY,erk=None, ras=None, state='up') +
             RAS(sos=None, pi3k=None, raf=None,state='active_gtp') <> erb(cpp=None, gap=ANY, comp='pm') % GAP() %
             SOS(grb=ANY,erk=None, ras=1, state='up') % RAS(sos=1, pi3k=None, raf=None,state='gtp'), k20, kd20)
            
        Rule('Ras_gtp_gdp_'+erb.name, erb(cpp=None, gap=ANY, comp='pm') % GAP() % SOS(grb=ANY,erk=None, ras=None, state='up') +
                 RAS(sos=None, pi3k=None, raf=None, state='gtp') <> erb(cpp=None, gap=ANY, comp='pm') % GAP() %
                 SOS(grb=ANY, erk=None, ras=1, state='up') % RAS(sos=1, pi3k=None,raf=None,state='gdp'), k19, kd19)

########################################################
def RTK_phos():
    " v650-v663 "
    " Erb dimers that are phoshproylated, but not bound to GAP can bind RTK (ES), followed by release of Enzyme RTK and \
    dephoshphorylated dimer product "
    "Reaction occurs ony for receptor dimers in the endo compartment "
    " In the case of ErbB1 type dimers, EGF returns in the unphosphorylated product"
    " Also in ErbB1 homodimers, ATP returns to the fold as well :) "
    
    # Initial amount
    # ==============
    Parameter('RTK_0', 7e+4)        # c280
    # Rate constants
    # ==============
    Parameter('k94', 5e-05)     # k94  a subset have rc k94b and kd94b ..the pattern is ErbB1/3 and ErbB1/4 dimers
    Parameter('k94b', 5e-05)    # Unnecesary parameter
    Parameter('kd94', 0.01)      # kd94
    Parameter('kd95', 33)           # kd95
    
    alias_model_components()
    
    # Initial conditions
    # ==================
    Initial(RTK(erb=None), RTK_0)
    
    # Rules
    # ======
    # check if catalyze can be used for these reactions
    
    # ERbB1~P:ErbB[1-2]~P + RTK_phosphatase <-> ERbB1~P:ErbB[1-2]~P : RTK_phosphatase
    for erb in receptors[:2]:
        Rule('bind_RTK_ErbB1_'+erb.name, RTK(erb=None) + MatchOnce(ErbB1(d=1, state='p', rtk=None, gap=None, comp='endo')
             % erb(d=1, state='p', gap=None,rtk=None, comp='endo')) <>
             RTK(erb=2)% ErbB1(d=1, state='p', rtk=2, gap=None, comp='endo') %
             erb(d=1, gap=None, state='p',rtk=None, comp='endo'), k94, kd94)
    
    # ERbB1~P:ErbB[3-4]~P + RTK_phosphatase <-> ERbB1~P:ErbB[3-4]~P : RTK_phosphatase
    for erb in receptors[2:]:
        Rule('bind_RTK_ErbB1_'+erb.name, RTK(erb=None) + ErbB1(d=1, state='p', rtk=None, gap=None, comp='endo')
             % erb(d=1, state='p', gap=None,rtk=None, comp='endo') <>
             RTK(erb=2)% ErbB1(d=1, state='p', rtk=2, gap=None, comp='endo') %
             erb(d=1, gap=None, state='p',rtk=None, comp='endo'), k94b, kd94)
    
    " Check for redundancy in rules for  ErbB1 homodimers "
    # ERbB1~P:ErbB[2-4]~P : RTK_phosphatase >> EGF:ERbB1:ErbB[2-4] + RTK+ptase
    for erb in receptors[1:]:
        Rule('kcat_RTK_ErbB1_'+erb.name, RTK(erb=2) % ErbB1(lig=None, d=1, state='p', rtk=2, gap=None, comp='endo') %
             erb(d=1, gap=None, state='p',rtk=None, comp='endo') >> RTK(erb=None) + EGF(rec=3, comp='endo') %
             ErbB1(lig=3,d=1, state='up', rtk=None, gap=None, comp='endo') % erb(d=1, state='up', gap=None,rtk=None, comp='endo'), kd95)

    # Special rule for the return of ATP in EGF
    # ERbB1~P:ErbB1~P : RTK_phosphatase >> 2(EGF:ERbB1:ATP) + RTK+ptase
    Rule('kcat_RTK_ErbB1_ErbB1', RTK(erb=2) % ErbB1(lig=ANY, d=1, state='p', rtk=2, gap=None, atp=None, comp='endo') %
         ErbB1(lig=ANY, d=1, gap=None, state='p',rtk=None, atp=None, comp='endo')>>
         RTK(erb=None) + ATP(erb=2, gab1=None) % ATP(erb=3, gab1=None) %
         ErbB1(lig=ANY,d=1, state='up', rtk=None, atp=2, gap=None, comp='endo') %
         ErbB1(lig=ANY, d=1, state='up', atp=3, gap=None,rtk=None, comp='endo'), kd95)

    for erb in receptors[1:]:
            #ERbB2~P:ErbB[2-4]~P + RTK_phosphatase <-> ERbB2~P:ErbB[2-4]~P : RTK_phosphatase
        Rule('bind_RTK_ErbB2_'+erb.name, RTK(erb=None) + MatchOnce(ErbB2(d=1, state='p', rtk=None, gap=None, comp='endo') %
         erb(d=1, state='p',gap=None,rtk=None, comp='endo')) <>
         RTK(erb=2)% ErbB2(d=1, state='p', rtk=2, gap=None, comp='endo') %
         erb(d=1, gap=None, state='p',rtk=None, comp='endo'), k94, kd94)
        # ERbB2~P:ErbB[2-4]~P : RTK_phosphatase >> ERbB2:ErbB[2-4] + RTK+ptase
        Rule('kcat_RTK_ErbB2_'+erb.name, RTK(erb=2)% ErbB2(d=1, state='p', rtk=2, gap=None, comp='endo') %
             erb(d=1, gap=None, state='p',rtk=None, comp='endo') >>
             RTK(erb=None) + ErbB2(d=1, state='up', rtk=None, gap=None, comp='endo') %
             erb(d=1, state='up', gap=None,rtk=None, comp='endo'), kd95)

########################################################
def bind_cPP():
    " v43-v161 "
    " cPP binds any Grb2 bound dimer. Once bound by cPP, ErbB dimers cannot participate in any reactions \
    (for eg: release Grb, bind SOS, RAS) unless cPP is released "
    " binding of cPP(any compartment) and receptor dimer(any compartment) gives a bound product in the plasma compartment "
    
    # Initial value
    Parameter('cPP_0', 4498.73)     # c12
    # Rate constants
    # ==============
    Parameter('k4', 6.73e-6)                    # forward rate constant for ErbB1 homodimers (plasma mem) binding cPP
    Parameter('k4b', 0)                         # forward rate constant for Erb heterodimers (plasma mem) binding cPP
    Parameter('k5', 0)                          # forward rate constant for ErbB1 dimers and ErbB2 homodimers (endo) binding cPP
    Parameter('k5b', 0)                         # forward rate constant for ErbB2 heterodimers binding cPP
    Parameter('kd4', 1.66e-4)                   # generic reverse rate constant for cPP release into plasma compartment
    Parameter('kd5', 0.80833)                   # reverse rate constant for ErbB1 homodimers cPP release into endo compartment
    Parameter('kd5b', 0.0080833)                # reverse rate constant for cPP release into endo compartment
    
    alias_model_components()
    
    # Initial condition
    Initial(cPP(erb=None, comp='pm'), cPP_0)
    
    # Rules
    # =====
    place = 'pm'
    # EGF:ErbB1 dimers require separate rules since bound ligand is also transported to pm compartment
    # 2(EGF:ErbB1):Grb2(){compartment = ANY) + cPP <-> 2(EGF:ErbB1):Grb2(): CPP { compartment = plasma membrane}
    Rule(place+'_bind_cPP_ErbB1_ErbB1_noSoS', cPP(erb=None,comp=place) + Grb2(sos=None, gab1=None) % EGF(rec=3, comp=place) %
         ErbB1(lig=3, d=1, cpp=None, gap=ANY, comp=place) % EGF(rec=4, comp=place) % ErbB1(lig=4, d=1,gap=None, comp=place) <>
         cPP(erb=2, comp='pm') % Grb2(sos=None, gab1=None) % EGF(rec=3, comp='pm') % ErbB1(lig=3, d=1,cpp=2, gap=ANY, comp='pm') %
         EGF(rec=4, comp='pm') % ErbB1(lig=4, d=1,gap=None, comp='pm'), k4, kd4)
        
    Rule(place+'_bind_cPP_ErbB1_ErbB1', cPP(erb=None,comp=place) + SOS(erk=None, state='up') % Grb2(sos= ANY,gab1=None) % EGF(rec=3, comp=place) %
          ErbB1(lig=3, d=1, cpp=None, gap=ANY, comp=place) % EGF(rec=4, comp=place) % ErbB1(lig=4, d=1,gap=None, comp=place) <>
          cPP(erb=2, comp='pm') % SOS(erk=None, state='up') % Grb2(sos= ANY,gab1=None) % EGF(rec=3, comp='pm') % ErbB1(lig=3, d=1,cpp=2, gap=ANY, comp='pm') %
          EGF(rec=4, comp='pm') % ErbB1(lig=4, d=1,gap=None, comp='pm'), k4, kd4)
         
    # ErbB1:Erb(2/3/4) : Grb2(){compartment = ANY) + cPP <-> ErbB1:ErbB2/3/4:Grb2(): CPP { compartment = plasma membrane}

    def make_dimer_bound(dimer_free):
        cp = dimer_free()
        cp.monomer_patterns[1].site_conditions.update(comp='pm', cpp=9)
        cp.monomer_patterns[2].site_conditions.update(comp='pm')
        return cp
    cpp_free = cPP(erb=None, comp=place)
    cpp_bound = cPP(erb=9, comp='pm')

    for erb in receptors[1:]:
        # Reactions involving dimers WITHOUT Shc and WITH Ras:GTP (v70-v72) use
        # a different rate from the rest, which might have been an error. The
        # multiple rules here are needed to reproduce this, as rules don't let
        # us make the kind of negative assertion we'd need to describe this case
        # concisely.

        dimer_free = Grb2(gab1=None, shc=ANY) % ErbB1(d=1, cpp=None, gap=ANY, comp=place) % erb(d=1, gap=None, comp=place)
        dimer_bound = make_dimer_bound(dimer_free)
        Rule(place+'_bind_cPP_ErbB1_'+erb.name+'_SHC', cpp_free + dimer_free <> cpp_bound % dimer_bound, k4b, kd4)

        dimer_free = dimer_free(shc=None, sos=None)
        dimer_bound = make_dimer_bound(dimer_free)
        Rule(place+'_bind_cPP_ErbB1_'+erb.name+'_noSOS', cpp_free + dimer_free <> cpp_bound % dimer_bound, k4b, kd4)

        # We'll need to save this pattern to build the next two. RAS.state
        # conflicts with ErbB*.state so we can't incrementally build them by
        # overwriting the same variable.
        dimer_free_noras = dimer_free(sos=2) % SOS(grb=2, ras=None)
        dimer_bound = make_dimer_bound(dimer_free_noras)
        Rule(place+'_bind_cPP_ErbB1_'+erb.name+'_noRAS', cpp_free + dimer_free_noras <> cpp_bound % dimer_bound, k4b, kd4)

        dimer_free = dimer_free_noras(ras=3) % RAS(sos=3, state='gdp')
        dimer_bound = make_dimer_bound(dimer_free)
        Rule(place+'_bind_cPP_ErbB1_'+erb.name+'_RAS_GDP', cpp_free + dimer_free <> cpp_bound % dimer_bound, k4b, kd4)

        # Finally, this is the special case.
        dimer_free = dimer_free_noras(ras=3) % RAS(sos=3, state='gtp')
        dimer_bound = make_dimer_bound(dimer_free)
        Rule(place+'_bind_cPP_ErbB1_'+erb.name+'_RAS_GTP', cpp_free + dimer_free <> cpp_bound % dimer_bound, k4, kd4)

    # ErbB2:Erb(2/3/4) : Grb2(){compartment = ANY) + cPP <-> ErbB2:ErbB2/3/4:Grb2(): CPP { compartment = plasma membrane}
    for erb in receptors[1:]:
        Rule(place+'_bind_cPP_ErbB2_'+erb.name, cPP(erb=None, comp=place) + Grb2(gab1=None) % ErbB2(d=1, cpp=None, gap=ANY, comp=place) %
         erb(d=1,gap=None, comp=place) <> cPP(erb=2, comp='pm') % Grb2(gab1=None) % ErbB2(d=1,cpp=2, gap=ANY, comp='pm') %
         erb(d=1,gap=None, comp='pm'), k5, kd4)
        
    place = 'endo'
    # EGF:ErbB1 dimers require separate rules since bound ligand is also transported to pm compartment
    # 2(EGF:ErbB1):Grb2(){compartment = ANY) + cPP <-> 2(EGF:ErbB1):Grb2(): CPP { compartment = plasma membrane}
    Rule(place+'_bind_cPP_ErbB1_ErbB1_noSOS', cPP(erb=None,comp=place) + Grb2(gab1=None, sos=None) % EGF(rec=3, comp=place) %
          ErbB1(lig=3, d=1, cpp=None, gap=ANY, comp=place) % EGF(rec=4, comp=place) % ErbB1(lig=4, d=1,gap=None, comp=place) <>
          cPP(erb=2, comp='pm') % Grb2(gab1=None, sos=None) % EGF(rec=3, comp='pm') % ErbB1(lig=3, d=1,cpp=2, gap=ANY, comp='pm') %
          EGF(rec=4, comp='pm') % ErbB1(lig=4, d=1,gap=None, comp='pm'), k5, kd5)
     
    Rule(place+'_bind_cPP_ErbB1_ErbB1', cPP(erb=None,comp=place) + SOS(erk=None, state='up')% Grb2(gab1=None, sos=ANY) % EGF(rec=3, comp=place) %
          ErbB1(lig=3, d=1, cpp=None, gap=ANY, comp=place) % EGF(rec=4, comp=place) % ErbB1(lig=4, d=1,gap=None, comp=place) <>
          cPP(erb=2, comp='pm') % SOS(erk=None, state='up') % Grb2(gab1=None, sos=ANY) % EGF(rec=3, comp='pm') % ErbB1(lig=3, d=1,cpp=2, gap=ANY, comp='pm') %
          EGF(rec=4, comp='pm') % ErbB1(lig=4, d=1,gap=None, comp='pm'), k5, kd5)
         
    # ErbB1:Erb(2/3/4) : Grb2(){compartment = ANY) + cPP <-> ErbB1:ErbB2/3/4:Grb2(): CPP { compartment = plasma membrane}
    for erb in receptors[1:]:
             Rule(place+'_bind_cPP_ErbB1_'+erb.name, cPP(erb=None,comp=place) + Grb2(gab1=None) % ErbB1(d=1, cpp=None, gap=ANY, comp=place) %
                  erb(d=1,gap=None, comp=place) <> cPP(erb=2, comp='pm') % Grb2(gab1=None) % ErbB1(d=1,cpp=2, gap=ANY, comp='pm') %
                  erb(d=1,gap=None, comp='pm'), k5, kd5b)

    # ErbB2:Erb(2/3/4) : Grb2(){compartment = ANY) + cPP <-> ErbB2:ErbB2/3/4:Grb2(): CPP { compartment = plasma membrane}
    for erb in receptors[1:]:
        Rule(place+'_bind_cPP_ErbB2_'+erb.name, cPP(erb=None, comp=place) + Grb2(gab1=None) % ErbB2(d=1, cpp=None, gap=ANY, comp=place) %
         erb(d=1,gap=None, comp=place) <> cPP(erb=2, comp='pm') % Grb2(gab1=None) % ErbB2(d=1,cpp=2, gap=ANY, comp='pm') %
         erb(d=1,gap=None, comp='pm'), k5, kd5b)

#########################
def R_deg_v2():
    Parameter('k61', 5.7e-4)
    Parameter('k60', 0.00266742)
    Parameter('k60b', 0.0471248)
    Parameter('k60c', 5.2e-4)
    Parameter('k62b', 4.16e-4)
    Parameter('k116', 0.0150356)
    alias_model_components()
    
    " The rules below tanslate to 92 reaction, instead of the 80 in the model. The remaining 12 appear to be missing degradation rxns in clusters 34, 35, 36, and 38 (Refer JHM diagram) "
    
    # degrade RAS bound receptor dimers in endo comaprtment
    Rule('degrade_ras_bound_ErbB1_homodimers', ErbB1(d=1, gap=ANY, comp='endo') % ErbB1(d=1, gap=None, comp='endo') % RAS() >> None, k60)
    
    for erb in receptors[1:]:
        Rule('degrade_ras_bound_ErbB1_heterodimers_'+erb.name, ErbB1(d=1, gap=ANY, comp='endo') %
             erb(d=1, gap=None, comp='endo') % RAS()>> None, k60b)
    
    Rule('degrade_ras_bound_ErbB2_homodimers', ErbB2(d=1, gap=ANY, comp='endo') % ErbB2(d=1, gap=None, comp='endo') % RAS() >> None, k60b)
    
    for erb in receptors[2:]:
        Rule('degrade_ras_bound_ErbB2_heterodimers_'+erb.name, ErbB2(d=1, gap=ANY, comp='endo') %
             erb(d=1, gap=None, comp='endo')% RAS() >> None, k60c)
    
    # degrade Grb2 bound dimers in endo compartment
    # Frustratingly, the reactions involving ErbB1 homodimers use one rate,
    # ErbB1 heterodimers and ErbB2 homodimers another rate, and the rest a
    # third rate. This structure is an attempt at implementing this scheme
    # without repeating ourselves too much.
    for erb in receptors[:2]:
        names_patterns = (
            ('degrade_grb2noSos_bound_dimers',
             erb(comp='endo', gs=1) % Grb2(erb=1, shc=None, sos=None)),
            ('degrade_grb2_bound_dimers',
             erb(comp='endo', cpp=None, gs=1) % Grb2(erb=1, shc=None, sos=2) %
             SOS(grb=2, ras=None, state='up', erk=None)),
            ('degrade_shc_p_grb2noSos_bound_dimers',
             erb(comp='endo', gs=1) % SHC(erb=1, grb=2) % Grb2(shc=2, sos=None))
            )
        for name, pattern in names_patterns:
            for other_erb in receptors:
                rule_name = name + '_' + erb.name + '_' + other_erb.name
                pattern_dimer = pattern(d=99) % other_erb(d=99)
                rate = k60c
                if erb is ErbB1 and other_erb is ErbB1:
                    rate = k60
                elif ((erb is ErbB1 and other_erb is not ErbB1) or
                      (erb is ErbB2 and other_erb is ErbB2)):
                    rate = k60b
                Rule(rule_name, pattern_dimer >> None, rate)
    
    # The next section has degradation reactions missing in the sbml model
    # ====================================================================
    # Degrade endo receptors of the form Shc~P:Grb2:SOS. Missing component, ErbB2 homodimers
    for erb in receptors:
        rate = k60b
        if erb is ErbB1:
            rate = k60
        Rule('degrade_ErbB1_'+erb.name+'_shcpgrb2sos', ErbB1(d=1, comp='endo', gap=ANY) % erb(d=1, comp='endo', gap=None) %
             Grb2(erb=None, shc=ANY, sos=ANY)% SOS(ras=None, state='up', erk=None) >> None, rate)
    for erb in receptors[2:]:
        Rule('degrade_ErbB2_'+erb.name+'_shcpgrb2sos', ErbB2(d=1, comp='endo', gap=ANY) % erb(d=1, comp='endo', gap=None) %
             Grb2(erb=None, shc=ANY, sos=ANY)% SOS(ras=None, state='up', erk=None) >> None, k60c)

    # Degrade endo receptors of the form Shc~P / Shc. Missing component, ErbB1 heterodimers
    Rule('degrade_ErbB1_ErbB1_shc', ErbB1(d=1, comp='endo', gap=ANY) % ErbB1(d=1, comp='endo', gap=None) %
         SHC(erb=ANY, grb=None) >> None, k60)
    for erb in receptors[1:]:
             rate = k60c
             if erb is ErbB2:
                 rate = k60b
             Rule('degrade_ErbB2_'+erb.name+'_shc', ErbB2(d=1, comp='endo', gap=ANY) % erb(d=1, comp='endo', gap=None) %
                  SHC(erb=ANY, grb=None)  >> None, rate)

    # Degrade endo receptors fo the form GAP. Only homodimers are degraded
    for erb in receptors[:2]:
        rate = k60b
        if erb is ErbB1:
            rate = k60
        Rule('degrade_gapbound_'+erb.name+'homodimers', erb(d=1, comp='endo', gap=ANY, gs=None) %
         erb(d=1, comp='endo', gap=None, gs=None) % GAP() >> None, rate)
         # ====================================================================================================
         
    # Degrade 2(EGF:ErbB1:ATP) v524
    Rule('degrade_ATP_bound_ErbB1_homodimers', ErbB1(d=1, gap=None, atp=ANY,comp='endo', state='up', rtk=None)%
              ErbB1(d=1, gap=None,  atp=ANY,comp='endo', state='up', rtk=None) >> None, k60)
         
    # degrade receptor dimers in endo that are not bound to ATP, rtk or phosphorylated
    for erb in receptors[1:]:
             Rule('degrade_noATP_bound_ErbB1_'+ erb.name, ErbB1(d=1, gap=None, comp='endo', atp=None, state='up', rtk=None) %
                  erb(d=1, gap=None, comp='endo', atp=None, state='up', rtk=None) >> None, k62b)
         
    for erb in receptors[1:]:
        rate = k62b
        if erb is ErbB2:
            rate = k60b
        Rule('degrade_noATP_bound_ErbB2_'+ erb.name, ErbB2(d=1, gap=None, comp='endo', atp=None, state='up', rtk=None) %
             erb(d=1, gap=None, comp='endo', atp=None, state='up', rtk=None) >> None, rate)
    
    
    # degrade ErbB monomers that are not phosphoryated or bound to ligand
    
    Rule('degrade_ErbB1_monomer', ErbB1(d=None, lig=None, atp=ANY, comp='endo', state='up') >> None, k60)
    
    for erb in receptors[1:]:
        Rule('degrade_monomers_' +erb.name, erb(d=None, lig=None, comp='endo', state='up') >> None, k60b)
    
    # Degrade Pase3 (V769)
    #Rule('degrade_Pase3', Pase3(erk=None) >> None, k116)
    
    Rule('degrade_endo_EGF', EGF(comp='endo', rec=None) >> None, k61)

def declare_observables():
    alias_model_components()
    
    # Original model missed the ErbB1/2..Gab1#P#P species when enumerating the
    # phosphorylated ErbB1 dimer species. We'll account for that by creating an
    # observable with the correct pattern, observables for the missed species,
    # and an expression with the difference of the two.
    Observable('pErbB1_total', ErbB1(d=ANY, state='p'))
    erbb1p_grb2_gab1pp = (ErbB1(state='p', gs=1) % Grb2(erb=1, gab1=2) %
                          Gab1(grb2=2, state='pp'))
    Observable('pErbB11_exceptions', erbb1p_grb2_gab1pp(d=3) % ErbB1(d=3))
    Observable('pErbB12_exceptions', erbb1p_grb2_gab1pp(d=3) % ErbB2(d=3))
    alias_model_components()
    # ErbB1:1 observable pattern isn't symmetric, so we need to explicitly
    # multiply it by 2 to account for the two ErbB1 molecules in those species.
    Expression('pErbB1',
               pErbB1_total - pErbB11_exceptions * 2 - pErbB12_exceptions)

