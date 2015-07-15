""" Model for ErbB receptor activation """
# enforce no free SOS:RAS
# compartments added

from pysb import*
from pysb.util import alias_model_components
import pysb.core
import pysb.macros
from pysb.macros import bind_table
from pysb.macros import catalyze
from pysb.macros import catalyze_state

##################
# Monomers of the model
def declare_monomers():
    global receptors
    global comprtmnts
    
    # ErbB receptors bind ligands at site lig, other ligands at side d, ATP at side atp, GAP on gs, SHC or Grb2 on gs, RTK on rtk, CPP on cpp
    Monomer('Erb1', ['lig', 'd', 'atp', 'state', 'gap', 'gs', 'rtk', 'cpp', 'comp'], {'state':['up', 'p', 'pm', 'endo'], 'comp': ['pm', 'endo']})
    Monomer('Erb2', ['lig', 'd', 'atp', 'state', 'gap', 'gs', 'rtk', 'cpp', 'comp'], {'state':['up', 'p'], 'comp': ['pm', 'endo']})
    Monomer('Erb3', ['lig', 'd', 'atp', 'state', 'gap', 'gs', 'rtk', 'cpp','comp'], {'state':['up', 'p'], 'comp': ['pm', 'endo']})
    Monomer('Erb4', ['lig', 'd', 'atp', 'state', 'gap', 'gs', 'rtk', 'cpp','comp'], {'state':['up', 'p'], 'comp': ['pm', 'endo']})
    # EGF binds ErbB1 on site lig
    Monomer('EGF', ['rec', 'comp'], {'comp': ['pm', 'endo']})
    # HRG binds ErbB3/4 on site lig
    Monomer('HRG', ['rec', 'comp'], {'comp': ['pm', 'endo']})
    Monomer('ATP', ['erb','gab1'])
    Monomer('RTK', ['erb'])

    Monomer('GAP', ['rec'])
    Monomer('Grb2', ['shc', 'erb', 'sos', 'gab1'])
    Monomer('SOS', ['grb', 'ras', 'erk','state'],{'state':['up','p']})
    Monomer('SHC', ['erb','grb', 'state'], {'state':['up','p']})
    Monomer('RAS', ['sos', 'pi3k', 'raf', 'state'], {'state':['gdp', 'gtp', 'active_gtp']})

    Monomer('cPP', ['erb', 'comp'], {'comp':['pm', 'endo']})

    Monomer('Gab1',['atp', 'grb2', 'shp2', 'state','erk','pase', 'pi3k'], {'state':['up','p', 'pp']})
    Monomer('Shp2',['gab1'])

    Monomer('Pase_9t', ['gab1'])
    Monomer('PI3K', ['gab1', 'pip2', 'ras'])
    Monomer('PIP2', ['pi3k'])
    Monomer('PIP3', ['akt', 'pdk','bnd'])
    Monomer('AKT', ['pip', 'pase', 'raf', 'state'], {'state':['up','p', 'pp']})
    Monomer('PDK1',['pip'])
    Monomer('Pase4', ['akt'])
    Monomer('Shp', ['pip'])
    Monomer('PTEN', ['pip'])
    Monomer('RAF', ['akt', 'pase1', 'mek','state', 'ras'], {'state':['up','p', 'p_ser']})
    Monomer('MEK', ['raf', 'pase2', 'erk', 'state'], {'state':['up','p', 'pp']})
    Monomer('ERK', ['mek', 'pase3', 'gab1','sos', 'state'], {'state':['up','p', 'pp']})
    Monomer('Pase1', ['raf'])
    Monomer('Pase2', ['mek'])
    Monomer('Pase3', ['erk'])
##
    alias_model_components()
    
    receptors = [Erb1, Erb2, Erb3, Erb4]
    comprtmnts = ['pm', 'endo']

##############################
# Non-zero initial conditions (in molecules per cell)
def declare_initial_conditions():
    Parameter('Erb1_0',	1.08e6)     # c531
    Parameter('Erb2_0', 4.62e5)     # c141
    Parameter('Erb3_0', 6230)       # 140
    Parameter('Erb4_0', 794)        # c143
    Parameter('EGF_0', 5e-9)        # c1
    Parameter('HRG_0', 0)           # c514 zero ic
    Parameter('ATP_0', 1.2e9)       # c105
    Parameter('RTK_0', 7e+4)        # c280

    Parameter('GAP_0', 534751)
    Parameter('Grb2_0', 1264.91)    # c22 //// Grb2:SOS has non-zero initial condition (c30)
    Parameter('SOS_0', 0)           # c24 // zero ic
    Parameter('Grb2SOS_0', 8.89e7)  # c30
    Parameter('SHC_0', 11e+5)       # c31
    Parameter('RAS_gdp_0', 58095.2) # c26
    Parameter('dummy_init', 1e-10)  # to initialize RAS-GDP and RAS_active GTP
    
    Parameter('cPP_0', 4498.73)     # c12

    Parameter('Gab1_0',94868.3)     # c426
    Parameter('Shp2_0', 1e+6)       # c463
    Parameter('Pase9t_0', 0)        # c521 // zero ic
    Parameter('PI3K_0', 3.55e+7)    # c287
    Parameter('PIP2_0', 393639)     # c444
    Parameter('PIP3_0', 0)          # c106 // zero ic
    Parameter('AKT_0', 905000)      # c107
    Parameter('PDK1_0', 3e8)        # c109
    Parameter('Pase4_0', 14.5e+5)   # c113
    Parameter('PTEN_0', 56100.9)   # c279
    Parameter('Shp_0', 2213.59)     # c461
    Parameter('RAF_0', 71131.2)     # c41
    Parameter('MEK_0', 3020000)     # c47
    Parameter('ERK_0', 695000)      # c55
    Parameter('Pase1_0', 5e+4)      # c44
    Parameter('Pase2_0', 124480)    # c53
    Parameter('Pase3_0', 16870.2)   # c60

    alias_model_components()
    
    Initial(Erb1(lig=None, d=None, atp=None, gap=None, gs=None, rtk=None, cpp=None, state='up', comp ='pm'), Erb1_0)
    Initial(Erb2(lig=None, d=None, atp=None, gap=None, gs=None, rtk=None, cpp=None, state='up', comp ='pm'), Erb2_0)
    Initial(Erb3(lig=None, atp=None, d=None, gap=None, gs=None, rtk=None,cpp=None, state='up', comp ='pm'), Erb3_0)
    Initial(Erb4(lig=None, atp=None, d=None, gap=None, gs=None, rtk=None,cpp=None, state='up', comp ='pm'), Erb4_0)
    Initial(EGF(rec=None, comp ='pm'), EGF_0)
    Initial(HRG(rec=None, comp ='pm'), HRG_0)
    Initial(ATP(erb=None, gab1=None), ATP_0)
    Initial(RTK(erb=None), RTK_0)
#
    Initial(GAP(rec=None), GAP_0)
    Initial(Grb2(shc=None, erb=None, sos=None, gab1=None), Grb2_0)
    Initial(SOS(grb=None, ras=None, erk=None, state='up'), SOS_0)
    Initial(SHC(erb=None, grb=None, state='up'), SHC_0)
    Initial(RAS(sos=None, pi3k=None, raf=None,state='gdp'), RAS_gdp_0)
    # Initialized temporarily to check topology
    #Initial(RAS(sos=None, pi3k=None, raf=None, state='gtp'), dummy_init)
    #Initial(RAS(sos=None, pi3k=None, raf=None,state='active_gtp'), dummy_init)
    #####
    Initial(cPP(erb=None, comp='pm'), cPP_0)

    Initial(Gab1(atp=None, grb2=None, shp2=None, erk=None, pase=None, pi3k=None, state='up'), Gab1_0)
    Initial(Shp2(gab1=None), Shp2_0)

    Initial(Pase_9t(gab1=None), Pase9t_0)
    Initial(PI3K(gab1=None, pip2=None, ras=None), PI3K_0)
    Initial(PIP2(pi3k=None), PIP2_0)
    Initial(PIP3(akt=None, pdk=None, bnd=None), PIP3_0)
    #  Initial(HRG(rec=None, comp ='endo'), dummy_init)
    Initial(AKT(pip=None, pase=None,raf=None, state='up'), AKT_0)
    Initial(PDK1(pip=None), PDK1_0)
    Initial(Pase4(akt=None), Pase4_0)
    Initial(PTEN(pip=None), PTEN_0)
    Initial(Shp(pip=None), Shp_0)
    Initial(RAF(akt=None, pase1=None, mek=None, ras=None, state='p'), RAF_0)
    Initial(MEK(raf=None, pase2=None, erk=None,state='up'), MEK_0)
    Initial(ERK(mek=None, pase3=None, sos=None, state='up', gab1=None), ERK_0)
    #Initial(ERK(mek=None, pase3=None,  state='pp', gab1=None), dummy_init)

    Initial(Pase1(raf=None), Pase1_0)
    Initial(Pase2(mek=None), Pase2_0)
    Initial(Pase3(erk=None), Pase3_0)

    Initial(Grb2(shc=None, erb=None, sos=1, gab1=None) % SOS(grb=1, ras=None, erk=None,state='up'), Grb2SOS_0)


########################################################
def ErbB1_priming():
    " v828 "
    " Unlike other Erb monomers, ErbB1 can participate in reactions only after binding with ATP "
    
    Parameter('kf_erb1_atp', 1.87e-8)  # k122
    Parameter('kr_erb1_atp', 1)        # kd122
    
    alias_model_components()
    
    # ErbB1 + ATP <-> ErbB1:ATP
    Rule('Erb1_atp', Erb1(lig=None,d=None,atp=None, comp='pm') + ATP(erb=None,gab1=None) <>
         Erb1(lig=None,d=None,atp=1, comp='pm') % ATP(erb=1,gab1=None), kf_erb1_atp, kr_erb1_atp)

########################################################
def transport():
    
    " Free ErbB monomers (v164,v176-178), phosphorylated dimers(v165-175,179-193) and cPP(v211) can shuttle between plasma membrane and endosomol membrane "
    " *** Need to check for other transport reactions "
    Parameter('kf_cPP_pm_endo', 0)              # k15
    Parameter('kr_cPP_pm_endo', 1.667e-8)       # kd15
    Parameter('k7', 5e-05)                      # k7
    Parameter('kd7', 1.38e-4)                   # kd7
    Parameter('k6', 0.013)                      # k6; k6b = 0
    Parameter('kd6', 5e-05)                     # kd6; kd6b = 0
    
    alias_model_components()
    
    # Transport Erb1:ATP v164
    Rule('Transport_Erb1_ATP', Erb1(d=None, lig=None, atp=ANY, state='up', comp='pm') <>
         Erb1(d=None, lig=None, atp=ANY, state='up', comp='endo'), k6, kd6)
    
    " V176-178 and v181-184 have rate constants set to zero in both directions "
    # Transport ErbB 2,3,4 monomers
#    for erb in receptors[1:]:
#        Rule('transport_' + erb.name, erb(lig=None, d=None, atp=None, state='up', comp='pm') <>
#             erb(lig=None, d=None, atp=None, state='up', comp='endo'), k6b, kd6b) # rate constants are k6b = kd6b = 0

    " v165 - 175 "
    Rule('Erb11_pm_to_endo', EGF(rec=2, comp='pm') % EGF(rec=3, comp='pm') % Erb1(lig=2, d=1, rtk=None, cpp=None, gap=None, state='p', comp='pm')
         % Erb1(lig=3, d=1, rtk=None, state='p', gap=None, comp='pm') <>
         EGF(rec=2, comp='endo') % EGF(rec=3, comp='endo') % Erb1(lig=2, d=1, rtk=None, cpp=None, gap=None, state='p', comp='endo') %
         Erb1(lig=3, d=1, rtk=None, state='p', gap=None, comp='endo'), k6, kd6)

    " v185-187 "
    for erb in receptors[1:]:
        Rule('Erb1_pm_to_endo'+erb.name, Erb1(d=1, rtk=None, cpp=None, gap=None, state='p', comp='pm') % erb(d=1, rtk=None, state='p', cpp=None,gap=None, comp='pm') <>
             Erb1(d=1, rtk=None, cpp=None, state='p', gap=None, comp='endo') % erb(d=1, rtk=None, state='p', cpp=None, gap=None, comp='endo'),
             k7, kd7)

    " v188-190 "
    for erb in receptors[1:]:
        Rule('Erb2_pm_to_endo'+erb.name, Erb2(d=1, rtk=None,atp=None, state='p', gap=None, comp='pm') % erb(d=1, atp=None, rtk=None, state='p', gap=None, comp='pm') <>
        Erb2(d=1, atp=None, rtk=None,state='p', gap=None, comp='endo') % erb(d=1, atp=None, rtk=None, state='p', gap=None, comp='endo'),
        k7, kd7)

    Rule('transport_cPP', cPP(erb=None, comp='pm') <> cPP(erb=None, comp='endo'), kf_cPP_pm_endo, kr_cPP_pm_endo)

    " v180, v191-193 have non-zero rate constants, and still need to be implimented "

########################################################
def ligand_binding():
    
    " EGF binds ErbB1 (v1), HRG bind ErbB3 and ErbB4 (v784-785)"
    
    Parameter('EGF_bind', 1e+7)            # k1
    Parameter('HRG_bind', 1e+7)            # k119
    Parameter('EGF_unbind', 0.033)      # kd1
    Parameter('HRG_unbind', 0.0103115)  # kd119
        
    alias_model_components()
    
    for place in comprtmnts:
        # ErbB1:ATP + EGF -> EGF:Erb1:ATP
        Rule('Erb1_EGF_'+place, Erb1(lig=None,d=None,atp=ANY, state='up', comp=place) + EGF(rec=None, comp=place) <>
             Erb1(lig=1,d=None,atp=ANY, state='up', comp=place) % EGF(rec=1, comp=place), EGF_bind, EGF_unbind)
        
        # ErbB3 + HRG -> ErbB3:HRG
        Rule('Erb3_HRG' + place, Erb3(lig=None,d=None,atp=None,state='up', comp=place) + HRG(rec=None, comp=place) <>
             Erb3(lig=1,d=None,atp=None, state='up', comp=place) % HRG(rec=1, comp=place), HRG_bind, HRG_unbind)
        
        # ErbB4 + HRG -> ErbB4:HRG
        Rule('Erb4_HRG'+place, Erb4(lig=None,d=None,atp=None, state='up', comp=place) + HRG(rec=None, comp=place) <>
             Erb4(lig=1,d=None,atp=None, state='up', comp=place) % HRG(rec=1, comp=place), HRG_bind, HRG_unbind)

########################################################
def receptor_dimerization():
    
    " Heterodimers form when a ligand bound monomer binds a naked monomer, as shown below "
    " Lig:monomer:ATP + monomer <-> Lig:dimer *** Note that ATP is missing from the right hand side product "

    " ErbB1 homodimers form when two ligand bound, ATP bound monomers dimerize. ErbB2 homodimers cannot be formed by this route \
    since they do not bind ligand, instead by secondary dimerization described in a separate module "
    
    Parameter('k2', 7.4e-6)
    Parameter('kd2', 0.16)
    Parameter('k2b', 3.73e-8)
    Parameter('kd2b', 0.016)
    Parameter('k120', 1.48e-8)
    Parameter('kd120', 0.1)
    Parameter('k120b', 5.9e-11)
    Parameter('kd120b', 0.1)
    
    alias_model_components()
    
    for place in comprtmnts:
    
        # EGF:Erb1:ATP + Erb2/3/4 <-> Erb1:Erb2/3/4
        for erb in receptors[1:]:
            Rule(place+'_bind_egf_E1_atp_' + erb.name, ATP(erb=2,gab1=None) % Erb1(lig=ANY,atp=2,d=None, state='up',comp=place) +
                 erb(lig=None,d=None,atp=None, state='up', comp=place) <>
                 Erb1(lig=ANY, atp=None, d=1, state='up', comp=place) % erb(lig=None,d=1, atp=None, state='up', comp=place),
                 k2b, kd2b)

        # EGf:Erb1:ATP + EGF:Erb1:ATP <-> 2(EGF:Erb1)
        Rule(place+'_bind_egf_E1_atp_egf_erb1', ATP(erb=2,gab1=None) % Erb1(lig=ANY,atp=2,d=None,state='up', comp=place) +
             ATP(erb=3,gab1=None) % Erb1(lig=ANY, d=None,atp=3,state='up', comp=place) <>
             Erb1(lig=ANY, atp=None, d=1,state='up', comp=place) % Erb1(lig=ANY, d=1, atp=None,state='up', comp=place) +
             ATP(erb=None,gab1=None) , k2, kd2)
    
        # HRG:Erb(3/4) + Erb2 <-> Erb2:Erb(3/4)
        bind_table([[   Erb2(atp=None, state='up', comp=place)],
                [Erb3(lig=ANY,atp=None, comp=place), (k120, kd120)],
                [Erb4(lig=ANY,atp=None, comp=place), (k120, kd120)]],'d', 'd')

        # HRG:Erb(3/4) + Erb1:ATP <-> Erb(3/4):Erb1
        Rule(place+'_E1_atp__H_E3', Erb3(lig=ANY, atp=None, d=None, comp=place) + ATP(erb=2,gab1=None) % Erb1(lig=None,atp=2, d=None, comp=place)
             <> Erb3(lig=ANY,atp=None, d=1, comp=place) % Erb1(lig=None, atp=None, d=1, comp=place), k120b, kd120b)

        Rule(place+'_E1_atp__H_E4', Erb4(lig=ANY, atp=None, d=None, comp=place) + ATP(erb=2,gab1=None) % Erb1(lig=None,atp=2, d=None,comp=place)
             <> Erb4(lig=ANY,atp=None, d=1, comp=place) % Erb1(lig=None, atp=None, d=1, comp=place), k120b, kd120b)

########################################################
def secondary_dimerization():
    
    " v673-v681, When phosphorylated dimers dissociate, dissociation products can pair up in new combinations. Erb2 homodimers are formed as a\
    result of one such secondary dimerization process "
    
    Parameter('k102', 5e-07)
    Parameter('kd102', 5.6)
    Parameter('k103', 8.3e-9)
    Parameter('kd103', 0.016)
    
    alias_model_components()
    
    " secondary dimerization does not take place in endo compartment. DO NOT loop over compartments"
        #for place in comprtmnts:
    
    # 2(EGF:Erb1)~P <-> EGF:Erb1
    Rule('Unbind_EGF_Erb1_P_', Erb1(lig=ANY, state='p',atp=None,d=1,gap=None,rtk=None, comp='pm')
             % Erb1(lig=ANY, state='p',atp=None,d=1,gap=None, rtk=None, comp='pm') <>
             Erb1(lig=ANY, state='p',atp=None,d=None,gap=None, rtk=None, comp='pm') +
             Erb1(lig=ANY, state='p',atp=None,d=None, gap=None, rtk=None, comp='pm'), k102, kd102)
    
    # EGF:Erb1~P + Erb2/3/4~P <-> (Erb1:Erb2/3/4)~P
    for erb in receptors[1:]:
            Rule('pbind_erb1_'+ erb.name, EGF(rec=2, comp='pm') % Erb1(lig=2, state='p', d=None, atp=None,gap=None, rtk=None, comp='pm') +
                 erb(lig=None, state='p', d=None, atp=None,gap=None, rtk=None, comp='pm') <>
                 Erb1(lig=None, state='p', d=1,atp=None,gap=None, rtk=None, comp='pm') %
                 erb(lig=None, state='p', d=1, atp=None,gap=None, rtk=None, comp='pm'), k102, kd102)
    
    # Erb2~P + Erb(1:4)~P <-> (Erb2:Erb1/2/3/4)~P
    bind_table([[	Erb1(lig=ANY,state='p',gap=None,rtk=None, comp='pm'),	Erb2(state='p',gap=None, rtk=None, comp='pm'),	Erb3(lig=None,state='p',gap=None, rtk=None, comp='pm'),	Erb4(lig=None,state='p',gap=None, rtk=None, comp='pm')],
                [Erb2(state='p',gap=None, rtk=None, comp='pm'),	None,	(k102, kd102),	(k103, kd103),	(k103, kd103)]],'d', 'd')

########################################################
def lateral_signaling():
    
    Parameter('bind_Erb2_Erb2up_kf', 8.3e-9)    # k103
    Parameter('bind_Erb2_Erb2up_kr', 0.016)     # kd103
    
    alias_model_components()
    
    for place in comprtmnts:
    
    # Erb2~P + Erb2 <-> Erb2~P:Erb2
        Rule('bind_Erb2_Erb2_up_'+place, Erb2(d=None, state='p', gap=None, rtk=None, comp=place) +
             Erb2(d=None, state='up',gap=None, rtk=None, comp=place) <>
             Erb2(d=1, state='p',gap=None, rtk=None, comp=place) % Erb2(d=1,state='up',gap=None, rtk=None, comp=place),
             bind_Erb2_Erb2up_kf,bind_Erb2_Erb2up_kr)

########################################################
def trans_phosphorylation():
    "v21-v29, v37-v41, v806-814, v822-827 "
    " Lig:dimer + ATP <-> Lig:dimer:ATP -> dimer~P + ATP *** Note that ligand dissapears in the ES -> E + P step "
    " check v29 since product end sup in the endosome"
    
    Parameter('kon', 1.87e-8)           # k122
    Parameter('koff', 1)                # kd122
    Parameter('kcat_phos', 0.177828)    # kd123
    Parameter('krcat_phos', 0)          # k123
    
    alias_model_components()
    
    for place in comprtmnts:
    
        # HRG:Erb(3/4):ERB(1/2) -> Erb(3/4)~P:Erb(1/2)~P
        for s in receptors[2:]:
            for r in receptors[:2]:
                Rule(place + '_' + s.name+'_tp_bind_'+ r.name, s(lig=ANY,d=1, state='up', atp=None, comp=place) %
                     r(lig=None,d=1, state='up', atp=None) + ATP(erb=None,gab1=None) <>
                     s(lig=ANY,d=1, state='up', atp=2, comp=place) % r(lig=None,d=1, state='up', atp=None) % ATP(erb=2, gab1=None), kon, koff)
                     
                Rule(place + '_' + s.name+'_tp_cat_'+ r.name, HRG(rec=3, comp=place) % s(lig=3,d=1, state='up', atp=2, comp=place, gap=None) %
                     r(lig=None,d=1, state='up', atp=None, gap=None) % ATP(erb=2,gab1=None) >>
                     s(lig=None,d=1, state='p', atp=None, comp=place, gap=None) % r(lig=None,d=1, state='p', atp=None, gap=None) +
                     ATP(erb=None,gab1=None), kcat_phos)

        # EGF:Erb1::Erb(2/3/4) -> Erb1~P:Erb(2/3/4)~P
        for r in receptors[1:]:
            Rule(place+'_Erb1_tp_bind_'+ r.name, Erb1(lig=ANY,d=1, state='up', atp=None, comp=place) %
                 r(lig=None,d=1, state='up', atp=None) + ATP(erb=None,gab1=None) <>
                 Erb1(lig=ANY,d=1, state='up', atp=2, comp=place) % r(lig=None,d=1, state='up', atp=None)% ATP(erb=2,gab1=None), kon, koff)
                 
            Rule(place +'_Erb1_tp_cat_'+ r.name, EGF(rec=3)%Erb1(lig=3,d=1, state='up', atp=2, comp=place) %
                 r(lig=None,d=1, state='up', atp=None)% ATP(erb=2,gab1=None) >>
                 Erb1(lig=None,d=1, state='p', atp=None, comp=place) % r(lig=None,d=1, state='p', atp=None) + ATP(erb=None,gab1=None), kcat_phos)

        # 2(EGF:Erb1) -> 2(EGF:Erb1)~P
        Rule('Erb1_tp_bind_Erb1_'+ place, Erb1(lig=ANY,d=1, state='up', atp=3, comp=place)% ATP(erb=3, gab1=None) %
             Erb1(lig=ANY,d=1, state='up', atp=None) + ATP(erb=None, gab1=None) <>
             ATP(erb=3, gab1=None) % ATP(erb=2,gab1=None) % Erb1(lig=ANY,d=1, state='up', atp=3, comp=place) % Erb1(lig=ANY,d=1, state='up', atp=2) , kon, koff)
            
        Rule('Erb1_tp_cat_Erb1_'+ place, Erb1(lig=ANY,d=1, state='up', atp=2, comp=place) %
             Erb1(lig=ANY,d=1, state='up', atp=3)% ATP(erb=2,gab1=None)% ATP(erb=3, gab1=None) >>
             Erb1(lig=ANY,d=1, state='p', atp=None, comp=place) % Erb1(lig=ANY,d=1, state='p', atp=None) + ATP(erb=None,gab1=None) , kcat_phos)

        # Erb2~P:Erb2 -> 2(Erb2)~P
        Rule('Erb2_tp_Erb2_up_'+ place, Erb2(d=1, state='p', comp=place) % Erb2(d=1, state='up') >>
             Erb2(d=1, state='p', comp=place) % Erb2(d=1, state='p'), kcat_phos)

        # Erb3/4~P:Erb2~P -> HRG:Erb3/4:Erb2:ATP

        for s in receptors[2:]:
            Rule(place+'krcat_Erb2_'+s.name, s(lig=None, d=1, state='p', gap=None, atp=None, rtk=None, comp=place) % Erb2(state='p', d=1, comp=place,rtk=None, cpp=None, gap=None)
                 >> HRG(rec=3, comp=place) % ATP(erb=2, gab1=None) % s(lig=3, d=1, state='up', gap=None, atp=2, rtk=None, comp=place) % Erb2(state='up', d=1, comp=place,rtk=None, cpp=None, gap=None), krcat_phos)

########################################################
def GAP_binding():
    
    "v194-v207"
    " In the Chen 2009 model, each dimer binds to a single unit of downstream substrates like GAP, Grb2, SHC. To reproduce this \
        reaction, GAP (and Grb2, SHc) bind only the Erb1 receptor in all its dimers or Erb2 in its dimers. In Erb1:Erb2 dimers, GAP binds Erb1 "
    
    Parameter('kf_gap', 5.9e-7)   # k8
    Parameter('kr_gap', 0.2)      # kd8
    
    alias_model_components()
    
    for place in comprtmnts:
    
        # Erb1:Erb1/2/3/4 + GAP <->
        for erb in receptors:
            Rule(place+'_gap_bind_erb1_'+erb.name, Erb1(state='p', d=1, gap=None,gs=None, rtk=None, comp=place) %
                 erb(state='p',d=1, gap=None,gs=None, rtk=None) + GAP(rec=None) <>
                 Erb1(state='p', d=1, gap=2,gs=None, rtk=None, comp=place) % erb(state='p',d=1, gap=None,gs=None, rtk=None) %
                 GAP(rec=2), kf_gap, kr_gap)
    
        # Erb2:Erb2/3/4 + GAP <->
        for erb in receptors[1:]:
            Rule(place+'_gap_bind_erb2_'+erb.name, Erb2(state='p', d=1, gap=None,gs=None, rtk=None, comp=place) %
                 erb(state='p',d=1,gap=None,gs=None, rtk=None) + GAP(rec=None) <>
                 Erb2(state='p', d=1, gap=2,gs=None, rtk=None, comp=place) % erb(state='p',d=1,gap=None,gs=None, rtk=None) %
                 GAP(rec=2), kf_gap, kr_gap)

########################################################
def Grb2_binding_v2():
    " Grb2 competes with SHC to bind GAP bound receptors, and dissocaites under all conditions except when SOS is bound to RAS "
    
    Parameter('kf_bind_erb_grb2', 1.67e-05)     # k16
    Parameter('kr_bind_erb_grb2', 0.55)         # kd24 = 0.55 or kd63 = 0.275; kd63 appears in 6 of the 29 backward reactions, no apparent pattern observed
    Parameter('kf_grb2_sos',  1.67e-05)         # k17
    Parameter('kr_grb2_sos', 0.06)              # kd17
    
    alias_model_components()
    
    for place in comprtmnts:
        
        # Grb2 binds Erb1:Erb(1/2/3/4) receptor dimers
        for erb in receptors:
            # Bind
            Rule(place+'_Grb2_bind_erb1_'+erb.name, Erb1(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) + Grb2(erb=None,shc=None,gab1=None, sos=None) <>
                 Grb2(erb=1,shc=None,gab1=None, sos=None) % Erb1(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),kf_bind_erb_grb2, kr_bind_erb_grb2)
                
            Rule(place+'_Grb2sos_bind_erb1_'+erb.name, Erb1(gap=ANY,gs=None,d=2,  cpp=None, comp=place) % erb(d=2,gs=None) + SOS(grb=3, ras=None, erk=None,state='up') % Grb2(erb=None,shc=None,gab1=None, sos=3) <>
                     SOS(grb=3, ras=None,erk=None, state='up') % Grb2(erb=1,shc=None,gab1=None, sos=3) % Erb1(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),kf_bind_erb_grb2, kr_bind_erb_grb2)
    
        # Grb2 binds Erb2:Erb2/3/4 receptor dimers
        for erb in receptors[1:]:
            Rule(place+'_Grb2_bind_erb2_'+ erb.name, Erb2(gap=ANY,gs=None,d=2,  cpp=None, comp=place) % erb(d=2,gs=None) +
                 Grb2(erb=None,shc=None,gab1=None, sos =None) <>
                 Grb2(erb=1,shc=None,gab1=None, sos=None) % Erb2(gap=ANY,gs=1,d=2,  cpp=None, comp=place) % erb(d=2,gs=None),kf_bind_erb_grb2, kr_bind_erb_grb2)
                
            Rule(place+'_Grb2sos_bind_erb2_'+ erb.name, Erb2(gap=ANY,gs=None,d=2,  cpp=None, comp=place) % erb(d=2,gs=None) +
             SOS(grb=3, ras=None,erk=None, state='up') %Grb2(erb=None,shc=None,gab1=None, sos=3) <>
             SOS(grb=3, ras=None, erk=None, state='up') % Grb2(erb=1,shc=None,gab1=None, sos=3) % Erb2(gap=ANY,gs=1,d=2,  cpp=None, comp=place) % erb(d=2,gs=None),kf_bind_erb_grb2, kr_bind_erb_grb2)
             
        # Grb2 binds SOS
        Rule(place+'_bind_free_Grb2_SOS', Grb2(sos=None, erb=None, gab1=None) + SOS(grb=None,ras=None,erk=None, state='up') <>
             Grb2(sos=1, erb=None,gab1=None) % SOS(grb=1, ras=None, erk=None, state='up'), kf_grb2_sos, kr_grb2_sos)
            
        for erb in receptors[:2]:
                 Rule(place+'_bind_Grb2_SOS_'+ erb.name, erb(cpp=None, comp=place, gap=ANY) % Grb2(sos=None,gab1=None) + SOS(grb=None,erk=None, ras=None, state='up') <>
                      erb(cpp=None, gap=ANY, comp=place) % Grb2(sos=1,gab1=None) % SOS(grb=1, ras=None, erk=None, state='up'), kf_grb2_sos, kr_grb2_sos)

        # Release SOS~P from Grb2
    Rule('release_SOS_P', SOS(grb=ANY, state='p') >> SOS(grb=None, state='p'), kr_grb2_sos)

########################################################
def SHC_binding_v2():
    
    " v367-v380"
    " SHC competes with Grb2 to bind GAP primed receptors "
    " SHC can dissociate from receptors under any condition other than when cPP is bound to the receptor, or RAS is bound to SOS, \
    giving rise to dissocation products like SHC, SHC~P. SHC~P:Grb2, SHC~P:Grb2:SOS"
    
    Parameter('kf_bind_erb_shc', 1.39e-7)   # k22 ** Note, separate rate constants k32, k37 for Shc~P binding
    Parameter('kr_bind_erb_shc', 0.1)       # kd22
    Parameter('kshc_kin', 6)                # k23
    Parameter('kshc_phos', 0.06)            # kd23
    
    alias_model_components()
    
    for place in comprtmnts:
        
        # SHC binds Erb1:Erb(1/2/3/4) receptor dimers
        for erb in receptors:
            Rule(place+'_SHC_bind_erb1_'+erb.name, Erb1(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) + SHC(erb=None, grb=None) <>
                 SHC(erb=1, grb=None) % Erb1(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),kf_bind_erb_shc, kr_bind_erb_shc)
    
        # SHC binds Erb2:Erb(2/3/4) receptor dimers
        for erb in receptors[1:]:
            Rule(place+'_SHC_bind_erb2_'+erb.name, Erb2(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) + SHC(erb=None, grb=None) <>
             SHC(erb=1, grb=None) % Erb2(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),kf_bind_erb_shc, kr_bind_erb_shc)
        
        # SHC binds Erb1:Erb(1/2/3/4) receptor dimers
        for erb in receptors:
            Rule(place+'_SHC_grb_bind_erb1_'+erb.name, Erb1(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) + Grb2(shc=3, sos=None) % SHC(erb=None, grb=3) <>
                 Grb2(shc=3, sos=None) % SHC(erb=1, grb=3) % Erb1(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),kf_bind_erb_shc, kr_bind_erb_shc)
        
        # SHC binds Erb2:Erb(2/3/4) receptor dimers
        for erb in receptors[1:]:
            Rule(place+'_SHC_grb_bind_erb2_'+erb.name, Erb2(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) + Grb2(shc=3, sos=None) %SHC(erb=None, grb=3) <>
                 Grb2(shc=3, sos=None) % SHC(erb=1, grb=3) % Erb2(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),kf_bind_erb_shc, kr_bind_erb_shc)
#
#        # SHC binds Erb1:Erb(1/2/3/4) receptor dimers
        for erb in receptors:
            Rule(place+'_SHC_grbsos_bind_erb1_'+erb.name, Erb1(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) + SOS(ras=None,erk=None, state='up') % SHC(erb=None, grb=ANY) <>
         SOS(ras=None, erk=None, state='up') % SHC(erb=1, grb=ANY) % Erb1(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),kf_bind_erb_shc, kr_bind_erb_shc)
        
        # SHC binds Erb2:Erb(2/3/4) receptor dimers
        for erb in receptors[1:]:
            Rule(place+'_SHC_grbsps_bind_erb2_'+erb.name, Erb2(gap=ANY,gs=None,d=2, cpp=None, comp=place) % erb(d=2,gs=None) + SOS(ras=None, erk=None, state='up') % SHC(erb=None, grb=ANY) <>
                 SOS(ras=None, erk=None, state='up') % SHC(erb=1, grb=ANY) % Erb2(gap=ANY,gs=1,d=2, cpp=None, comp=place) % erb(d=2,gs=None),kf_bind_erb_shc, kr_bind_erb_shc)

        # SHC <-> SHC~P
        Rule(place+'_shc_phos', SHC(erb=ANY, grb=None, state='up') <> SHC(erb=ANY, grb=None, state='p'), kshc_kin, kshc_phos)

########################################################
def secondary_Grb2_binding():
    
    " Grb2 binds SHC~P, and dissocaties under all conditions except when RAS is bound to SOS "
    
    Parameter('kf_bind_shc_grb2',  1.67e-05)     # k16
    Parameter('kr_bind_shc_grb2',  0.55)         # kd24 = 0.55 or kd63 = 0.275; kd63 appears in 6 of the 29 backward reactions, no apparent pattern observed

    alias_model_components()
  
    
#    # Grb2 binds SHC
#    Rule('shc_grb_binding', SHC(erb=ANY, grb=None, state='p') + Grb2(erb=None, shc=None,gab1=None) >>
#         SHC(erb=ANY, grb=1, state='p') % Grb2(erb=None, shc=1,gab1=None), kf_bind_shc_grb2)
#
    for place in comprtmnts:
        # Grb2 unbinds from SHC
        for erb in receptors[:2]:
            Rule(place+'_shc_grb_unbinding_'+erb.name, erb(cpp=None,gs=2,comp=place) % SHC(grb=1,erb=2, state='p') %
                 Grb2(erb=None,sos=None, shc=1, gab1=None) <> erb(cpp=None,gs=2,comp=place) % SHC(grb=None,erb=2, state='p') +
                 Grb2(erb=None,sos=None, shc=None,gab1=None), kr_bind_shc_grb2, kf_bind_shc_grb2)
                 
            Rule(place+'_shc_grbsos_unbinding_'+erb.name, erb(cpp=None,gs=3,comp=place) % SHC(grb=1, erb=3, state='p') %
                 Grb2(erb=None,sos=2, shc=1, gab1=None) % SOS(grb=2, ras=None, erk=None, state='up') <> erb(cpp=None,gs=3,comp=place) %
                 SHC(grb=None, erb=3, state='p') + Grb2(erb=None,sos=2, shc=None,gab1=None) % SOS(grb=2, erk=None, ras=None, state='up'), kr_bind_shc_grb2, kf_bind_shc_grb2)

########################################################
def RAS_binds_sos():
    
    "4 typical reactions:(1) dimer + RAS_GDP <-> dimer:RAS_GDP (2) dimer + RAS_GTP <-> dimer:RAS_GDP \
    (3) dimer + RAS_GDP <-> dimer:RAS_GTP (4) dimer + RAS_active_GTP <-> dimer:RAS_GTP "
    
    Parameter('kf_gdp_gdp',2.e-5)    # k18
    Parameter('kr_gdp_gdp', 1.3)        # kd18

    Parameter('kf_gtp_gdp',1.6e-7)    # k19
    Parameter('kr_gtp_gdp', 0.5)        # kd19

    Parameter('kf_gdp_gtp',3.6e-7)    # k21
    Parameter('kr_gdp_gtp', 0.23)       # kd21

    Parameter('kf_agtp_gtp', 1.1e-5)  # k20
    Parameter('kr_agtp_gtp', 0.4)       # kd20
    
    alias_model_components()

    # RAS binds SOS
    for place in comprtmnts:
        for erb in receptors[:2]:
            Rule(place+'_Ras_gdp_gdp_'+erb.name, erb(cpp=None,gap=ANY,comp=place) % GAP() % SOS(grb=ANY,erk=None, ras=None, state='up') +
                 RAS(sos=None, pi3k=None, raf=None, state='gdp') <>  erb(cpp=None, gap=ANY, comp=place) % GAP() % SOS(grb=ANY,erk=None, ras=1, state='up') %
                 RAS(sos=1, pi3k=None, raf=None, state='gdp'), kf_gdp_gdp, kr_gdp_gdp)

            Rule(place+'_Ras_gtp_gdp_'+erb.name, erb(cpp=None, gap=ANY, comp=place) % GAP() % SOS(grb=ANY,erk=None, ras=None, state='up') +
                 RAS(sos=None, pi3k=None, raf=None, state='gtp') <> erb(cpp=None, gap=ANY, comp=place) % GAP() % SOS(grb=ANY, erk=None, ras=1, state='up') %
                 RAS(sos=1, pi3k=None,raf=None,state='gdp'), kf_gtp_gdp, kr_gtp_gdp)

            Rule(place+'_Ras_gdp_gtp_'+erb.name, erb(cpp=None, gap=ANY, comp=place) % GAP() % SOS(grb=ANY,erk=None, ras=None, state='up') +
                 RAS(sos=None, pi3k=None, raf=None, state='gdp') <> erb(cpp=None, gap=ANY, comp=place) % GAP() % SOS(grb=ANY,erk=None, ras=1, state='up') %
                 RAS(sos=1, pi3k=None, raf=None, state='gtp'), kf_gdp_gtp, kr_gdp_gtp)

            Rule(place+'_Ras_agtp_gtp'+erb.name, erb(cpp=None, gap=ANY,comp=place) % GAP() % SOS(grb=ANY,erk=None, ras=None, state='up') +
                 RAS(sos=None, pi3k=None, raf=None,state='active_gtp') <> erb(cpp=None, gap=ANY, comp=place) % GAP() % SOS(grb=ANY,erk=None, ras=1, state='up') %
                 RAS(sos=1, pi3k=None, raf=None,state='gtp'), kf_agtp_gtp, kr_agtp_gtp)

########################################################
def RTK_phos():
    " v650-v663 "
    " Erb dimers that are phoshproylated, but not bound to GAP can bind RTK (ES), followed by release of Enzyme RTK and \
    dephoshphorylated dimer product "
    "Reaction occurs ony for receptor dimers in the endo compartment "
    " In the case of Erb1 type dimers, EGF returns in the unphosphorylated product"

    Parameter('kf_bind_erb_rtk', 5e-05)     # k94
    Parameter('kr_bind_erb_rtk', 0.01)      # kd94
    Parameter('kcat_erb_rtk', 33)           # kd95
    Parameter('Erb2_dimer_assoc', 8.3e-9)   # k103
    Parameter('Erb2_dimer_disocc', 0.016)   # kd103

    alias_model_components()

    # check if catalyze can be used for these reactions
    for erb in receptors:
        Rule('bind_RTK_erb1_'+erb.name, RTK(erb=None) + Erb1(d=1, state='p', rtk=None, gap=None, comp='endo')
             % erb(d=1, state='p', gap=None,rtk=None, comp='endo') <> RTK(erb=2)% Erb1(d=1, state='p', rtk=2, gap=None, comp='endo') %
             erb(d=1, gap=None, state='p',rtk=None, comp='endo'), kf_bind_erb_rtk, kr_bind_erb_rtk)
             
        Rule('kcat_RTK_erb1_'+erb.name, RTK(erb=2) % Erb1(lig=None, d=1, state='p', rtk=2, gap=None, comp='endo') %
             erb(d=1, gap=None, state='p',rtk=None, comp='endo') >> RTK(erb=None) + EGF(rec=3, comp='endo') %
             Erb1(lig=3,d=1, state='up', rtk=None, gap=None, comp='endo') % erb(d=1, state='up', gap=None,rtk=None, comp='endo'), kcat_erb_rtk)
             
    for erb in receptors[1:]:
        Rule('bind_RTK_erb2_'+erb.name, RTK(erb=None) + Erb2(d=1, state='p', rtk=None, gap=None, comp='endo') %
             erb(d=1, state='p',gap=None,rtk=None, comp='endo') <> RTK(erb=2)% Erb2(d=1, state='p', rtk=2, gap=None, comp='endo') %
             erb(d=1, gap=None, state='p',rtk=None, comp='endo'), kf_bind_erb_rtk, kr_bind_erb_rtk)
             
        Rule('kcat_RTK_erb2_'+erb.name, RTK(erb=2)% Erb2(d=1, state='p', rtk=2, gap=None, comp='endo') %
             erb(d=1, gap=None, state='p',rtk=None, comp='endo') >> RTK(erb=None) + Erb2(d=1, state='up', rtk=None, gap=None, comp='endo') %
                erb(d=1, state='up', gap=None,rtk=None, comp='endo'), kcat_erb_rtk)
#        
#        # disassociation of Erb2 dimers (v678, v679)
#        Rule('dissoc_Erb2_'+erb.name, Erb2(d=1, state='up', lig=None, rtk=None, gap=None, comp='endo') %
#             erb(d=1, state='up', lig=None, gap=None,rtk=None, comp='endo') >> Erb2(d=None,lig=None, state='up', rtk=None, gap=None, comp='endo') +
#             erb(d=None, lig=None, state='up', gap=None,rtk=None, comp='endo'), Erb2_dimer_disocc)

########################################################
def bind_cPP():
    " v43-v161 "
    " cPP binds any Grb2 bound dimer. Once bound by cPP, ErbB dimers cannot participate in any reactions \
        (for eg: release Grb, bind SOS, RAS) unless cPP is released "
    " binding of cPP(any compartment) and receptor dimer(any compartment) gives a bound product in the plasma compartment "

#Parameter('kf_bind_cpp_erb', 0)             # k4 = 6.73e-6; 4b = 5 =  5b = 0
    Parameter('kr_bind_cpp_erb', 0.0080833)     # kd4b = 1.66e-4; kd5b = 0.0080833

    Parameter('k4', 6.73e-6)                    # forward rate constant for Erb1 homodimers (plasma mem) binding cPP
    Parameter('k4b', 0)                         # forward rate constant for Erb heterodimers (plasma mem) binding cPP
    Parameter('k5', 0)                          # forward rate constant for Erb1 dimers and Erb2 homodimers (endo) binding cPP
    Parameter('k5b', 0)                         # forward rate constant for Erb2 heterodimers binding cPP
    
    Parameter('kd4', 1.66e-4)                   # generic reverse rate constant for cPP release into plasma compartment
    Parameter('kd5b', 0.0080833)                # generic reverse rate constant for cPP release into endo compartment
    
    alias_model_components()
    
#    
#    # EGF:Erb1 dimers require separate rules since bound ligand is also transported to pm compartment
#    # 2(EGF:Erb1):Grb2(){compartment = ANY) + cPP <-> 2(EGF:Erb1):Grb2(): CPP { compartment = plasma membrane}
#    Rule(place+'_bind_cPP_erb1_erb1', cPP(erb=None,comp='pm') + Grb2(gab1=None) % EGF(rec=3, comp='pm') %
#        Erb1(lig=3, d=1, cpp=None, gap=ANY, comp='pm') % EGF(rec=4, comp='pm') % Erb1(lig=4, d=1,gap=None, comp=place) <>
#        cPP(erb=2, comp='pm') % Grb2(gab1=None) % EGF(rec=3, comp='pm') % Erb1(lig=3, d=1,cpp=2, gap=ANY, comp='pm') %
#        EGF(rec=4, comp='pm') % Erb1(lig=4, d=1,gap=None, comp='pm'), kf_bind_cpp_erb, kr_bind_cpp_erb)
#    
#     # Erb1:Erb(2/3/4) : Grb2(){compartment = ANY) + cPP <-> Erb1:Erb2/3/4:Grb2(): CPP { compartment = plasma membrane}
#     for erb in receptors[1:]:
#         Rule(place+'_bind_cPP_erb1_'+erb.name, cPP(erb=None,comp=place) + Grb2(gab1=None) % Erb1(d=1, cpp=None, gap=ANY, comp=place) %
#              erb(d=1,gap=None, comp=place) <> cPP(erb=2, comp='pm') % Grb2(gab1=None) % Erb1(d=1,cpp=2, gap=ANY, comp='pm') %
#              erb(d=1,gap=None, comp='pm'), kf_bind_cpp_erb, kr_bind_cpp_erb)
#     
#     # Erb2:Erb(2/3/4) : Grb2(){compartment = ANY) + cPP <-> Erb2:Erb2/3/4:Grb2(): CPP { compartment = plasma membrane}
#     for erb in receptors[1:]:
#         Rule(place+'_bind_cPP_erb2_'+erb.name, cPP(erb=None, comp=place) + Grb2(gab1=None) % Erb2(d=1, cpp=None, gap=ANY, comp=place) %
#              erb(d=1,gap=None, comp=place) <> cPP(erb=2, comp='pm') % Grb2(gab1=None) % Erb2(d=1,cpp=2, gap=ANY, comp='pm') %
#              erb(d=1,gap=None, comp='pm'), kf_bind_cpp_erb, kr_bind_cpp_erb)
#for place in comprtmnts:
#        if place == 'pm':
#            kr_bind_cpp_erb.value = kd4.value
#        #k4.value = 6.73e-06
#        elif place == 'endo':
#            kr_bind_cpp_erb.value = kd5b.value
#            k4.value = k4b.value

    place = 'pm'
    # EGF:Erb1 dimers require separate rules since bound ligand is also transported to pm compartment
    # 2(EGF:Erb1):Grb2(){compartment = ANY) + cPP <-> 2(EGF:Erb1):Grb2(): CPP { compartment = plasma membrane}
    Rule(place+'_bind_cPP_erb1_erb1_noSoS', cPP(erb=None,comp=place) + Grb2(sos=None, gab1=None) % EGF(rec=3, comp=place) %
         Erb1(lig=3, d=1, cpp=None, gap=ANY, comp=place) % EGF(rec=4, comp=place) % Erb1(lig=4, d=1,gap=None, comp=place) <>
         cPP(erb=2, comp='pm') % Grb2(sos=None, gab1=None) % EGF(rec=3, comp='pm') % Erb1(lig=3, d=1,cpp=2, gap=ANY, comp='pm') %
         EGF(rec=4, comp='pm') % Erb1(lig=4, d=1,gap=None, comp='pm'), k4, kd4)
         
    Rule(place+'_bind_cPP_erb1_erb1', cPP(erb=None,comp=place) + SOS(erk=None, state='up') % Grb2(sos= ANY,gab1=None) % EGF(rec=3, comp=place) %
            Erb1(lig=3, d=1, cpp=None, gap=ANY, comp=place) % EGF(rec=4, comp=place) % Erb1(lig=4, d=1,gap=None, comp=place) <>
            cPP(erb=2, comp='pm') % SOS(erk=None, state='up') % Grb2(sos= ANY,gab1=None) % EGF(rec=3, comp='pm') % Erb1(lig=3, d=1,cpp=2, gap=ANY, comp='pm') %
            EGF(rec=4, comp='pm') % Erb1(lig=4, d=1,gap=None, comp='pm'), k4, kd4)
    
    # Erb1:Erb(2/3/4) : Grb2(){compartment = ANY) + cPP <-> Erb1:Erb2/3/4:Grb2(): CPP { compartment = plasma membrane}
    for erb in receptors[1:]:
        Rule(place+'_bind_cPP_erb1_'+erb.name, cPP(erb=None,comp=place) + Grb2(gab1=None) % Erb1(d=1, cpp=None, gap=ANY, comp=place) %
             erb(d=1,gap=None, comp=place) <> cPP(erb=2, comp='pm') % Grb2(gab1=None) % Erb1(d=1,cpp=2, gap=ANY, comp='pm') %
             erb(d=1,gap=None, comp='pm'), k5, kd4)
    
    # Erb2:Erb(2/3/4) : Grb2(){compartment = ANY) + cPP <-> Erb2:Erb2/3/4:Grb2(): CPP { compartment = plasma membrane}
    for erb in receptors[1:]:
        Rule(place+'_bind_cPP_erb2_'+erb.name, cPP(erb=None, comp=place) + Grb2(gab1=None) % Erb2(d=1, cpp=None, gap=ANY, comp=place) %
             erb(d=1,gap=None, comp=place) <> cPP(erb=2, comp='pm') % Grb2(gab1=None) % Erb2(d=1,cpp=2, gap=ANY, comp='pm') %
             erb(d=1,gap=None, comp='pm'), k5, kd4)

    place = 'endo'
    # EGF:Erb1 dimers require separate rules since bound ligand is also transported to pm compartment
    # 2(EGF:Erb1):Grb2(){compartment = ANY) + cPP <-> 2(EGF:Erb1):Grb2(): CPP { compartment = plasma membrane}
    Rule(place+'_bind_cPP_erb1_erb1_noSOS', cPP(erb=None,comp=place) + Grb2(gab1=None, sos=None) % EGF(rec=3, comp=place) %
         Erb1(lig=3, d=1, cpp=None, gap=ANY, comp=place) % EGF(rec=4, comp=place) % Erb1(lig=4, d=1,gap=None, comp=place) <>
         cPP(erb=2, comp='pm') % Grb2(gab1=None, sos=None) % EGF(rec=3, comp='pm') % Erb1(lig=3, d=1,cpp=2, gap=ANY, comp='pm') %
         EGF(rec=4, comp='pm') % Erb1(lig=4, d=1,gap=None, comp='pm'), k5, kd5b)

    Rule(place+'_bind_cPP_erb1_erb1', cPP(erb=None,comp=place) + SOS(erk=None, state='up')% Grb2(gab1=None, sos=ANY) % EGF(rec=3, comp=place) %
         Erb1(lig=3, d=1, cpp=None, gap=ANY, comp=place) % EGF(rec=4, comp=place) % Erb1(lig=4, d=1,gap=None, comp=place) <>
         cPP(erb=2, comp='pm') % SOS(erk=None, state='up') % Grb2(gab1=None, sos=ANY) % EGF(rec=3, comp='pm') % Erb1(lig=3, d=1,cpp=2, gap=ANY, comp='pm') %
         EGF(rec=4, comp='pm') % Erb1(lig=4, d=1,gap=None, comp='pm'), k5, kd5b)
    
         # Erb1:Erb(2/3/4) : Grb2(){compartment = ANY) + cPP <-> Erb1:Erb2/3/4:Grb2(): CPP { compartment = plasma membrane}
    for erb in receptors[1:]:
             Rule(place+'_bind_cPP_erb1_'+erb.name, cPP(erb=None,comp=place) + Grb2(gab1=None) % Erb1(d=1, cpp=None, gap=ANY, comp=place) %
                  erb(d=1,gap=None, comp=place) <> cPP(erb=2, comp='pm') % Grb2(gab1=None) % Erb1(d=1,cpp=2, gap=ANY, comp='pm') %
                  erb(d=1,gap=None, comp='pm'), k5, kd5b)

    # Erb2:Erb(2/3/4) : Grb2(){compartment = ANY) + cPP <-> Erb2:Erb2/3/4:Grb2(): CPP { compartment = plasma membrane}
    for erb in receptors[1:]:
        Rule(place+'_bind_cPP_erb2_'+erb.name, cPP(erb=None, comp=place) + Grb2(gab1=None) % Erb2(d=1, cpp=None, gap=ANY, comp=place) %
         erb(d=1,gap=None, comp=place) <> cPP(erb=2, comp='pm') % Grb2(gab1=None) % Erb2(d=1,cpp=2, gap=ANY, comp='pm') %
         erb(d=1,gap=None, comp='pm'), k5, kd5b)

#
################### ErbB_PI3K pathway ##############################
def bind_Gab1():
    " v688-694 "
    " Gab1 binds receptor dimers via Grb2 only when (1) dimers are in plasma mem (2) Grb2 (no SOS) is directly bound to receptor "

    Parameter('kf_bind_Gab1', 6.67e-05)         # k105
    Parameter('kr_bind_Gab1', 0.1)              # kd105
    Parameter('kf_bind_ATP_Gab1', 1.87e-8)      # k122
    Parameter('kr_bind_ATP_Gab1', 1.0)          # kd122
    Parameter('kcat_Gab1_p', 0.177828)          # kd123

    alias_model_components()
    
    for erb in receptors[:2]:
        Rule('Gab1_binds_Grb2_'+erb.name, Gab1(grb2=None,atp=None,shp2=None,erk=None, state='up') +
             Grb2(sos=None, erb=2,gab1=None) % erb(gs=2, comp='pm',cpp=None) <>
             Gab1(grb2=1,atp=None,state='up', erk=None, shp2=None) % Grb2(sos=None, erb=2, gab1=1) % erb(gs=2, comp='pm',cpp=None),
             kf_bind_Gab1, kr_bind_Gab1)

    ## v30-v36, v815-v821
    ## Gab1 + ATP <-> Gab1:ATP -> Gab1~P + ATP
    catalyze(ATP(erb=None), 'gab1', Gab1(state='up',grb2=ANY, shp2=None,erk=None), 'atp', Gab1(state='p',grb2=ANY, shp2=None,erk=None),
             (kf_bind_ATP_Gab1, kr_bind_ATP_Gab1, kcat_Gab1_p))

########################################################
def Shp2_catalysis():
    "v707-v720"
    " Shp2 binds phoshphoryated Gab1 and dephosphorylates it in a  2-step reaction "

    Parameter('kf_bind_Shp2', 3.3e-5)   # k107
    Parameter('kr_bind_Shp2',0.1)       # kd107
    Parameter('kcat_Shp2', 5)           # kd108

    alias_model_components()

    catalyze_state(Shp2(), 'gab1', Gab1(erk=None, pi3k=None), 'shp2', 'state', 'p', 'up', (kf_bind_Shp2, kr_bind_Shp2, kcat_Shp2))

########################################################
def Erk_catalysis():
    "v723-v750"
    " Erk~PP binds and phosphorylates Gab1~P via 2-step catalysis"

    Parameter('kf_bind_Erk', 3.3e-4)    # k110
    Parameter('kr_bind_Erk', 0.1)       # kd110
    Parameter('kcat_Erk', 6.57)         # kd111

    alias_model_components()

    catalyze_state(ERK(pase3=None, mek=None, sos=None, state='pp'), 'gab1', Gab1(shp2=None, pi3k=None), 'erk', 'state', 'p', 'pp', (kf_bind_Erk, kr_bind_Erk, kcat_Erk))

########################################################
def Pase_9t_catalysis():
    "v770-v783"
    "Pase_9t binds and dephoshhorylates Gab1~PP to Gab1~P in a 2 step catalysis"

    Parameter('kf_bind_Pase9t', 8.3e-8)     # k117
    Parameter('kr_bind_Pase9t', 0.1)        # kd117
    Parameter('kcat_Pase9t', 0.03)          # k118

    alias_model_components()
              
    catalyze_state(Pase_9t(),'gab1', Gab1(), 'pase', 'state', 'pp', 'p', (kf_bind_Pase9t, kr_bind_Pase9t, kcat_Pase9t))

########################################################
def bind_PI3K():
    " v621-v627 "
    " PI3K binds scaffolding protein Gab1~P "

    Parameter('kf_bind_pi3k', 1.5e-5)   # k66
    Parameter('kr_bind_pi3k', 0.2)      # kd66

    alias_model_components()

    Rule('bind_gab1_pI3k', Gab1(state='p', shp2=None, erk=None, pi3k=None) + PI3K(gab1=None, ras=None, pip2=None) <>
         Gab1(state='p', shp2=None, erk=None, pi3k=1) % PI3K(gab1=1, ras=None, pip2=None), kf_bind_pi3k, kr_bind_pi3k)

########################################################
def PIP2_PIP3():
    " v695-v701 "
    " v628-v638 "
    " Pip2 binds PI3K "
    " ***  Note Erb2:Erb4 reactions (v701) needs to be checked in the model "
    " **** Erb2:Erb3:PI3K can bind PIP2 repeatedly until 6 Pip2 are bound. Not included in the model yet "

    Parameter('kf_bind_Pip2', 1.33e-5)  # k106
    Parameter('kr_bind_Pip2', 0.1)      # kd106
    Parameter('kcat_Pip2', 0.2)         # kd68

    alias_model_components()

    catalyze(PI3K(gab1=ANY, ras=None), 'pip2', PIP2(), 'pi3k', PIP3(akt=None, pdk=None, bnd=None), (kf_bind_Pip2, kr_bind_Pip2, kcat_Pip2))

    pip2_number = 6;
        #   def pip2_pip3_binding(pip2_number):
    for i in range(1, 1 + pip2_number -1):
        #print ""
        num_pip = str(i)
        pip2_common = PIP2(pi3k=1)
        for d in range(2, i+1):
            pip2_common %= PIP2(pi3k=d)
        reactant = Erb2() % Erb3() % PI3K(ras=None, gab1=999, pip2=range(1,i+1)) % pip2_common
        product = Erb2() % Erb3() %PI3K(ras=None, gab1=999, pip2=range(1, i+2)) % pip2_common % PIP2(pi3k=i+1)
        #print reactant
        #print product
        Rule('catalye_pip2_pip3'+num_pip, reactant + PIP2(pi3k=None) >> product, kf_bind_Pip2)

########################################################
def PI3K_binds_RAS():
    "v751-764 , check v764 "
    " Ras GDP/GTP + dimer:PI3K <> PI3K:Ras_GDP "

    Parameter('kf_ras_gdp_gdp',0.0047)    # k112
    Parameter('kr_ras_gdp_gdp', 0.1) # kd112
    Parameter('kf_ras_gtp_gdp', 0)    # k113 *****
    Parameter('kr_ras_gtp_gdp', 177) # kd113

    alias_model_components()

    Rule('pi3k_gdp_gdp', PI3K(gab1=ANY, pip2=None, ras=None) + RAS(sos=None,raf=None, pi3k=None, state='gdp') <>
         PI3K(gab1=ANY, ras=1, pip2=None) % RAS(sos=None,pi3k=1, raf=None, state='gdp'), kf_ras_gdp_gdp, kr_ras_gdp_gdp)

    Rule('pi3k_gtp_gdp', PI3K(gab1=ANY, pip2=None, ras=None) + RAS(sos=None,pi3k=None, raf=None, state='gtp') <>
         PI3K(gab1=ANY, ras=1, pip2=None) % RAS(sos=None,pi3k=1, raf=None, state='gdp'), kf_ras_gtp_gdp, kr_ras_gtp_gdp)

########################################################
def AKT_rxns():
    "v639-v649, v765-768"
    
    Parameter('kf_pip3_akt',3.3e-5)   # k69
    Parameter('kr_pip3_akt',0.1)   # kd69
    
    Parameter('kf_pip3_pdk',6.6e-7)   # k70
    Parameter('kr_pip3_pdk',0.1)   # kd70
    
    Parameter('kd71', 25)
    
    Parameter('kd72', 5)
    
    Parameter('kdissoc_pip3_pdk1', 142)   # kd76
    
    Parameter('kf_aktpp_Pase4',6.36e-7)    #k74
    Parameter('kr_aktpp_Pase4',0.355)    # kd74
    Parameter('kf_aktp_Pase4',0.00374845)     # k73
    Parameter('kr_aktp_Pase4',0.5)     # kd73
    Parameter('kcat_Pase4',0.00633957)        # kd75
    Parameter('kf_pip3_bnd', 5e-6)     # k109
    Parameter('kr_pip3_bnd',0.1)    # kd109
    Parameter('kcat_pip3_pip2', 0.2)      # kd104
    Parameter('kf_akt_raf', 4.98e-6)        # k114
    Parameter('kr_akt_raf',0.1)             # kd114
    Parameter('kcat_akt_raf', 1.0)      # kd115

    alias_model_components()
    
    pstates = ['up', 'p', 'pp']
    
    for ps in pstates[:2]:
        # PIP3 + AKT (up/p) <-> PIP3:AKT (up/p)
        Rule('Pip3_binds_AKT_'+ps, PIP3(akt=None, pdk=None, bnd=None) + AKT(pip=None, pase=None,raf=None, state=ps) <>
             PIP3(akt=1, pdk=None, bnd=None) % AKT(pip=1, pase=None,raf=None, state=ps), kf_pip3_akt, kr_pip3_akt)
        
        # PIp3:AKT + PDK1
        Rule('Pip3AKT_bind_PDK_'+ps, PIP3(akt=1, pdk=None, bnd=None) % AKT(pip=1, pase=None,raf=None, state=ps) +PDK1(pip=None) <>
             PIP3(akt=1, pdk=2, bnd=None) % AKT(pip=1, pase=None,raf=None, state=ps) % PDK1(pip=2), kf_pip3_pdk, kr_pip3_pdk)

    # PIp3:AKT:PDK1 -> PIp3:PDK1 + AKT~P
    Rule('AKT_first_phosphoryl', PIP3(akt=1, pdk=2, bnd=None) % AKT(pip=1, pase=None,raf=None, state='up') % PDK1(pip=2) >>
         PIP3(akt=None, pdk=2, bnd=None) % PDK1(pip=2) + AKT(pip=None,raf=None, pase=None, state='p'), kd71)

    # PIp3:AKT~P:PDK1 -> PIp3:PDK1 + AKT~PP
    Rule('AKT_second_phosphoryl', PIP3(akt=1, pdk=2, bnd=None) % AKT(pip=1,raf=None, pase=None, state='p') % PDK1(pip=2) >>
     PIP3(akt=None, pdk=2, bnd=None) % PDK1(pip=2) +  AKT(pip=None,raf=None, pase=None, state='pp') , kd72)

    # Dissoc pip3:PDk1
    Rule('dissoc_pip3_pdk1', PIP3(akt=None, pdk=1, bnd=None) % PDK1(pip=1) >> PIP3(akt=None, pdk=None, bnd=None) + PDK1(pip=None), kdissoc_pip3_pdk1)

    # Catalyse dephoshporylation of ATK by Pase4
    catalyze_state(Pase4(),'akt', AKT(pip=None,raf=None), 'pase', 'state', 'pp', 'p', (kf_aktpp_Pase4, kr_aktpp_Pase4, kcat_Pase4))

    catalyze_state(Pase4(),'akt', AKT(pip=None,raf=None), 'pase', 'state', 'p', 'up', (kf_aktp_Pase4, kr_aktp_Pase4, kcat_Pase4))

    # Pip3 binds SHP/PTEN to give Pip2
    for binder in Shp, PTEN:
        Rule('bind_Pip3_'+binder.name, binder(pip=None) + PIP3(akt=None, pdk=None, bnd=None) <> binder(pip=1) % PIP3(akt=None, pdk=None, bnd=1), kf_pip3_bnd, kr_pip3_bnd)
        Rule('Pip3_2_'+binder.name, binder(pip=1) % PIP3(akt=None, pdk=None,bnd=1) >> binder(pip=None) + PIP2(pi3k=None), kcat_pip3_pip2)

    # Akt convers Raf#p to Raf#p:Ser
    catalyze_state(AKT(pip=None, pase=None,state='pp'),'raf', RAF(pase1=None, mek=None), 'akt', 'state', 'p', 'p_ser', (kf_akt_raf, kr_akt_raf, kcat_akt_raf))


################### ErbB_MAPK pathway ##############################
def MAPK_pathway():
    " v409-412, v487-512 "
    Parameter('kf_raf_ras', 5e-06 ) #k28
    Parameter('kr_raf_ras', 0.0053) #kd28
    Parameter('kcat_raf_ras', 3.1) # kd29
    Parameter('krev_raf_ras', 1.17e-06) # k29

    Parameter('kf_raf_pase1', 6e-5)    # k42
    Parameter('kr_raf_pase1', 0.014)    # kd42
    Parameter('kcat_raf_pase1', 31.6228)  # kd43
    
    Parameter('kf_raf_mek', 1.07e-5)    # k44
    Parameter('kr_raf_mek', 0.033)    # kd52
    Parameter('kcat_raf_mek', 1.9)  # kd45
    
    Parameter('kf_mek_pase2', 2.37e-5)    # k48
    Parameter('kr_mek_pase2', 0.79)    # kd48
    Parameter('kcat_mek_pase2', 0.11)  # kd49
    
    Parameter('kf_erk_pase3', 3.97e-4)    # k56
    Parameter('kr_erk_pase3', 5.0)    # kd56
    Parameter('kcat_erk_pase3', 0.0076)  # kd57
    
    Parameter('k58', 8.33e-7)    # k56
    Parameter('kd58', 56.7)    # kd56
    
    Parameter('kf_erk_sos',1)
    Parameter('kr_erk_sos', 0.1)
    Parameter('kcat_erk_sos', 1)
    
    alias_model_components()

    # RAS:GTP + RAF -> RasGTP:RAF
    # RAS:active_GTP + RAF~P -> RasGTP:RAF
    Rule('RAS_binds_RAF', RAS(sos=None, pi3k=None, raf=None, state='gtp') +
         RAF(akt=None, mek=None, pase1=None, ras=None,state='up') <>
         RAS(sos=None, pi3k=None, raf=1, state='gtp') %
         RAF(akt=None, mek=None, pase1=None, ras=1,state='up'), kf_raf_ras, kr_raf_ras)
         
    Rule('RAS_RAF_cat', RAS(sos=None, pi3k=None, raf=1, state='gtp')
         % RAF(akt=None, mek=None, pase1=None, ras=1,state='up') <>
         RAS(sos=None, pi3k=None, raf=None, state='active_gtp') +
         RAF(akt=None, mek=None, pase1=None, ras=None,state='p'),kcat_raf_ras, krev_raf_ras)

    # RAF1~P dephoshporylated by Pase1
    catalyze_state(Pase1(), 'raf', RAF(ras=None, akt=None, mek=None), 'pase1', 'state', 'p', 'up', (kf_raf_pase1, kr_raf_pase1, kcat_raf_pase1))

    # RAF~P phoshprylates MEK
    catalyze_state(RAF(akt=None, ras=None, pase1=None, state='p'), 'mek', MEK(pase2=None, erk=None), 'raf', 'state', 'up', 'p', (kf_raf_mek, kr_raf_mek, kcat_raf_mek))

    catalyze_state(RAF(akt=None, ras=None, pase1=None, state='p'), 'mek', MEK(pase2=None,erk=None), 'raf', 'state', 'p', 'pp', (kf_raf_mek, kr_raf_mek, kcat_raf_mek))
    
    
    # MEK~P dephosphorylated by Pase2
    catalyze_state(Pase2(), 'mek', MEK(raf=None,erk=None), 'pase2', 'state', 'p', 'up', (kf_mek_pase2, kr_mek_pase2, kcat_mek_pase2))
    
    catalyze_state(Pase2(), 'mek', MEK(raf=None,erk=None), 'pase2', 'state', 'pp', 'p', (kf_mek_pase2, kr_mek_pase2, kcat_mek_pase2))

    # MEK~P phosphorylates ERk
    catalyze_state(MEK(raf=None, pase2=None, state='pp'), 'erk', ERK(pase3=None, sos=None, gab1=None), 'mek', 'state', 'up', 'p', (kf_raf_mek, kr_raf_mek, kcat_raf_mek))

    catalyze_state(MEK(raf=None, pase2=None, state='pp'), 'erk', ERK(pase3=None, sos=None, gab1=None), 'mek', 'state', 'p', 'pp', (kf_raf_mek, kr_raf_mek, kcat_raf_mek))
    

    # ERK~P dephosphorylated by Pase3
    catalyze_state(Pase3(), 'erk', ERK(mek=None, sos=None, gab1=None), 'pase3', 'state', 'p', 'up', (k58, kd58, kcat_erk_pase3))
    
    catalyze_state(Pase3(), 'erk', ERK(mek=None,sos=None, gab1=None), 'pase3', 'state', 'pp', 'p', (kf_erk_pase3, kr_erk_pase3, kcat_erk_pase3))

    # ERK#PP phosphorylates SOS in ErbB1 homodimers
    Rule('sos_binds_erk', ERK(mek=None, gab1=None, pase3=None, state='pp', sos=None) +  Erb1(cpp=None, rtk=None) % Erb1(cpp=None, rtk=None) % SOS(ras=None, erk=None, state='up') <>
        ERK(mek=None, gab1=None, pase3=None, state='pp', sos=1) %  Erb1(cpp=None, rtk=None) % Erb1(cpp=None,rtk=None) % SOS(ras=None, erk=1, state='up'), kf_erk_sos, kr_erk_sos)
         
    # ERK#PP phosphorylaes free SOS
    Rule('freeSos_binds_erk', ERK(mek=None, gab1=None, pase3=None, state='pp', sos=None) + SOS(ras=None, erk=None, grb=None, state='up') <>
        ERK(mek=None, gab1=None, pase3=None, state='pp', sos=1) % SOS(ras=None, erk=1, grb=None, state='up'), kf_erk_sos, kr_erk_sos)

    Rule('Sos_p_catalysis', ERK(sos=1, state='pp') % SOS(erk=1,state='up') >> ERK(sos=None, state='pp') + SOS( erk=None,state = 'p'), kcat_erk_sos)


##########################################
def R_deg():
    Parameter('kdeg', 1)    # k60, k60b, k60c
    Parameter('k116', 0.015)
    alias_model_components()
    
    " The rules below tanslate to 92 reaction, instead of the 80 in the model. The remaining 12 appear to be missing degradation rxns in clusters 34, 35, 36, and 38 {Refer JHM diagram "
    
    for erb in receptors:
        # degrade GAP bound receptor dimers in endo comaprtment
        Rule('degrade_gap_bound_' + erb.name, erb(d=ANY, gap=ANY, comp='endo') >> None, kdeg)

        # degrade receptor dimers in endo that are not bound to ATP, rtk or phosphorylated
        Rule('degrade_noATP_bound_' + erb.name, erb(d=ANY, gap=None, comp='endo', atp=None, state='up', rtk=None)
             >> None, kdeg)

        # degrade ErbB monomers that are not phosphoryated or bound to ligand
        Rule('degrade_monomers_' +erb.name, erb(d=None, lig=None, comp='endo', state='up') >> None, kdeg)

    # Degrade Pase3 (V769)
    Rule('degrade_Pase3', Pase3(erk=None) >> None, k116)

##########################################
def R_deg_v2():
    Parameter('kdeg', 1)
    Parameter('k60', 0.0026)
    Parameter('k60b', 0.047)
    Parameter('k60c', 5.2e-4)
    Parameter('k62b', 4.1e-4)
    Parameter('k116', 0.015)
    alias_model_components()
    
    " The rules below tanslate to 92 reaction, instead of the 80 in the model. The remaining 12 appear to be missing degradation rxns in clusters 34, 35, 36, and 38 {Refer JHM diagram "
    
    # degrade GAP bound receptor dimers in endo comaprtment
    Rule('degrade_gap_bound_Erb1_homodimers', Erb1(d=1, gap=ANY, comp='endo') % Erb1(d=1, gap=None, comp='endo') % RAS() >> None, k60)
    
    for erb in receptors[1:]:
        Rule('degrade_gap_bound_Erb1_heterodimers_'+erb.name, Erb1(d=1, gap=ANY, comp='endo') % erb(d=1, gap=None, comp='endo') % RAS()>> None, k60b)
    
    Rule('degrade_gap_bound_Erb2_homodimers', Erb2(d=1, gap=ANY, comp='endo') % Erb2(d=1, gap=None, comp='endo') % RAS()>> None, k60b)
    
    for erb in receptors[2:]:
        Rule('degrade_gap_bound_Erb2_heterodimers_'+erb.name, Erb2(d=1, gap=ANY, comp='endo') % erb(d=1, gap=None, comp='endo')% RAS() >> None, k60c)
        
    # degrade receptor dimers in endo that are not bound to ATP, rtk or phosphorylated
    Rule('degrade_noATP_bound_Erb1_homodimers', Erb1(d=1, gap=None, comp='endo', atp=None, state='up', rtk=None)% Erb1(d=1, gap=None, comp='endo', atp=None, state='up', rtk=None)
             >> None, k62b)

    for erb in receptors[1:]:
        Rule('degrade_noATP_bound_Erb1_'+ erb.name, Erb1(d=1, gap=None, comp='endo', atp=None, state='up', rtk=None) %
                                                       erb(d=1, gap=None, comp='endo', atp=None, state='up', rtk=None) >> None, k62b)
             
    Rule('degrade_noATP_bound_Erb2_homodimers', Erb2(d=1, gap=None, comp='endo', atp=None, state='up', rtk=None) % Erb2(d=1, gap=None, comp='endo', atp=None, state='up', rtk=None)  >> None, k62b)
             
    for erb in receptors[1:]:
             Rule('degrade_noATP_bound_Erb2_'+ erb.name, Erb2(d=1, gap=None, comp='endo', atp=None, state='up', rtk=None) %
                                                            erb(d=1, gap=None, comp='endo', atp=None, state='up', rtk=None) >> None, k62b)
                  
         
    # degrade ErbB monomers that are not phosphoryated or bound to ligand
                  
    Rule('degrade_Erb1_monomer', Erb1(d=None, lig=None, atp=ANY, comp='endo', state='up') >> None, k60)
                  
    for erb in receptors[1:]:
                  Rule('degrade_monomers_' +erb.name, erb(d=None, lig=None, comp='endo', state='up') >> None, k60b)
    
    # Degrade Pase3 (V769)
    Rule('degrade_Pase3', Pase3(erk=None) >> None, k116)

    Rule('degrade_endo_EGF', EGF(comp='endo', rec=None) >> None, kdeg)

########################################################
def declare_observables():
    alias_model_components()
    
    Observable('pErbB1', Erb1(d=ANY,state='p'))
    Observable('pERK', ERK(state='pp'))
    Observable('pAKT', AKT(state='pp'))


if __name__ == '__main__':
    print __doc__
    print "NOTE: This model code is designed to be imported and " \
        "programatically manipulated,\nnot executed directly. The above " \
            " output is merely a diagnostic aid."


#Notes: Following parameters are not set as constants
#kd1 //
#k1c //
#kd1c //
#k1d
#kd1d //
#k6
#kd6 //
#k35
#kd35
#k48
#kd48
#k56
#kd56
#k58
#kd58
#kd73
#kd74
#kd75
#k74
#k75
#k97 //
#kd97
#k97c //
#kd97c
#k99
#k103
#kd103
#k116 //
#k119 //
#k120
#kd129
#k120b
#kd120b
#kd122 //
#k1
#k122 //