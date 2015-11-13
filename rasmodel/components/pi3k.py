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
def PI3K_monomers():
    """ Declare monomers in the PI3K arm of the pathway, namely, 
        Gab1, ERK, PI3K, Shp2, PIP2, PIP3, PTEN, Shp, AKT, PDK1 Pase9t, and Pase3
        
        A description of sites on the monomers is given below
        ======================================================
        Gab1 sites: 'atp' is a site to bind ATP
                    'grb2' is a site to bind Grb2
                    'shp2' is a site to bind Shp2
                    'erk' is a site to bind ERK
                    'pase' is a site to bind Pase9t
                    'pi3k' is a site to bind PI3K
                    'state' denotes the phoshorylation status of the species, 
                    with 'up' denoting unphosphorylated,'p' donating singly-phoshorylated and 'pp' denoting doubly phoshporylated
       'gab1' are sites on Shp2, Pase9t, and PI3K to bind Gab1 
       'pip' are sites on AKT, Shp, PTEN, PDK1 to bind PIP3
    """

    Monomer('Gab1',['atp', 'grb2', 'shp2', 'state','erk','pase', 'pi3k'], {'state':['up','p', 'pp']})
    Monomer('Shp2',['gab1'])
    Monomer('Pase_9t', ['gab1'])
    Monomer('PI3K', ['gab1', 'pip2', 'ras'])  # 'pip2' and 'ras' are sites on PI3K to bind PIP2 and RAS
    Monomer('PIP2', ['pi3k']) # 'pi3k' is a site on PIP2 to bind PI3K
    Monomer('PIP3', ['akt', 'pdk','bnd']) # 'akt' and 'pdk' are sites on PIP3 to bind AKT and PDK respectively
    Monomer('AKT', ['pip', 'pase', 'raf', 'state'], {'state':['up','p', 'pp']})
    Monomer('PDK1',['pip'])
    Monomer('Shp', ['pip'])
    Monomer('PTEN', ['pip'])
    Monomer('Pase4', ['akt']) # 'akt' is a site on Pase4 to bind AKT

    alias_model_components()

    global receptors
    receptors = [ErbB1, ErbB2, ErbB3, ErbB4]

def bind_Gab1():
    " v688-694 "
    " Gab1 binds receptor dimers via Grb2 only when (1) dimers are in plasma mem (2) Grb2 (no SOS) is directly bound to receptor "
    
    # Initial amount
    # ==============
    Parameter('Gab1_0',94868.3)     # c426
    # Rate constant
    # ==============
    Parameter('k105', 6.67e-05)         # k105
    Parameter('kd105', 0.1)              # kd105
    Parameter('k122_gab', 1.8704e-8)      # k122
    Parameter('kd122_gab', 1.0)          # kd122
    Parameter('kd123_gab', 0.177828)          # kd123

    alias_model_components()
    
    # Initial conditions
    # ==============
    Initial(Gab1(atp=None, grb2=None, shp2=None, erk=None, pase=None, pi3k=None, state='up'), Gab1_0)
    
    # Rules
    # =====
    for erb in receptors[:2]:
        Rule('Gab1_binds_Grb2_'+erb.name, Gab1(grb2=None,atp=None,shp2=None,erk=None, state='up') +
             Grb2(sos=None, erb=2,gab1=None) % erb(gs=2, comp='pm',cpp=None) <>
             Gab1(grb2=1,atp=None,state='up', erk=None, shp2=None) % Grb2(sos=None, erb=2, gab1=1) % erb(gs=2, comp='pm',cpp=None),
             k105, kd105)

    ## v30-v36, v815-v821
    ## Gab1 + ATP <-> Gab1:ATP -> Gab1~P + ATP
    catalyze(ATP(erb=None), 'gab1', Gab1(state='up',grb2=ANY, shp2=None,erk=None), 'atp', Gab1(state='p',grb2=ANY, shp2=None,erk=None),
             (k122_gab, kd122_gab, kd123_gab))

########################################################
def Shp2_catalysis():
    "v707-v720"
    " Shp2 binds phoshphoryated Gab1 and dephosphorylates it in a  2-step reaction "
    
    # Initial amount
    # ==============
    Parameter('Shp2_0', 1e+6)       # c463
    # Rate constant
    # ==============
    Parameter('k107', 3.33e-5)   # k107
    Parameter('kd107',0.1)       # kd107
    Parameter('kd108', 5)           # kd108

    alias_model_components()
    
    # Initial conditions
    # ==============
    Initial(Shp2(gab1=None), Shp2_0)

    catalyze_state(Shp2(), 'gab1', Gab1(erk=None, pi3k=None), 'shp2', 'state', 'p', 'up', (k107, kd107, kd108))

########################################################
def Erk_catalysis():
    "v723-v750"
    " Erk~PP binds and phosphorylates Gab1~P via 2-step catalysis"
    
    Parameter('k110', 3.33e-4)    # k110
    Parameter('kd110', 0.1)       # kd110
    Parameter('kd111', 6.57)         # kd111

    alias_model_components()
    
    # Initial amount
    # ==============

    catalyze_state(ERK(pase3=None, mek=None, sos=None, state='pp'), 'gab1', Gab1(shp2=None, pi3k=None), 'erk', 'state', 'p', 'pp',
                   (k110, kd110, kd111))

########################################################
def Pase_9t_catalysis():
    "v770-v783"
    "Pase_9t binds and dephoshhorylates Gab1~PP to Gab1~P in a 2 step catalysis"
    
    # Initial amount
    # ==============
    Parameter('Pase9t_0', 0)        # c521 // zero ic
    # Rate constant
    # ==============
    Parameter('k117', 8.33e-8)     # k117
    Parameter('kd117', 0.1)        # kd117
    Parameter('kd118', 0.03)          # k118

    alias_model_components()
    
    # Initial conditions
    # ==============
    Initial(Pase_9t(gab1=None), Pase9t_0)
    
    # Rules
    # =====
    catalyze_state(Pase_9t(),'gab1', Gab1(), 'pase', 'state', 'pp', 'p', (k117, kd117, kd118))

########################################################
def bind_PI3K():
    " v621-v627 "
    " PI3K binds scaffolding protein Gab1~P "
    " yet to account for reactions with rc k67, kd67 "
    
    # Initial amount
    # ==============
    Parameter('PI3K_0', 3.55656e+7)    # c287
    # Rate constants
    # ==============
    Parameter('k66', 1.5e-5)   # k66  # K67 and kd67 also used for a subset of rxns, no pattern observed :(
    Parameter('kd66', 0.2)
    Parameter('k67', 5e-5)
    Parameter('kd67', 0.02)

    alias_model_components()
    
    # Initial conditions
    # ==============
    Initial(PI3K(gab1=None, pip2=None, ras=None), PI3K_0)

    # Rules
    # ======
    # This was almost a nice single rule, except a few of the species need a
    # different rate.
    for erb, other_erbs in (ErbB1, receptors), (ErbB2, receptors[1:]):
        for other_erb in other_erbs:
            rates = (k66, kd66)
            if ((erb is ErbB2 and (other_erb is ErbB2 or other_erb is ErbB3)) or
                (erb is ErbB1 and other_erb is ErbB3)):
                rates = (k67, kd67)
            Rule('_'.join((erb.name, other_erb.name, 'bind_gab1_pI3k')),
                 erb() % other_erb() % Gab1(state='p', shp2=None, erk=None, pi3k=None) + PI3K(gab1=None, ras=None, pip2=None) <>
                 erb() % other_erb() % Gab1(state='p', shp2=None, erk=None, pi3k=1) % PI3K(gab1=1, ras=None, pip2=None),
                 *rates)

########################################################
def PIP2_PIP3():
    " v695-v701 "
    " v628-v638 "
    " Pip2 binds PI3K "
    " ***  Note ErbB2:ErbB4 reactions (v701) needs to be checked in the model "
    " An issue with the model is that the ErbB2:Erbb4: Pi3k:Pip2 complex is mislabeled as PI3k(c455), resulting in the ES complex \
    never breaking down to give E + P"
    
    # Initial amount
    # ==============
    Parameter('PIP2_0', 393639)     # c444
    Parameter('PIP3_0', 0)          # c106 // zero ic
    # Rate constants
    # ==============
    Parameter('k106', 1.33e-5)
    Parameter('kd106', 0.1)
    Parameter('kd68', 0.2)
    Parameter('k106b', 2.63418e-8)
    Parameter('kd106b', 0.1)
    Parameter('kd68b', 20.5)
    
    alias_model_components()
    
    # Initial conditions
    # ==============
    Initial(PIP2(pi3k=None), PIP2_0)
    Initial(PIP3(akt=None, pdk=None, bnd=None), PIP3_0)
    
    # Rules
    # =====
    catalyze(ErbB1(gap=ANY) % PI3K(gab1=ANY, ras=None), 'pip2', PIP2(), 'pi3k', PIP3(akt=None, pdk=None, bnd=None), (k106b, kd106b, kd68))
    
    " Unlike ErbB1, could not write generic rule for ErbB2 dimers since ErbB2:ErbB4 ES complex is not broken down to E + P "
    for erb in receptors[1:]:
        Rule('bind_Pip2_Erb2_'+erb.name,ErbB2(gap=ANY) % erb(gap=None) % PI3K(gab1=ANY, ras=None, pip2=None) + PIP2(pi3k=None) <>
             ErbB2(gap=ANY) % erb(gap=None) % PI3K(gab1=ANY, ras=None, pip2=1) % PIP2(pi3k=1), k106, kd106)
    for erb in receptors[1:3]:
        forward_rate = kd68 if erb is ErbB2 else kd68b
        Rule('catalyze_pip3_erb2_'+erb.name, ErbB2(gap=ANY) % erb(gap=None) % PI3K(gab1=ANY, ras=None, pip2=1) % PIP2(pi3k=1) >>
             ErbB2(gap=ANY) % erb(gap=None) % PI3K(gab1=ANY, ras=None, pip2=None) + PIP3(akt=None, pdk=None, bnd=None), forward_rate)

    pip2_number = 6;
    for i in range(1, 1 + pip2_number -1):
        num_pip = str(i)
        pip2_common = PIP2(pi3k=1)
        for d in range(2, i+1):
            pip2_common %= PIP2(pi3k=d)
        reactant = ErbB2() % ErbB3() % PI3K(ras=None, gab1=999, pip2=range(1,i+1)) % pip2_common
        product = MatchOnce(ErbB2() % ErbB3() %PI3K(ras=None, gab1=999, pip2=range(1, i+2)) % pip2_common % PIP2(pi3k=i+1))

        Rule('bind_pip2_pi3k'+num_pip, reactant + PIP2(pi3k=None) <> product, k106, kd106)
        Rule('release_pip2'+num_pip, product  >> reactant + PIP3(akt=None, pdk=None, bnd=None), kd68b)

########################################################
def PI3K_binds_RAS():
    "v751-764 , check v764 "
    " Ras GDP/GTP + dimer:PI3K <> PI3K:Ras_GDP "
    
    # Initial amount
    # ==============

    Parameter('k112',0.0047067)    # k112
    Parameter('kd112', 0.1) # kd112
    Parameter('k113', 0)    # k113 *****
    Parameter('kd113', 177.828) # kd113

    alias_model_components()
    
    # Initial amount
    # ==============

    Rule('pi3k_gdp_gdp', PI3K(gab1=ANY, pip2=None, ras=None) + RAS(sos=None,raf=None, pi3k=None, state='gdp') <>
     PI3K(gab1=ANY, ras=1, pip2=None) % RAS(sos=None,pi3k=1, raf=None, state='gdp'), k112, kd112)
    
    " To account for the incorrect reaction v764, we need to specify receptor subtype in the rule below "
    for erb in receptors:
        Rule('ErbB1_'+erb.name+'_pi3k_gtp_gdp', ErbB1()%erb()% PI3K(gab1=ANY, pip2=None, ras=None) +
             RAS(sos=None,pi3k=None, raf=None, state='gtp') <>
             ErbB1()%erb()% PI3K(gab1=ANY, ras=1, pip2=None) % RAS(sos=None,pi3k=1, raf=None, state='gdp'), k113, kd113)
    
    for erb in receptors[1:3]:
        Rule('ErbB2_'+erb.name+'_pi3k_gtp_gdp', ErbB2()%erb()% PI3K(gab1=ANY, pip2=None, ras=None) +
             RAS(sos=None,pi3k=None, raf=None, state='gtp') <>
             ErbB2()%erb()% PI3K(gab1=ANY, ras=1, pip2=None) % RAS(sos=None,pi3k=1, raf=None, state='gdp'), k113, kd113)

    # Since Shpt does not have an additional bidning site for RAS, I use a hack of bind RAS with GRb2's empty SOS bidning site "
    Rule('magic_v764', ErbB4()%ErbB2()% Grb2(sos=None) % Gab1(state='p', pi3k=None, shp2=3) % Shp2(gab1=3) +
     RAS(sos=None,raf=None, pi3k=None, state='gtp') <>
     ErbB4()%ErbB2()% Grb2(sos=None)% Gab1(state='p',pi3k=2, shp2=None)% PI3K(gab1=2, ras=1, pip2=None) %
     RAS(sos=None,raf=None, pi3k=1, state='gdp'), k113, kd113)
#####################
def AKT_rxns():
    "v639-v649, v765-768"
    
    # Initial amount
    # ==============
    Parameter('AKT_0', 905000)      # c107
    Parameter('PDK1_0', 3.00416e8)        # c109
    Parameter('Pase4_0', 4.5e+5)   # c113
    Parameter('PTEN_0', 56100.9)   # c279
    Parameter('Shp_0', 2213.59)     # c461
    Parameter('RAF_0', 71131.2)     # c41
    # Rate constants
    # ==============
    Parameter('k69',3.33e-5)   # k69
    Parameter('kd69',0.1)   # kd69
    Parameter('k70',6.67e-7)   # k70
    Parameter('kd70',0.1)   # kd70
    Parameter('kd71', 25.2)
    Parameter('kd72', 5.01187)
    Parameter('kd76', 142.262)   # kd76
    Parameter('k74',6.36184e-7)    #k74
    Parameter('kd74',0.355656)    # kd74
    Parameter('k73',0.00374845)     # k73
    Parameter('kd73',0.5)     # kd73
    Parameter('kd75',0.00633957)        # kd75
    Parameter('k109', 5e-6)     # k109
    Parameter('kd109',0.1)    # kd109
    Parameter('kd104', 0.2)      # kd104
    Parameter('k114', 4.98816e-6)        # k114
    Parameter('kd114',0.1)             # kd114
    Parameter('kd115', 1.0)      # kd115
    
    alias_model_components()
    
    # Initial conditions
    # ==================
    Initial(AKT(pip=None, pase=None,raf=None, state='up'), AKT_0)
    Initial(PDK1(pip=None), PDK1_0)
    Initial(Pase4(akt=None), Pase4_0)
    Initial(PTEN(pip=None), PTEN_0)
    Initial(Shp(pip=None), Shp_0)
    Initial(RAF(akt=None, pase1=None, mek=None, ras=None, state='up'), RAF_0)
    
    # Rules
    # =====
    pstates = ['up', 'p', 'pp']
    
    for ps in pstates[:2]:
        # PIP3 + AKT (up/p) <-> PIP3:AKT (up/p)
        Rule('Pip3_binds_AKT_'+ps, PIP3(akt=None, pdk=None, bnd=None) + AKT(pip=None, pase=None,raf=None, state=ps) <>
             PIP3(akt=1, pdk=None, bnd=None) % AKT(pip=1, pase=None,raf=None, state=ps), k69, kd69)
            
             # PIp3:AKT + PDK1
        Rule('Pip3AKT_bind_PDK_'+ps, PIP3(akt=1, pdk=None, bnd=None) % AKT(pip=1, pase=None,raf=None, state=ps) +PDK1(pip=None) <>
                  PIP3(akt=1, pdk=2, bnd=None) % AKT(pip=1, pase=None,raf=None, state=ps) % PDK1(pip=2), k70, kd70)

    # PIp3:AKT:PDK1 -> PIp3:PDK1 + AKT~P
    Rule('AKT_first_phosphoryl', PIP3(akt=1, pdk=2, bnd=None) % AKT(pip=1, pase=None,raf=None, state='up') % PDK1(pip=2) >>
     PIP3(akt=None, pdk=2, bnd=None) % PDK1(pip=2) + AKT(pip=None,raf=None, pase=None, state='p'), kd71)
    
    # PIp3:AKT~P:PDK1 -> PIp3:PDK1 + AKT~PP
    Rule('AKT_second_phosphoryl', PIP3(akt=1, pdk=2, bnd=None) % AKT(pip=1,raf=None, pase=None, state='p') % PDK1(pip=2) >>
         PIP3(akt=None, pdk=2, bnd=None) % PDK1(pip=2) +  AKT(pip=None,raf=None, pase=None, state='pp') , kd72)
        
    # Dissoc pip3:PDk1
    Rule('dissoc_pip3_pdk1', PIP3(akt=None, pdk=1, bnd=None) % PDK1(pip=1) >> PIP3(akt=None, pdk=None, bnd=None) + PDK1(pip=None), kd76)
         
    # Catalyse dephoshporylation of ATK by Pase4
    catalyze_state(Pase4(),'akt', AKT(pip=None,raf=None), 'pase', 'state', 'pp', 'p', (k74, kd74, kd75))
         
    catalyze_state(Pase4(),'akt', AKT(pip=None,raf=None), 'pase', 'state', 'p', 'up', (k73, kd73, kd75))
         
    # Pip3 binds SHP/PTEN to give Pip2
    for binder in Shp, PTEN:
        Rule('bind_Pip3_'+binder.name, binder(pip=None) + PIP3(akt=None, pdk=None, bnd=None) <>
             binder(pip=1) % PIP3(akt=None, pdk=None, bnd=1), k109, kd109)
        Rule('Pip3_2_'+binder.name, binder(pip=1) % PIP3(akt=None, pdk=None,bnd=1) >> binder(pip=None) + PIP2(pi3k=None), kd104)

    # Akt convers Raf#p to Raf#p:Ser
    catalyze_state(AKT(pip=None, pase=None,state='pp'),'raf', RAF(pase1=None, mek=None), 'akt', 'state', 'p', 'p_ser',
                   (k114, kd114, kd115))

def declare_observables():
    alias_model_components()
    
    Observable('pAKT', AKT(state='pp'))
