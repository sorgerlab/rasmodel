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


def MAPK_monomers():
    """ Declare monomers in the MAPK arm of the pathway
        'state' denotes the phoshorylation status of the species, with 'up' denoting unphosphorylated,
        'p' denoting phoshorylated', 'pp' denoting doubly-phosphorylated
    """
    
    Monomer('Pase3', ['erk']) # 'erk' is a site on Pase3 to bind ERK
    Monomer('Pase1', ['raf']) # 'raf' is a site on Pase1 to bind RAF
    Monomer('Pase2', ['mek']) # 'mek' is a site on Pase2 to bind MEK
    Monomer('RAF', ['akt', 'pase1', 'mek','state', 'ras'], {'state':['up','p', 'p_ser']})
    Monomer('MEK', ['raf', 'pase2', 'erk', 'state'], {'state':['up','p', 'pp']})
    
    alias_model_components()

def MAPK_pathway():
    " v409-412, v487-512 "
    
    # Initial amount
    # ==============
    Parameter('MEK_0', 3020000)     # c47
    Parameter('ERK_0', 695000)      # c55
    Parameter('Pase1_0', 5e+4)      # c44
    Parameter('Pase2_0', 124480)    # c53
    Parameter('Pase3_0', 16870.2)   # c60
    # Rate constants
    # ==============
    Parameter('k28', 5e-06 ) #k28
    Parameter('kd28', 0.0053) #kd28
    Parameter('kd29', 3.1) # kd29
    Parameter('k29', 1.17e-06) # k29
    
    Parameter('k42', 6e-5)    # k42
    Parameter('kd42', 0.014)    # kd42
    Parameter('kd43', 31.6228)  # kd43
    
    Parameter('k44', 1.07e-5)    # k44
    Parameter('kd52', 0.033)    # kd52
    Parameter('kd45', 1.9)  # kd45
    Parameter('kd47', 0.8)  # kd45
    
    Parameter('k48', 2.37e-5)    # k48
    Parameter('kd48', 0.79)    # kd48
    Parameter('k50', 4.74e-8)    # k48
    Parameter('kd50', 0.2529)    # kd48
    Parameter('kd49', 0.11)  # kd49
    
    Parameter('k52', 2.37e-5)    # k48
    Parameter('kd44', 0.79)    # kd48
    Parameter('kd53', 0.11)  # kd49
    Parameter('kd55', 70.16)  # kd49
    
    
    Parameter('k56', 3.97e-4)    # k56
    Parameter('kd56', 5.0)    # kd56
    Parameter('kd57', 0.0076)  # kd57
    
    Parameter('k58', 8.33e-7)    # k56
    Parameter('kd58', 56.7)    # kd56
    
    Parameter('k64',1)
    Parameter('kd64', 0.1)
    Parameter('kd65', 1)
    
    alias_model_components()
    
    # Initial conditions
    # ==================
    Initial(MEK(raf=None, pase2=None, erk=None,state='up'), MEK_0)
    Initial(ERK(mek=None, pase3=None, sos=None, state='up', gab1=None), ERK_0)
    Initial(Pase1(raf=None), Pase1_0)
    Initial(Pase2(mek=None), Pase2_0)
    Initial(Pase3(erk=None), Pase3_0)
    
    # Rules
    # =====
    " v409-412, v487-512 "

    
    alias_model_components()
    
    # RAS:GTP + RAF -> RasGTP:RAF
    # RAS:active_GTP + RAF~P -> RasGTP:RAF
    Rule('RAS_binds_RAF', RAS(sos=None, pi3k=None, raf=None, state='gtp') +
         RAF(akt=None, mek=None, pase1=None, ras=None,state='up') <>
         RAS(sos=None, pi3k=None, raf=1, state='gtp') %
         RAF(akt=None, mek=None, pase1=None, ras=1,state='up'), k28, kd28)
        
    Rule('RAS_RAF_cat', RAS(sos=None, pi3k=None, raf=1, state='gtp')
              % RAF(akt=None, mek=None, pase1=None, ras=1,state='up') <>
              RAS(sos=None, pi3k=None, raf=None, state='active_gtp') +
              RAF(akt=None, mek=None, pase1=None, ras=None,state='p'),kd29, k29)
         
     # RAF1~P dephoshporylated by Pase1
    catalyze_state(Pase1(), 'raf', RAF(ras=None, akt=None, mek=None), 'pase1', 'state', 'p', 'up', (k42, kd42, kd43))
     
     # RAF~P phoshprylates MEK
    catalyze_state(RAF(akt=None, ras=None, pase1=None, state='p'), 'mek', MEK(pase2=None, erk=None), 'raf', 'state', 'up', 'p', (k44, kd52, kd45))
     
    catalyze_state(RAF(akt=None, ras=None, pase1=None, state='p'), 'mek', MEK(pase2=None,erk=None), 'raf', 'state', 'p', 'pp', (k44, kd52, kd47))
     
     
     # MEK~P dephosphorylated by Pase2
    catalyze_state(Pase2(), 'mek', MEK(raf=None,erk=None), 'pase2', 'state', 'p', 'up', (k50, kd50, kd49))
     
    catalyze_state(Pase2(), 'mek', MEK(raf=None,erk=None), 'pase2', 'state', 'pp', 'p', (k48, kd48, kd49))
     
     # MEK~P phosphorylates ERk
    catalyze_state(MEK(raf=None, pase2=None, state='pp'), 'erk', ERK(pase3=None, sos=None, gab1=None), 'mek', 'state', 'up', 'p', (k52, kd44, kd53))
     
    catalyze_state(MEK(raf=None, pase2=None, state='pp'), 'erk', ERK(pase3=None, sos=None, gab1=None), 'mek', 'state', 'p', 'pp', (k52, k44, kd55))
     
     
     # ERK~P dephosphorylated by Pase3
    catalyze_state(Pase3(), 'erk', ERK(mek=None, sos=None, gab1=None), 'pase3', 'state', 'p', 'up', (k58, kd58, kd57))
     
    catalyze_state(Pase3(), 'erk', ERK(mek=None,sos=None, gab1=None), 'pase3', 'state', 'pp', 'p', (k56, kd56, kd57))
     
    # ERK#PP phosphorylates SOS in ErbB1 homodimers "only
    Rule('sos_binds_erk', ERK(mek=None, gab1=None, pase3=None, state='pp', sos=None) +  ErbB1(cpp=None, rtk=None, comp='pm') % ErbB1(cpp=None, rtk=None, comp='pm') % SOS(ras=None, erk=None, state='up') <>
          ERK(mek=None, gab1=None, pase3=None, state='pp', sos=1) %  ErbB1(cpp=None, rtk=None, comp='pm') % ErbB1(cpp=None,rtk=None, comp='pm') % SOS(ras=None, erk=1, state='up'), k64, kd64)
     
     # ERK#PP phosphorylaes free SOS
    Rule('freeSos_binds_erk', ERK(mek=None, gab1=None, pase3=None, state='pp', sos=None) + SOS(ras=None, erk=None, grb=None, state='up') <>
          ERK(mek=None, gab1=None, pase3=None, state='pp', sos=1) % SOS(ras=None, erk=1, grb=None, state='up'), k64, kd64)
     
    Rule('free_Sos_p_catalysis', ERK(sos=1, state='pp') % SOS(erk=1,state='up',grb=None) >> ERK(sos=None, state='pp') + SOS(grb=None, erk=None,state='p'), kd65)
     
    Rule('Sos_p_catalysis', ERK(sos=1, state='pp') % SOS(erk=1,state='up') % ErbB1(cpp=None, rtk=None, comp='pm') % ErbB1(cpp=None, rtk=None, comp='pm')   >> ERK(sos=None, state='pp') + SOS(erk=None,state='p') % ErbB1(cpp=None, rtk=None, comp='pm') % ErbB1(cpp=None, rtk=None, comp='pm'), kd65)

    # Degrade Pase3 (V769)
    Rule('degrade_Pase3', Pase3(erk=None) >> None, k116)

def declare_observables():
    alias_model_components()
    
    Observable('pERK', ERK(state='pp'))