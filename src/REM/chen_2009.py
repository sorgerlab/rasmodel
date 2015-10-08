""" The RAS Executable model (REM).
    The model currently does not account for the following:
     - species ending in '_inh', '_i', 'half_active'
     - Events (referred to as rules in simbiology)
    When instantiated with all modules, the model has 436 species, 713 reactions and 212 parameters 
"""
from pysb import *
from REM import ErbB_modules
from REM import PI3K_modules
from REM import MAPK_modules

Model()

# Declare monomers
ErbB_modules.receptor_monomers()
ErbB_modules.ligand_monomers()
ErbB_modules.adapter_monomers()
PI3K_modules.PI3K_monomers()
MAPK_modules.MAPK_monomers()

# Generate ErbB (upstream) pathway
ErbB_modules.transport()
ErbB_modules.ErbB1_priming()
ErbB_modules.ligand_binding()
ErbB_modules.receptor_dimerization()
ErbB_modules.lateral_signaling()
ErbB_modules.secondary_dimerization()
ErbB_modules.trans_phosphorylation()
ErbB_modules.GAP_binding()
ErbB_modules.Grb2_binding()
ErbB_modules.SHC_binding()
ErbB_modules.secondary_Grb2_binding()
ErbB_modules.RAS_binds_sos()
ErbB_modules.RTK_phos()
ErbB_modules.bind_cPP()
ErbB_modules.R_deg_v2()
ErbB_modules.ErbB2_magic_rxns()

# Generate PI3K/AKT pathway
PI3K_modules.bind_Gab1()
PI3K_modules.Shp2_catalysis()
PI3K_modules.Erk_catalysis()
PI3K_modules.Pase_9t_catalysis()
PI3K_modules.bind_PI3K()
PI3K_modules.PIP2_PIP3()
PI3K_modules.PI3K_binds_RAS()
PI3K_modules.AKT_rxns()

# Generate MAPK pathway
MAPK_modules.MAPK_pathway()

# Declare
ErbB_modules.declare_observables()
PI3K_modules.declare_observables()
MAPK_modules.declare_observables()



