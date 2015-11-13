""" Recapitulation of the Chen 2009 ErbB signaling model.
    The model currently does not account for the following:
     - species ending in '_inh', '_i', 'half_active'
     - Events (referred to as rules in simbiology)
    When instantiated with all modules, the model has 436 species, 713 reactions and 212 parameters 
"""
from pysb import Model
from ..components import erbb, pi3k, mapk
from ..utils import merge_parameters

# Create the model.
Model()

# Declare monomers.
erbb.receptor_monomers()
erbb.ligand_monomers()
erbb.adapter_monomers()
pi3k.PI3K_monomers()
mapk.MAPK_monomers()

# Generate ErbB (upstream) pathway.
erbb.transport()
erbb.ErbB1_priming()
erbb.ligand_binding()
erbb.receptor_dimerization()
erbb.lateral_signaling()
erbb.secondary_dimerization()
erbb.trans_phosphorylation()
erbb.GAP_binding()
erbb.Grb2_binding()
erbb.SHC_binding()
erbb.secondary_Grb2_binding()
erbb.RAS_binds_sos()
erbb.RTK_phos()
erbb.bind_cPP()
erbb.R_deg_v2()
erbb.ErbB2_magic_rxns()

# Generate PI3K/AKT pathway.
pi3k.bind_Gab1()
pi3k.Shp2_catalysis()
pi3k.Erk_catalysis()
pi3k.Pase_9t_catalysis()
pi3k.bind_PI3K()
pi3k.PIP2_PIP3()
pi3k.PI3K_binds_RAS()
pi3k.AKT_rxns()

# Generate MAPK pathway.
mapk.MAPK_pathway()

# Define observables.
erbb.declare_observables()
pi3k.declare_observables()
mapk.declare_observables()

# Merge parameters that were redundantly specified in multiple component
# modules.
merge_parameters(model, 'k103', [k103_ls, k103_magic])
merge_parameters(model, 'kd103', [kd103_ls, kd103_magic])
merge_parameters(model, 'k122', [k122_gab, k122_ls, k122_priming])
merge_parameters(model, 'kd122', [kd122_gab, kd122_ls, kd122_priming])
merge_parameters(model, 'k16', [k16_scndry])
merge_parameters(model, 'kd123', [kd123_gab, kd123_ls])
merge_parameters(model, 'kd24', [kd24_scndry])
