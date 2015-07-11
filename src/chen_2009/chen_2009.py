from pysb import *

Model()
#from ErbB_rules_post_gap import*
#from rules_with_compartments import*
from rules_PI3K_MAPK import *

declare_monomers()

declare_initial_conditions()
transport()
ErbB1_priming()
ligand_binding()
receptor_dimerization()
#lateral_signaling()
secondary_dimerization()
trans_phosphorylation()
GAP_binding()
Grb2_binding_v2()
SHC_binding_v2()
secondary_Grb2_binding()
RAS_binds_sos()
RTK_phos()
bind_cPP()
bind_Gab1()
Shp2_catalysis()
Erk_catalysis()
Pase_9t_catalysis()
bind_PI3K()
PIP2_PIP3()
PI3K_binds_RAS()
AKT_rxns()
MAPK_pathway()
R_deg_v2()
declare_observables()


if __name__ == '__main__':
        print __doc__
        print "NOTE: This model code is designed to be imported and " \
                "programatically manipulated,\not executed directly. The above " \
                " output is merely a diagnostic aid."
