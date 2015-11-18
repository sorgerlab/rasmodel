from pysb import *

Model()

# import EGFR_to_RAS
import EGFR_to_RAS_dim as EGFR_to_RAS
EGFR_to_RAS.monomers()
EGFR_to_RAS.KRAS_activation()
EGFR_to_RAS.SOS_dephosphorylation()
EGFR_to_RAS.declare_observables()

import BRAF_module
BRAF_module.monomers()
BRAF_module.BRAF_dynamics()
BRAF_module.observables()

import ERK_phosphorylation
ERK_phosphorylation.monomers()
ERK_phosphorylation.by_BRAF_mut()
ERK_phosphorylation.DUSP_phospatase()
ERK_phosphorylation.ERK_feedback()
ERK_phosphorylation.declare_observables()

# import MEK_ERK_phosphorylation
# MEK_ERK_phosphorylation.monomers()
# MEK_ERK_phosphorylation.by_BRAF_mut()
# MEK_ERK_phosphorylation.MEK_phosphorylates_ERK()
# MEK_ERK_phosphorylation.DUSP_phospatase()
# MEK_ERK_phosphorylation.PP2A_phosphatase()
# MEK_ERK_phosphorylation.ERK_feedback()
# MEK_ERK_phosphorylation.declare_observables()

# import CRAF_module
# CRAF_module.monomers()
# CRAF_module.CRAF_binds_KRAS()


from pysb.export import export
matlab_output = export(model, 'matlab')

with open('matlab_files/run_timecourse.m', 'w') as f:
    f.write(matlab_output)
