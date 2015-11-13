from pysb import *

Model()

import EGFR_to_RAS
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


from pysb.export import export
matlab_output = export(model, 'matlab')

with open('matlab_files/run_timecourse.m', 'w') as f:
    f.write(matlab_output)
