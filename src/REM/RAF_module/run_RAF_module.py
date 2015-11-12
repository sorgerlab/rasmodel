# Run RAF module using different initial conditions
# for Ras and Vemurafenib levels
from pysb import *

Model()

import EGFR_to_RAS
EGFR_to_RAS.monomers()
EGFR_to_RAS.KRAS_activation()
EGFR_to_RAS.declare_observables()

import BRAF_module
BRAF_module.monomers()
BRAF_module.BRAF_dynamics()
BRAF_module.observables()

import ERK_phosphorylation
ERK_phosphorylation.monomers()
ERK_phosphorylation.by_BRAF_mut()
ERK_phosphorylation.DUSP_phospatase()
ERK_phosphorylation.declare_observables()

# from ERK_phosphorylation import by_BRAF_wt
import numpy as np
from pysb.integrate import Solver

Ras_range = 2 * np.linspace(0, 1000, num=5)
sos_range = 2 * np.linspace(0, 1000, num=10)
Vemurafenib_range = 2 * np.linspace(0, 1000, num=10)

RAF_WT_level = []
RAF_V600E_level = []
ERK_P_level = []
KRAS_level = []

for s in sos_range:
    wt_holder_row = []
    v600e_holder_row = []
    erk_holder_row = []
    kras_holder_row = []
    for v in Vemurafenib_range:

        model.parameters['SOS_0'].value = s
        model.parameters['Vem_0'].value = v

        ts = np.linspace(0, 100, 100)
        solver = Solver(model, ts)
        solver.run()

        wt_holder_row.append(solver.yobs['BRAF_WT_active'][-1])
        v600e_holder_row.append(solver.yobs['BRAF_V600E_active'][-1])
        erk_holder_row.append(solver.yobs['ERK_P'][-1])
        kras_holder_row.append(solver.yobs['active_KRAS'][-1])

    RAF_WT_level.append(wt_holder_row)
    RAF_V600E_level.append(v600e_holder_row)
    ERK_P_level.append(erk_holder_row)
    KRAS_level.append(kras_holder_row)

wt = np.array(RAF_WT_level, dtype='float')
v600e = np.array(RAF_V600E_level, dtype='float')
erkp = np.array(ERK_P_level, dtype='float')
kras = np.array(KRAS_level, dtype='float')

# Plots
# -----
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2)

# wt_c = ax1.pcolor(Vemurafenib_range, Ras_range, wt, vmin=0)
# mut_c = ax1.pcolor(Vemurafenib_range, sos_range, v600e, vmin=0)
kras_act = ax1.pcolor(Vemurafenib_range, sos_range, kras, vmin=0)
erk_p = ax2.pcolor(Vemurafenib_range, sos_range, erkp, vmin=0)

ax1.set_xlabel('Vemurafenib')
ax2.set_xlabel('Vemurafenib')
ax1.set_ylabel('SOS')
ax1.set_title('KRAS:GTP')
ax2.set_title('ERK_P')

plt.colorbar(kras_act, ax=ax1)
plt.colorbar(erk_p, ax=ax2)

plt.show()
