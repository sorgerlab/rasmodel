from rasmodel.scenarios.default import model

import numpy as np
from matplotlib import pyplot as plt
from pysb.integrate import Solver
from pysb import *

from tbidbaxlipo.util import fitting


KRAS = model.monomers['KRAS']
GTP = model.monomers['GTP']

total_pi = 50000

for mutant in KRAS.site_states['mutant']:
    Initial(KRAS(gtp=1, gap=None, gef=None, p_loop=None, s1s2='open', CAAX=None,
             mutant=mutant) % GTP(p=1, label='n'),
             Parameter('KRAS_%s_GTP_0' % mutant, 0))

plt.ion()
plt.figure()

t = np.linspace(0, 1000, 1000) # 1000 seconds

for mutant in KRAS.site_states['mutant']:
    # Zero out all initial conditions
    for ic in model.initial_conditions:
        ic[1].value = 0
    model.parameters['KRAS_%s_GTP_0' % mutant].value = total_pi

    sol = Solver(model, t)
    sol.run()
    plt.plot(t, sol.yobs['Pi_'] / total_pi, label=mutant)
    plt.ylabel('GTP hydrolyzed (%)')
    plt.ylim(top=1)
    plt.xlabel('Time (s)')
    plt.title('Intrinsic hydrolysis')

plt.legend(loc='upper left', fontsize=11, frameon=False)

plt.figure()
for mutant in KRAS.site_states['mutant']:
    # Zero out all initial conditions
    for ic in model.initial_conditions:
        ic[1].value = 0
    model.parameters['RASA1_0'].value = 50000
    model.parameters['KRAS_%s_GTP_0' % mutant].value = total_pi

    sol = Solver(model, t)
    sol.run()
    plt.plot(t, sol.yobs['Pi_'] / total_pi, label=mutant)
    plt.ylabel('GTP hydrolyzed (%)')
    plt.ylim(top=1)
    plt.xlabel('Time (s)')
    plt.title('GAP-mediated hydrolysis')

plt.legend(loc='upper right', fontsize=11, frameon=False)


