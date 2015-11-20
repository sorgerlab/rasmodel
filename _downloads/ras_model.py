# Preliminaries
from pysb import *
from pysb.macros import bind, bind_table, equilibrate
Model()

# Define the site structure for various Ras family members.
# All of the Ras proteins have the following structural/regulatory features:
for ras_name in ['KRAS', 'NRAS', 'HRAS']:
    Monomer(ras_name,
            ['gtp', 'gef', 'p_loop', 's1s2', 'CAAX', 'oncogenic'],
            {'s1s2': ['closed', 'open'],
             'oncogenic': ['n', 'y']})

# Guanine nucleotides
Monomer('GTP', ['p', 'label'], {'label': ['n', 'y']})
Monomer('GDP', ['p', 'label'], {'label': ['n', 'y']})
Monomer('Pi', []) # inorganic phosphate, used for GDP/GTP recycling



def ras_binds_gxp(ras, gxp, klist):
    # Alias for Ras bound to GXP. Note that these equilibria are for
    # Ras that is NOT bound to a Ras-GEF.
    rasgxp = ras(gef=None, gtp=99) % gxp(p=99)
    # Get the rates from the list:
    (kf1, kr1, kf2, kr2) = klist

    bind(ras(gtp=None, s1s2='closed'), 'gtp', gxp(), 'p', [kf1, kr1])
#

    equilibrate(rasgxp(s1s2='closed'), rasgxp(s1s2='open'), [kf2, kr2])
#



ras_gdp_kf1 = 1e7   # M^-1 s^-1
ras_gdp_K1 = 5.7e4  # M^-1
ras_gdp_kr1 = ras_gdp_kf1 / ras_gdp_K1  # s^-1

ras_gdp_kf2 = 14.8   # s^-1
ras_gdp_kr2 = 1.8e-5 # s^-1

ras_gtp_kf1 = 1e7    # M^-1 s^-1
ras_gtp_K1 = 1.25e5  # M^-1
ras_gtp_kr1 = ras_gtp_kf1 / ras_gtp_K1  # s^-1

ras_gtp_kf2 = 16.7   # s^-1

ras_gdp_klist = [ras_gdp_kf1, ras_gdp_kr1, ras_gdp_kf2, ras_gdp_kr2]
ras_binds_gxp(HRAS, GDP, ras_gdp_klist)

ras_gtp_klist = [ras_gtp_kf1, ras_gtp_kr1, ras_gtp_kf2, ras_gdp_kr2]
ras_binds_gxp(HRAS, GTP, ras_gtp_klist)

#ras_binds_gxp(KRAS, GDP, ras_gdp_klist)
#ras_binds_gxp(KRAS, GTP, ras_gtp_klist)
#ras_binds_gxp(NRAS, GDP, ras_gdp_klist)
#ras_binds_gxp(NRAS, GTP, ras_gtp_klist)





def ras_converts_gtp_to_gdp(ras, kcat):
    k = Parameter('k_{0}_gtpase'.format(ras.name), 1.)
    # Instantiate the rule for both labeled and unlabeled GTP/GDP
    Rule('{0}_converts_GTP_GDP'.format(ras.name),
         ras(gef=None, gtp=1, s1s2='open') % GTP(p=1, label='n') >>
         ras(gef=None, gtp=1, s1s2='open') % GDP(p=1, label='n') + Pi(),
         k)
    Rule('{0}_converts_mGTP_mGDP'.format(ras.name),
         ras(gef=None, gtp=1, s1s2='open') % GTP(p=1, label='y') >>
         ras(gef=None, gtp=1, s1s2='open') % GDP(p=1, label='y') + Pi(),
         k)



# Convert 2.8e-2 min^-1 to units of s^-1
wt_ras_hydrolysis_rate = 2.8e-2 * 60

ras_converts_gtp_to_gdp(HRAS, wt_ras_hydrolysis_rate)



def recycle_gtp_from_gdp():
    k = Parameter('k_recycle_gtp_from_gdp', 1e7)
    # Note that only unbound GDP can be recycled!
    Rule('recycle_gtp_from_gdp_rule',
         GDP(p=None, label='n') + Pi() >> GTP(p=None, label='n'), k)
    Rule('recycle_mgtp_from_mgdp_rule',
         GDP(p=None, label='y') + Pi() >> GTP(p=None, label='y'), k)

recycle_gtp_from_gdp()





# A key thing to note here is that the mutations in G12, G15, and K16 appear
# to affect the affinity of Ras for GTP and GDP, not the catalytic rate.

# Unlike the mutations in G12 and its neighbors, which seem to affect
# activity by affecting GTP/GDP binding, the reduced activity resulting
# from mutations in Q61 appear to be attributed to an affect on the
# catalytic rate.

# As an implementation detail, note that the mutant rate should be
# constrained to be less than the wild type rate through the use of an
# Expression incorporating a scaling parameter between [0, 1].

# Add autophosphorylation of Ras A59T if it later turns out to be
# significant.








# Declare a list of RasGEFs along with their site structure.
# The names in the list below are HGNC standard names.
# (note: Cdc25Mm = RASGRF1)
ras_gef_names = ['SOS1', 'SOS2', 'RASGRF1', 'RASGRF2']
for ras_gef_name in ras_gef_names:
    Monomer(ras_gef_name, ['rasgef'])



def ras_gef_exchange_cycle(ras, rasgef, gxp,
                           k2_list, k3_list, k4a_list, k4b_list):
    # Alias for Ras bound to GXP
    rasgxp = ras(gef=None, gtp=99) % gxp(p=99)

    # Binding of RasGEF to nucleotide-free Ras (K2)
    bind(ras(gtp=None, s1s2='closed'), 'gef', rasgef(), 'rasgef', k2_list)

    # Binding of RasGEF to RasGXP (K3)
    bind(rasgxp(s1s2='open'), 'gef', rasgef(), 'rasgef', k3_list)

    # Binding of GXP to Ras/RasGEF complex
    bind(ras(s1s2='closed', gef=1) % rasgef(rasgef=1), 'gtp',
         gxp(), 'p', k4a_list)

    # Isomerization of Ras-RasGEF-GXP from loose to tight
    equilibrate(rasgxp(gef=1, s1s2='closed') % rasgef(rasgef=1),
                rasgxp(gef=1, s1s2='open') % rasgef(rasgef=1), k4b_list)



# Binding of RasGEF to nucleotide-free Ras
kf2 = 0.33e6        # M^-1 s^-1
kr2 = 1e-3          # s^-1

# Binding of RasGEF to RasGXP
KD3 = 0.6e-3        # M
kf3 = 3.4e4         # M^-1 s^-1 (lower limit)
kr3 = KD3 * kf3     # s^-1

# Binding of GXP to Ras/RasGEF complex
KD4a = 8.6e-6       # M
kf4a = 1e7          # M^-1 s^-1
kr4a = KD4a * kf4a  # s^-1

# Isomerization of Ras-RasGEF-GXP from loose to tight
kf4b = 20.4         # s^-1
kr4b = 3.9          # s^-1

#ras_gef_exchange_cycle(HRAS, RASGRF1, GTP, GDP)






#  p120GAP = RASA1.
# CaLB (calcium lipid binding domain) is also known as a C2 domain.
Monomer('RASA1', ['SH2_1', 'SH3', 'SH2_2', 'PH', 'C2', 'rasgap'])



Monomer('NF1', ['rasgap', 'CRALTRIO'])



# GAP1m = RASA2
Monomer('RASA2', ['C2_1', 'C2_2', 'rasgap', 'PH', 'BTK'])



# GAPIII = GAP1IP4BP = RASA3
Monomer('RASA3', ['C2_1', 'C2_2', 'rasgap', 'PH', 'BTK'])







# Iterate over every monomer
for m in model.monomers:
    states_dict = {}
    # Iterate over every site in the monomer
    for s in m.sites:
        # If it's in the site states dict, assign it the first of the
        # listed states
        if s in m.site_states:
            states_dict[s] = m.site_states[s][0]
        # Otherwise (e.g., the site is used only for binding) assign it
        # a state None, meaning unbound:
        else:
            states_dict[s] = None

    # Create the initial condition parameter based on the protein name
    initial_value = Parameter('{0}_0'.format(m.name), 1.0e-8)

    # Create the initial condition
    Initial(m(**states_dict), initial_value)



GTP_0.value = 468e-6

GDP_0.value = GTP_0.value / 10.

Parameter('mGTP_0', 0.)
Initial(GTP(p=None, label='y'), mGTP_0)

Parameter('mGDP_0', 0.)
Initial(GDP(p=None, label='y'), mGDP_0)






Observable('HRAS_GTP_', HRAS(gtp=1) % GTP(p=1))
Observable('HRAS_mGTP_', HRAS(gtp=1) % GTP(p=1, label='y'))
Observable('HRAS_GTP_closed_', HRAS(gtp=1, s1s2='closed') % GTP(p=1))
Observable('HRAS_GTP_open_', HRAS(gtp=1, s1s2='open') % GTP(p=1))
Observable('HRAS_GDP_', HRAS(gtp=1) % GDP(p=1))
Observable('HRAS_mGDP_', HRAS(gtp=1) % GDP(p=1, label='y'))
Observable('HRAS_GDP_closed_', HRAS(gtp=1, s1s2='closed') % GDP(p=1))
Observable('HRAS_GDP_open_', HRAS(gtp=1, s1s2='open') % GDP(p=1))
Observable('HRAS_nf_', HRAS(gtp=None))

Observable('GTP_', GTP())
Observable('GDP_', GDP())
Observable('Pi_', Pi())




