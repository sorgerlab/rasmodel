from pysb.macros import bind_table

for ras_name in ['KRAS', 'NRAS', 'HRAS']:
    Monomer(ras_name,
            ['gtp', 'p_loop', 'switch1', 'switch2', 'CAAX', 'oncogenic'],
            {'oncogenic': ['y', 'n']})
Monomer('GTP', ['p'])
Monomer('GDP', ['p'])



def ras_binds_gtp_and_gdp(ras):
    bind_table([[     GTP,  GDP],
                [ras,   1,    1]], 'gtp', 'p', kf=1e-3)

def ras_converts_gtp_to_gdp(ras):
    k = Parameter('k_{0}_gtpase'.format(ras.name), 1.)
    Rule('{0}_converts_GTP_GDP'.format(ras.name),
         ras(gtp=1) % GTP(p=1) >>
         ras(gtp=1) % GDP(p=1),
         k)

def ras_interactions(ras):
    ras_binds_gtp_and_gdp(ras)
    ras_converts_gtp_to_gdp(ras)

ras_interactions(KRAS)



# A key thing to note here is that the mutations in G12, G15, and K16 appear
# to affect the affinity of Ras for GTP and GDP, not the catalytic rate.

# Unlike the mutations in G12 and its neighbors, which seem to affect
# activity by affecting GTP/GDP binding, the reduced activity resulting
# from mutations in Q61 appear to be attributed to an affect on the catalytic
# rate.

# As an implementation detail, note that the mutant rate should be constrained
# to be less than the wild type rate through the use of an Expression
# incorporating a scaling parameter between [0, 1].

Parameter('k_mut_gtpase', 0.1)

# Mutant Ras has diminished GTPase activity:
for ras in [KRAS, HRAS, NRAS]:
    ras_mut = ras(oncogenic='y')

    Rule('{0.name}_mut_converts_GTP_GDP'.format(ras),
         ras_mut(gtp=1) % GTP(p=1) >>
         ras_mut(gtp=1) % GDP(p=1),
         k_mut_gtpase)



# Add autophosphorylation of Ras A59T if it later turns out to be significant.






