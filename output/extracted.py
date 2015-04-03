for ras_name in ['KRAS', 'NRAS', 'HRAS']:
    Monomer(ras_name,
            ['gtp', 'p_loop', 'switch1', 'switch2', 'CAAX', 'oncogenic'],
            {'oncogenic': ['y', 'n']})
Monomer('GTP', ['p'])
Monomer('GDP', ['p'])

Parameter('k_gtp_bind', 1.)
Parameter('k_gtp_unbind', 1.)
Parameter('k_gdp_bind', 1.)
Parameter('k_gdp_unbind', 1.)
Parameter('k_gtpase', 1.)
for ras in [KRAS, HRAS, NRAS]:
    ras_wt = ras(oncogenic='n')
    Rule('{0.name}_binds_GTP'.format(ras),
         ras_wt(gtp=None) + GTP(p=None) <>
         ras_wt(gtp=1) % GTP(p=1),
         k_gtp_bind, k_gtp_unbind)
    Rule('{0.name}_converts_GTP_GDP'.format(ras),
         ras_wt(gtp=1) % GTP(p=1) >>
         ras_wt(gtp=1) % GDP(p=1),
         k_gtpase)
    Rule('{0.name}_binds_GDP'.format(ras),
         ras_wt(gtp=None) + GDP(p=None) <>
         ras_wt(gtp=1) % GDP(p=1),
         k_gdp_bind, k_gdp_unbind)

# Get identities of other mutants from reference 46.

# How to indicate necessity vs. sufficiency?






