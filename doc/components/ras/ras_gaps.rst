Ras regulation by RasGAPs
=========================

GTP hydrolysis by Ras is slow but is accelerated by RasGAPs.

    [PMID9247124]_: The GTP-binding proteins return to the inactive state by
    virtue of the GTPase reaction, which is usually very slow but can be
    accelerated by the action of GAPs, in the case of the Ras/Ras-GAPs and
    Ran/Ran-GAP interactions by several orders of magnitude [PMID1569940]_
    [PMID8262937]_ [PMID7548002]_.


RasGAPs for the Ras subfamily
-----------------------------

There are several GAPs for the Ras subfamily.

    [PMID9247124]_: Five mammalian GAPs for Ras have been described. The first,
    p120GAP, is the prototype of this class of proteins and was the first one
    to be isolated (20, 21 and 22). Apart from being a regulator of Ras, its
    N-terminal domain contains a number of signalling modules such as SH2, SH3,
    PH, Calb domains and is believed to be a signal transduction molecule that
    may act independently of Ras (23 and 24).

RASA1
~~~~~

The domain information for RASA1 was taken from `Uniprot (ID: P20936)
<http://www.uniprot.org/uniprot/P20936>`_.

::

    def rasgap_monomers():

        #  p120GAP = RASA1.
        # CaLB (calcium lipid binding domain) is also known as a C2 domain.
        Monomer('RASA1', ['SH2_1', 'SH3', 'SH2_2', 'PH', 'C2', 'rasgap'])

NF1
~~~

    [PMID9247124]_: The second Ras-specific GAP is neurofibromin (NF1), which
    is the product of the neurofibromatosis gene [PMID25] and has also been
    shown to stimulate the GTPase of Ras 26, 27 and 28.  This gene has been
    found to be frequently mutated in patients with the disease
    neurofibromatosis type I 29, 30 and 31 but also, albeit less frequently, in
    solid tumors [PMID32].

`Uniprot (ID: P21359) <http://www.uniprot.org/uniprot/P21359>`_.

::

        Monomer('NF1', ['rasgap', 'CRALTRIO'])
    #

RASA2
~~~~~

    [PMID9247124]_: GAP1m, a mammalian homologue of the Drosophila GAP1 gene
    [PMID33], has been described, and a close homologue GAPIII [PMID34], both
    of which contain, in addition to the GRD (GAP-related domain), C2 domains
    and a PH domain.

`Uniprot (ID: Q15283) <http://www.uniprot.org/uniprot/Q15283>`_.

::

        # GAP1m = RASA2
        Monomer('RASA2', ['C2_1', 'C2_2', 'rasgap', 'PH', 'BTK'])
    #

RASA3
~~~~~

NOTE: As it turned out GAPIII and GAP1IP4BP are the same protein.

    [PMID9247124]_: Recently an inositol-4-phosphate (IP4) binding protein
    GAP1IP4BP has been purified, cloned, and found to contain a Ras-GAP
    catalytic domain.  In contrast to the other GAP mentioned, which are
    specific for Ras, GAP1IP4BP stimulates the GTPase of both Ras and Rap
    [PMID35].

`Uniprot (ID: Q14644) <http://www.uniprot.org/uniprot/Q14644>`_.

::

        # GAPIII = GAP1IP4BP = RASA3
        Monomer('RASA3', ['C2_1', 'C2_2', 'rasgap', 'PH', 'BTK'])
    #

Mechanism of accelerating hydrolysis
------------------------------------


    [9247124]_: Under saturating conditions of pl20GAP, the :math:`k_{cat}` of
    the GTPase reaction of Ras is 19 :math:`s^{-1}`, which means that the
    reaction is stimulated more than :math:`10^5`-fold, with a KD of 9.7 uM for
    the GAP/Ras-GTP interaction [1569940]_.

Because the co-crystal structure of Ras with RasGAP shows the occlusion of the
GTP binding pocket by RasGAP, we model this by asserting that the RasGAP binds
at the Switch 1/Switch 2 region on Ras. This has the effect of preventing the
binding of downstream effectors such as Raf.

Further, for the time being we assume that RasGAPs bind only to the GTP-bound
state, in which s1s2 are 'closed', i.e., once the nucleotide is in the tightly
bound state; and that a RasGAP and a RasGEF cannot simultaneously be bound to
the same Ras molecule.

These assumptions yield the following mechanism::

    def rasgap_mediated_hydrolysis(model, ras, rasgap, kf, kr, kcat):
        GTP = model.monomers['GTP']
        GDP = model.monomers['GDP']
        Pi = model.monomers['Pi']
        # RasGAP, with an unbound RasGAP domain, binds to GTP-bound Ras with
        # its Switch 1/Switch 2 region in the 'closed' state, with no RasGEF
        # bound:
        ras_gtp = ras(s1s2='closed', gef=None, gtp=2) % GTP(p=2)
        ras_gdp = ras(s1s2='closed', gef=None, gtp=2) % GDP(p=2)
        bind(rasgap, 'rasgap', ras_gtp, 'gap', [kf, kr])
        # Because a GDP is being created on the right-hand side of this rule,
        # it must be fully specified. This means that we need to handle
        # the case when the GTP/GDP is either labeled or unlabeled.
        for label in ['y', 'n']:
            _macro_rule('catalyze',
                        rasgap(rasgap=1) % ras_gtp(gap=1, label=label) >>
                        rasgap(rasgap=None) +
                            ras_gdp(gap=None, label=label) + Pi(),
                        [kcat], ['kcat'])

For the time being, we include only the hydrolysis reaction for KRAS,
mediated by RASA1 (P120GAP)::

    def kras_rasgaps(model):
        KRAS = model.monomers['KRAS']
        RASA1 = model.monomers['RASA1']

        # Diffusion-limited KRAS/RasGAP association:
        kf = 1e-2
        # Calculate dissociation rate based on roughly 10 uM KD reported above
        kr = 100. # KD = 100 / 0.01 = 10,000 nM = 10 uM

        hydrolysis_rates = {'WT': 4300e-5,
                            'G12A': 32e-5,
                            'G12C': 20e-5,
                            'G12D': 89e-5,
                            'G12R': 20e-5,
                            'G12V': 24e-5,
                            'G13D': 20e-5,
                            'Q61L': 12e-5,
                            'Q61H': 5e-5,}
        # Iterate over all of the mutants that we're considering
        for mutant in KRAS.site_states['mutant']:
            kras = KRAS(mutant=mutant)
            if mutant in hydrolysis_rates:
                mutant_kcat = hydrolysis_rates[mutant]
            else:
                mutant_kcat = hydrolysis_rates['WT']
            rasgap_mediated_hydrolysis(model, kras, RASA1, kf, kr, mutant_kcat)


