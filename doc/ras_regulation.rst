Ras activity and regulation
===========================

Anatomy of Ras
--------------

For a diagram of the primary structure of Ras, see
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915522/figure/F1/).::

    from pysb.macros import bind, bind_table, equilibrate

    # All of the Ras proteins have the following structural/regulatory features:
    for ras_name in ['KRAS', 'NRAS', 'HRAS']:
        Monomer(ras_name,
                ['gtp', 'gef', 'p_loop', 's1s2', 'CAAX', 'oncogenic'],
                {'s1s2': ['closed', 'open'],
                 'oncogenic': ['y', 'n']})

    # Guanine nucleotides
    Monomer('GTP', ['p'])
    Monomer('GDP', ['p'])

A note on alternative splice variants:

    [3304147]_: The three functional ras genes code for highly related proteins
    generically known as p21 (71). The p21 coding sequences of each of these
    genes are equally distributed in four exons except for the K-ras-2 gene,
    which possesses two alternative fourth coding exons (exons IVA and IVB)
    that allow the synthesis of two isomorphic p21 proteins of 188 and 189
    residues that differ in their carboxy terminal domains (54, 60, 72).

Ras proteins are GTPases
------------------------

    [3304147]_: ras proteins, independently of their phylogenetic origin, have
    been shown to bind guanine nucleotides (GTP and GDP) ([3304147_22]_
    [3304147_23]_ [3304147_24]_ [3304147_25]_) and possess intrinsic GTPase
    activity ([3304147_25]_ [3304147_26]_ [6147754]_ [6148703]_ [3304147_29]_)

::

    #def ras_interactions(ras):
    #    ras_binds_gtp_and_gdp(ras)
    #    ras_converts_gtp_to_gdp(ras)

Ras binds GTP and GDP
~~~~~~~~~~~~~~~~~~~~~

Mechanism and rates associated with GTP/GDP binding by Ras.

The following statements were taken from a kinetic analysis of Ras and
nucleotide interactions. All rates were measured at 20C.

    [9585556]_: the intrinsic dissociation rate of Ras for GTP (1 × 10-5 s-1) is
    2-fold lower than that for GDP (2 × 10-5 s-1)...

    [9585556]_: Numerically, it was more convenient to use the corresponding
    differential equations with the program FACSIMILE and to calculate for 1000
    s with the assumption of fast association rate constants (in all cases:
    10^7 M-1 s-1).

    [9585556]_: The equilibrium dissociation constant for Ras-3′mdGDP (KD1) had
    been determined independently as 9 pM from nucleotide association and
    dissociation experiments (Tables 2 and 3).

::

    # Make the same assumptions about fast association rate constants as
    # 9585556, and use their measured values for dissociation rates in the
    # absence of GEFs (association rates are in units of nM^-1 s^-1).
    def ras_binds_gtp_and_gdp(ras, k_gtp_diss, k_gdp_diss):
        kf = 1e-2
        bind_table([[                   GTP,               GDP],
                    [ras,  (kf, k_gtp_diss),  (kf, k_gdp_diss)]], 'gtp', 'p')

    # The data in Table 1 gives a value of 1.2e-5 for the dissociation rate with
    # GDP, whereas the text gives rates of 1e-5 and 2e-5 for GTP/GDP,
    # respectively.
    ras_binds_gtp_and_gdp(HRAS, 1e-5, 2e-5)

    # The rate for KRAS/GDP association is given in Table 1 as 1.6e-5, but the
    # KRAS/GTP rate is not measured.
    ras_binds_gtp_and_gdp(KRAS, 8e-6, 1.6e-5)

    # The rate for NRAS/GDP association is given in Table 1 as 1.0e-5, but the
    # NRAS/GTP rate is not measured.
    ras_binds_gtp_and_gdp(NRAS, 5e-6, 1.0e-5)

Ras converts GTP to GDP
~~~~~~~~~~~~~~~~~~~~~~~

GTP hydrolysis by wild-type Ras is very slow in the absence of RasGAPs.

    [1569940]_: It has been reported that the in vitro GTPase activity of wild-type
    p21, which proceeds at a rate of 0.028 min^-1 at 37°C ([2502546]_), is
    accelerated 100- to 200-fold by GAP, as measured under nonsaturating
    conditions.

::

    # Convert 2.8e-2 min^-1 to units of s^-1
    wt_ras_hydrolysis_rate = 2.8e-2 * 60

    def ras_converts_gtp_to_gdp(ras):
        k = Parameter('k_{0}_gtpase'.format(ras.name), 1.)
        Rule('{0}_converts_GTP_GDP'.format(ras.name),
             ras(gtp=1) % GTP(p=1) >>
             ras(gtp=1) % GDP(p=1),
             k)

Ras regulation by GEFs
----------------------

Identities of Ras-specific GEFs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several guanosine exchange factors (GEFs) that are specific to Ras
subfamily proteins.

    [9585556]_: Several genes have been isolated from different organisms
    encoding proteins that have a GEF activity specific for Ras (for which we
    use the general name RasGEFs throughout this paper): SOS1 and SOS2
    ([9585556_1]_ [9585556_2]_ [9585556_3]_ [9585556_4]_) ; Cdc25Mm, also
    called RasGrf ([9585556_5]_ [9585556_6]_ [9585556_7]_); and mRas-GRF2
    ([9585556_8]_).

RasGEFs contain a specific domain responsible for activating Ras proteins.

    [9585556]_: The RasGEFs are proteins of considerable length, 120 - 160 kDa,
    and contain several regions which are generally accepted to represent
    structural domains (12). A region of 200 - 300 amino acids, the RasGEF
    domain, is shared by all GEFs which act on members of the Ras subfamily,
    and their activity is specific toward either Ras, Ral, or Rap. The fact
    that truncated versions of various lengths, containing this RasGEF domain,
    have been shown to be active RasGEFs in vivo and in vitro (4 , 13 - 16)
    confirms that this region indeed represents the Ras-specific guanine
    nucleotide exchange domain.

::

    # Declare a list of RasGEFs along with their site structure.
    # The names in the list below are HGNC standard names.
    # (note: Cdc25Mm = RASGRF1)
    ras_gef_names = ['SOS1', 'SOS2', 'RASGRF1', 'RASGRF2']
    for ras_gef_name in ras_gef_names:
        Monomer(ras_gef_name, ['rasgef'])

Mechanism of GEFs
~~~~~~~~~~~~~~~~~

Ras binds RasGEFs in the absence of nucleotides.

    [9690470]_: Biochemical studies of Ras exchange factors have shown that the
    complex of Ras with these proteins is stable in the absence of nucleotides
    and is dissociated by the rebinding of either GDP or GTP ([9585556]_
    [9690470_17]_ [9690470_18]_ [9690470_21]_ [9690470_22]_) The principal role
    for the exchange factor is to facilitate nucleotide release, and it does
    not seem to control significantly the preferential rebinding of GTP over
    GDP ([9585556]_, [9690470_22]_, [9690470_23]_).  Cellular concentrations of
    GTP are 10-fold higher than GDP, which results in the loading of GTP onto
    Ras.

The following study used purified HRAS and mouse RASGRF1:

    [9690470]_: The mechanism of nucleotide release by the catalytic domain of
    murine Cdc25 (Cdc25Mm) has been investigated recently using fluorescently
    labelled nucleotides [9585556]_.  The affinity of Cdc25Mm for
    nucleotide-free Ras (Kd = 4.6 nM) is found to be several orders of
    magnitude higher than that for nucleotide-bound Ras, and the maximal
    acceleration by Cdc25Mm of the rate of dissociation of nucleotide is more
    than 10^5.

    [9585556]_: The best fit of our data resulted in similar quantum yields and
    a value of 4.6 nM for KD2 (NOTE: Kd between nucleotide-free H-Ras and
    RasGRF1). A variation in the value for KD2 of approximately 2-fold resulted
    in fits of comparable quality.

The activity of GEF (RASGRF1 in this case) does not depend on whether Ras
(HRAS) is loaded with GTP or GDP.

    [9585556]_: However, since the intrinsic dissociation rate of Ras for GTP
    (1 × 10-5 s-1) is 2-fold lower than that for GDP (2 × 10-5 s-1), the
    stimulatory action of Cdc25Mm285 is practically independent of the nature
    of the bound nucleotide.

    [9585556]_: Although we did not reach complete saturation at 600 μM
    Ras‚nucleotide, the data could be fitted to obtain a maximal rate of
    3′mdGDP release from Ras of 3.9 s-1 and an apparent Km value of 386 μM.
    Since the intrinsic dissociation rate of 3′mdGDP is 2 × 10-5 s-1 (Table 1),
    the acceleration of GDP dissociation from Ras by this GEF is approximately
    2 × 105-fold. An apparent Km of approximately 300 μM was obtained for the
    triphosphate-bound form of Ras, confirming that there is no pronounced
    specificity toward the nature of the Ras-bound nucleotide (data not shown).

.. warning:: GEF binding to GTP bound Ras?

    Can GEFs bind to Ras and cause ejection of nucleotide before the GTP/GDP
    conversion is complete? Moreover, if GEF binds to Ras-GTP, can the
    hydrolysis to GDP proceed while GEF is bound?

The RasGEF exchange cycle
~~~~~~~~~~~~~~~~~~~~~~~~~

The following reaction scheme for the GEF exchange cycle, along with the
associated rates, are drawn from [9585556]_.

.. image:: /images/9585556_rasgef_cycle.png
    :width: 600px

::

    def ras_gef_exchange_cycle(ras, rasgef, gxp):
        # An alias for Ras bound to GXP
        rasgxp = ras(gef=None, gtp=99) % gxp(p=99)

        # Nucleotide-free Ras binds GTP/GDP
        # KD1a is given as 11.8 uM; we calculate the off-rate assuming
        # a fast on rate of 1e7 M^-1 s^-1.
        KD1a = 11.8e-6
        kf1a = 1e7
        kr1a = KD1a * kf1a
        bind(ras(gtp=None, s1s2='closed'), 'gtp', gxp(), 'p', [kf1a, kr1a])

        # Isomerization/conformational change of Ras resulting from nucleotide
        # binding; also described as the conversion of the nucleotide from
        # loosely bound to tightly bound.
        kf1b = 26.8
        kr1b = 20e-6
        equilibrate(rasgxp(s1s2='closed'), rasgxp(s1s2='open'), [kf1b, kr1b])

        # Binding of RasGEF to nucleotide-free Ras
        kf2 = 0.33e6
        kr2 = 1e-3
        bind(ras(gtp=None, s1s2='closed'), 'gef', rasgef(), 'rasgef',
             [kf2, kr2])

        # Binding of RasGEF to RasGXP
        KD3 = 0.6e-3
        kf3 = 3.4e4 # Lower limit
        kr3 = KD3 * kf3
        bind(rasgxp(s1s2='open'), 'gef', rasgef(), 'rasgef', [kf3, kr3])

        # Binding of GXP to Ras/RasGEF complex
        KD4a = 8.6e-6
        kf4a = kf1a # on rate is insensitive to presence of GEF
        kr4a = KD4a * kf4a
        bind(ras(s1s2='closed', gef=50) % rasgef(rasgef=50), 'gtp',
             gxp(), 'p', [kf4a, kr4a])

        # Isomerization of Ras-RasGEF-GXP from loose to tight
        kf4b = 20.4
        kr4b = 3.9
        equilibrate(rasgxp(gef=1, s1s2='closed') % rasgef(rasgef=1),
                    rasgxp(gef=1, s1s2='open') % rasgef(rasgef=1), [kf4b, kr4b])

Instantiate the RasGEF cycle for HRAS and RASGRF1::

    ras_gef_exchange_cycle(HRAS, RASGRF1, GTP)
    ras_gef_exchange_cycle(HRAS, RASGRF1, GDP)

.. warning:: How does GTP hydrolysis fit into the cycle?

    Can Ras hydrolyze GTP to GDP at any point in this cycle? Or can this only
    happen when Ras is bound to GDP and GEF is not bound? Does it only happen
    when nucleotide is in the tightly bound conformation?

[9585556]_: Therefore, we tested the nucleotide specificity of the interaction
of Cdc25Mm285 (CdcMm285 is the fragment of CdcMm/RasGRF1 containing the RasGEF
domain) with Ras. Figure 1 shows the release of Ras-bound 3′mdGDP or 3′mdGTP (4
μM), in the presence of an excess of unlabeled nucleotide and in the presence
or absence of 1 μM Cdc25Mm285. The Cdc25Mm285-stimulated dissociation rate of
Ras-3′mdGDP is approximately twice that of Ras-3′mdGTP, with values of 0.0098
and 0.0046 s-1, respectively.  However, since the intrinsic dissociation rate
of Ras for GTP (1 × 10-5 s-1) is 2-fold lower than that for GDP (2 × 10-5 s-1),
the stimulatory action of Cdc25Mm285 is practically independent of the nature
of the bound nucleotide. The difference in stimulated dissociation rates is
somewhat smaller than the results of Jacquet et al. (16) but is similar to the
results with the yeast proteins CDC25 and RAS2 obtained by Haney and Broach
(28).

[9690470]_: Kinetic analysis of nucleotide association shows that the reaction
proceeds by the formation of a ternary complex of a loosely bound nucleotide
and Ras – Cdc25Mm followed by conversion to a form in which the nucleotide is
tightly bound to Ras [9585556]_. In light of the structure of the Ras–Sos
complex, the first step can be interpreted as the interaction of the base and
the ribose of the nucleotide with the part of the Ras binding site that is not
occluded by Sos. The second step would involve a conformational change in the
Switch 2 segment and release of Switch 1, resulting in the restructuring of a
competent binding site for phosphate and magnesium, and the subsequent
dissociation of Sos.

[9690470]_: As a nucleotide-exchange factor, Sos functions under two apparently
conflicting imperatives. The interaction between Sos and Ras must be strong
enough to dislodge the tightly bound nucleotide, but the Ras – Sos complex must
also be poised for subsequent displacement by incoming nucleotides. The
structure of the Ras – Sos complex shows that Ras and Sos meet these demands by
forming a tight complex that is anchored at one end of the nucleotide- binding
site, where phosphate and magnesium are normally bound. The interface between
Sos and Ras is mainly hydrophilic, suggesting a ready unzippering through
water-mediated displacements of the coordinating side chains. The main
interacting elements of Sos avoid direct occlusion of the nucleotide-binding
site, except the region where the terminal phosphate groups and the magnesium
ion are bound. This feature allows incoming nucleotides to reverse the process
by competing for the groups that ligate the phosphate and metal ion.

[9690470]_: The overall shape of the catalytic domain of Sos is that of an
oblong bowl (Fig. 2), with Ras bound at the centre of the bowl. The regions of
Ras that interact most closely with Sos include the phosphate-binding P-loop
(residues 10 – 17) and surrounding segments (including strand 􏰧1 and helix 􏰦1),
the Switch 1 region (defined here as residues 25–40) and the Switch 2 region
(defined here as residues 57 – 75). Additional interactions are seen with helix
3 (residues 95–105; Fig. 3a, b). The interface between Ras and Sos is primarily
hydrophilic and very extensive, with 3,600 A^2 of surface area buried in the
complex.

[9690470]_: The most obvious effect of Sos binding to Ras is the opening of the
nucleotide binding site as a result of the displacement of Switch 1 of Ras by
the insertion of the helical hairpin formed by aH and aI of Sos (Fig. 5)

Switch 1 and Switch 2 are the only regions of Ras in which structural changes
are directly induced by Sos.

The change in the Switch 1 region of Ras when bound to Sos is drastic...Switch
1 is completely removed from the nucleotide-binding site.

One important aspect of the insertion of the helical hairpin of Sos into the
Switch 1 region is that it does not result in a significant occlusion of the
guanine and ribose binding sites (Fig. 5d). Instead, this structural distortion
breaks the network of direct and water-mediated interactions between Switch 1
and the nucleotide. For example, in the nucleotide-bound forms of Ras, Phe 28
interacts with the guanine base through a perpendicular aromatic – aromatic
interaction (Fig. 5a). Mutation of Phe28 to leucine results in a significant
increase in the intrinsic rate of dissociation of nucleotide from Ras18. In the
Sos complex, the Calpha of Phe 28 moves 9.6 A and the side chain no longer
interacts with the nucleotide-binding site (Fig. 5b).

The Switch 2 region of Ras makes important interactions with GTP and not with
GDP (19,46). Nevertheless, structural changes that are induced in Switch 2 by
Sos result in the exclusion of both GDP and GTP, because they affect magnesium
binding as well as the conformation of Lys 16 in the P-loop, a crucial
phosphate ligand.

Specificity of RASGRF1 for Ras isoforms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[9585556]_: Three mammalian isoforms of Ras, H-, K-, and N-Ras, have been
identified which are highly conserved intheirprimarysequence.
Thesignificanceofhavingmore than one isoform is not understood at present,
although the isoforms may have different functions in different tissues, since
certain types of tumors have a preference for a particular activated Ras gene,
such as K-Ras for lung, colon and pancreas cancers and N-Ras for myeloid
leukemias (25). To see whether Cdc25Mm285 acts differently on the three
isoforms, we tested the GEF activity of Cdc25Mm285 on these proteins. As
summarized in Table 1, Cdc25Mm285 is active on all isoforms, being somewhat
more active on N-Ras, in accordance with the results of Leonardsen et al. (26).

Ras regulation by RasGAPs
-------------------------

GTP hydrolysis by Ras is slow but is accelerated by RasGAPs.

    [9247124]_: The GTP-binding proteins return to the inactive state by virtue
    of the GTPase reaction, which is usually very slow but can be accelerated
    by the action of GAPs, in the case of the Ras/Ras-GAPs and Ran/Ran-GAP
    interactions by several orders of magnitude [1569940]_ [8262937]_
    [7548002]_.

There are several GAPs for the Ras subfamily.

    [9247124]_: Five mammalian GAPs for Ras have been described. The first,
    p120GAP, is the prototype of this class of proteins and was the first one
    to be isolated 20, 21 and 22. Apart from being a regulator of Ras, its
    N-terminal domain contains a number of signalling modules such as SH2, SH3,
    PH, Calb domains and is believed to be a signal transduction molecule that
    may act independently of Ras 23 and 24.

RASA1
~~~~~

The domain information for RASA1 was taken from `Uniprot (ID: P20936)
<http://www.uniprot.org/uniprot/P20936>`_.

::

    #  p120GAP = RASA1.
    # CaLB (calcium lipid binding domain) is also known as a C2 domain.
    Monomer('RASA1', ['SH2_1', 'SH3', 'SH2_2', 'PH', 'C2', 'rasgap'])

..

NF1
~~~

    [9247124]_: The second Ras-specific GAP is neurofibromin (NF1), which is
    the product of the neurofibromatosis gene [25] and has also been shown to
    stimulate the GTPase of Ras 26, 27 and 28.  This gene has been found to be
    frequently mutated in patients with the disease neurofibromatosis type I
    29, 30 and 31 but also, albeit less frequently, in solid tumors [32].

`Uniprot (ID: P21359) <http://www.uniprot.org/uniprot/P21359>`_.

::

    Monomer('NF1', ['rasgap', 'CRALTRIO'])

RASA2
~~~~~

    [9247124]_: GAP1m, a mammalian homologue of the Drosophila GAP1 gene [33],
    has been described, and a close homologue GAPIII [34], both of which
    contain, in addition to the GRD (GAP-related domain), C2 domains and a PH
    domain.

`Uniprot (ID: Q15283) <http://www.uniprot.org/uniprot/Q15283>`_.

::

    # GAP1m = RASA2
    Monomer('RASA2', ['C2_1', 'C2_2', 'rasgap', 'PH', 'BTK'])

RASA3
~~~~~


NOTE: As it turned out GAPIII and GAP1IP4BP are the same protein.

    [9247124]_: Recently an inositol-4-phosphate (IP4) binding protein GAP1IP4BP
    has been purified, cloned, and found to contain a Ras-GAP catalytic domain.
    In contrast to the other GAP mentioned, which are specific for Ras,
    GAP1IP4BP stimulates the GTPase of both Ras and Rap [35].

`Uniprot (ID: Q14644) <http://www.uniprot.org/uniprot/Q14644>`_.

::

    # GAPIII = GAP1IP4BP = RASA3
    Monomer('RASA3', ['C2_1', 'C2_2', 'rasgap', 'PH', 'BTK'])

Oncogenic Ras mutants have reduced GTP binding and GTPase activity
-------------------------------------------------------------------

[18568040]_: In 1984, three groups reported that mutated Ras oncoproteins
differ functionally from their normal counterparts [6147754]_
[18568040_42]_ [6148703]_. The oncogenic forms of Ras exhibited impaired
GTPase activity, which suggested that the hydrolysis of GTP somehow terminates
the activated state of the protein, which is consistent with the presumed
analogy to the behaviour of G proteins...Furthermore, the link between the
much-studied Gly-to-Val substitution of residue 12 of H-Ras and GTP hydrolysis
was made the following year by Frank McCormick’s group, which noted that
antibodies that are specific to that region blocked GTP binding [18568040_44]_.

[3304147]_: Early studies have predicted that replacement of Gly12 by any other
amino acid residue (except proline) would disrupt the a-helical structure of
the amino terminal domain of ras proteins, causing a conformational change that
would prevent its proper folding (112-114). Thus, replacement or elimination of
Gly12 may create a rigid domain that cannot efficiently interact with the
phosphoryl region of the GTP molecule, reducing the GTPase activity of ras
proteins. Two additional residues in this domain, Glyl5 and Lysl6, are present
in other guanine nucleotide-bindingproteins(109, 111). Substitution of Lys16 by
Asn16 significantly reduces GTP/GDP affinity without affecting base
specificity, an observation consistent with the idea that these residues are
also part of the phosphoryl group (95)::

    # A key thing to note here is that the mutations in G12, G15, and K16 appear
    # to affect the affinity of Ras for GTP and GDP, not the catalytic rate.

[18568040]_: Other oncogenic mutations (such as Gln61leu in H-Ras) were
also shown to impair GTP hydrolysis [18568040_45]_ and other oncogenic forms of
Ras were later determined to be impaired in GTP hydrolysis (for example, REF.
[18568040_46]_).

[3304147]_: Substitution of Gln61 by 17 different amino acid residues
invariably results in decreased GTPase activity ([3304147_25]_, 117).

.. _FIG4a: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915522/figure/F4/
.. _FIG4b: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915522/figure/F4/

[18568040]_: The overall Ras structure was shown to consist of a
hydrophobic core of six stranded β-sheets and five α-helices that are
interconnected by a series of ten loops (FIG4a_). Five of these loops are
situated on one facet of the protein and have crucial roles in determining the
high affinity nucleotide interactions of Ras and in regulating GTP hydrolysis.
In particular, the GTP γ-phosphate is stabilized by interactions that are
established with the residues of loops 1, 2 and 4 (for example, lys16, Tyr32,
Thr35, Gly60 and Gln61; see FIG4b_). A prominent role is attributed to Gln61,
which stabilizes the transition state of GTP hydrolysis to GDP, in addition to
participating in the orientation of the nucleophilic attack that is necessary
for this reaction. As such, oncogenic mutations of Gln61 reduce the intrinsic
GTP hydrolysis rate, thereby placing the Ras protein in a constitutively active
state.::

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

Autophosphorylation of Ras A59T
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[3304147]_: In addition to GTP/GDP binding and GTPase activity, ras proteins
carrying an Ala59 -> Thr59 mutation exhibit an autophosphorylating activity of
an, as yet, unknown biological significance [3304147_23]_. In all cases, Thr59
has been found to be the phosphate receptor site (106). No transphosphorylating
activity has been detected with any ras protein, including those carrying Thr59
mutations::

    # Add autophosphorylation of Ras A59T if it later turns out to be significant.

Anatomy of Ras regulation
-------------------------

[18568040]_: The structural differences between the RasGDP and the RasGTP
conformations reside mainly in two highly dynamic regions, termed switch i
(residues 30–40) and switch ii (residues 60–76). Both regions are required for
the interactions of Ras with upstream as well as downstream partners (see also
FIG. 2a). The binding of GTP alters the conformation of switch i, primarily
through the inward reorientation of the side chain of Thr35, thereby enabling
its interactions with the GTP γ-phosphate as well as the Mg2+ ion. Similarly,
the γ-phosphate induces significant changes in the orientation of the switch ii
region through interactions it establishes with Gly60 (FIG. 4b).

Post-translational modifications of the C-terminus
--------------------------------------------------

An initial study in this area, published in 1982, showed that the mature form
of viral H-Ras localized to the cell membrane47. Several months later it was
demonstrated that viral H-Ras is palmitoylated at the C terminus; the resulting
attached lipid moiety facilitated its association with the membrane48. The
functional connection between this lipid modification and Ras function was made
by Douglas Lowy’s group in 1984, which showed that lipid binding and membrane
association were actually required for the transforming activity of the viral
H-Ras oncoprotein49,50.

working with cellular H-Ras, Stuart Aaronson’s group proceeded to demonstrate
that this C-terminal processing and membrane recruitment of Ras is a
prerequisite to its biochemical activation51.

The molecular mechanisms of Ras lipid processing were laid out over the
subsequent 5 years through a series of observations using yeast genetics,
protein biochemistry and in vitro cellular systems52–57 (FIGS 2,3).3).

Indeed, the C-terminal CAAX motif, previously found to be important for Ras
function, was found to be the target of a post-translational modification that
involved the addition of a farnesyl isoprenoid lipid, catalysed by the enzyme
farnesyl transferase (FTase).

Subsequent studies determined that this prenylation reaction is followed by the
proteolytic cleavage of the AAX sequence, catalysed by Ras-converting enzyme-1
(RCE1) and the carboxymethylation of the now terminal Cys residue by the
isoprenylcysteine carboxymethyltransferase-1 (ICMT1) enzyme.

Although these CAAX-signal modifications appeared to be essential for the
association of Ras with the plasma membrane, other studies identified the
requirement for a second C-terminal signal that facilitates full membrane
recruitment and hence full Ras function (for example, see REF. 57). For
K-Ras-4B, this second signal is a string of positively-charged lys residues
upstream of the C terminus that are sufficient to anchor the protein to the
membrane. However, prenylated H-Ras, N-Ras and K-Ras-4A require a further
palmitoylation step in which a palmitoyl moiety is attached to upstream
C-terminal Cys residues before their anchoring in the membrane is stabilized.


References
----------

.. [3304147] Barbacid M. **ras genes.** Annu Rev Biochem. 1987;56:779-827. Review. :pmid:`3304147`. :download:`PDF </pdf/3304147.pdf>`.

.. [3304147_22] Scolnick EM, Papageorge AG, Shih TY. **Guanine nucleotide-binding activity as an assay for src protein of rat-derived murine sarcoma viruses.** Proceedings of the National Academy of Sciences of the United States of America. 1979;76(10):5355-5359. :pmid:`228288`.

.. [3304147_23] Shih TY, Papageorge AG, Stokes PE, Weeks MO, Scolnick EM. **Guanine nucleotide-binding and autophosphorylating activities associated with the p21src protein of Harvey murine sarcoma virus.** Nature. 1980 Oct 23;287(5784):686-91. :pmid:`6253810`.

.. [3304147_24] Tamanoi, F., Walsh, M., Kataoka, T., & Wigler, M. (1984). **A product of yeast RAS2 gene is a guanine nucleotide binding protein.** Proceedings of the National Academy of Sciences of the United States of America, 81(22), 6924–6928. :pmid:`6438624`.

.. [3304147_25] Temeles GL, Gibbs JB, D'Alonzo JS, Sigal IS, Scolnick EM. **Yeast and mammalian ras proteins have conserved biochemical properties.** Nature. 1985 Feb 21-27;313(6004):700-3. :pmid:`3919305`.

.. [3304147_26] Gibbs JB, Sigal IS, Poe M, Scolnick EM. **Intrinsic GTPase activity distinguishes normal and oncogenic ras p21 molecules.** Proc Natl Acad Sci U S A. 1984 Sep;81(18):5704-8. :pmid:`6148751`.

.. [6147754] McGrath JP, Capon DJ, Goeddel DV, Levinson AD. **Comparative biochemical properties of normal and activated human ras p21 protein.** Nature. 1984 Aug 23-29;310(5979):644-9. :pmid:`6147754`.

.. [6148703] Sweet RW, Yokoyama S, Kamata T, Feramisco JR, Rosenberg M, Gross M. **The product of ras is a GTPase and the T24 oncogenic mutant is deficient in this activity.** Nature. 1984 Sep 20-26;311(5983):273-5. :pmid:`6148703`.

.. [3304147_29] Manne V, Bekesi E, Kung HF. **Ha-ras proteins exhibit GTPase activity: point mutations that activate Ha-ras gene products result in decreased GTPase activity.** Proc Natl Acad Sci U S A. 1985 Jan;82(2):376-80. :pmid:`2982154`.

.. [18568040] Karnoub AE, Weinberg RA. **Ras oncogenes: split personalities.** Nature reviews Molecular cell biology. 2008;9(7):517-531. :doi:`10.1038/nrm24381`. :pmid:`18568040`.

.. [18568040_42] Gibbs JB, Sigal IS, Poe M, Scolnick EM. **Intrinsic GTPase activity distinguishes normal and oncogenic ras p21 molecules.** Proc Natl Acad Sci USA. 1984;81:5704–5708. :pmid:`6148751`.

.. [18568040_44] Clark R, Wong G, Arnheim N, Nitecki D, McCormick F. **Antibodies specific for amino acid 12 of the ras oncogene product inhibit GTP binding.** Proc Natl Acad Sci USA. 1985;82:5280–5284.:pmid:`3927300`.

.. [18568040_45] Der CJ, Finkel T, Cooper GM. **Biological and biochemical properties of human rasH genes mutated at codon 61.** Cell. 1986;44:167–176. :pmid:`3510078`.

.. [18568040_46] Trahey M, McCormick F. **A cytoplasmic protein stimulates normal N-ras p21 GTPase, but does not affect oncogenic mutants.** Science.  1987;238:542–545. References 41–46 established that oncogenic mutation of ras affects its nucleotide cycle. :pmid:`2821624`.

.. [9690470] Boriack-Sjodin PA, Margarit SM, Bar-Sagi D, Kuriyan J. **The structural basis of the activation of Ras by Sos.** Nature. 1998 Jul 23;394(6691):337-43. :pmid:`9690470`. :download:`PDF </pdf/9690470.pdf>`.

.. [9585556] Lenzen C, Cool RH, Prinz H, Kuhlmann J, Wittinghofer A. **Kinetic analysis by fluorescence of the interaction between Ras and the catalytic domain of the guanine nucleotide exchange factor Cdc25Mm.** Biochemistry. 1998 May 19;37(20):7420-30. :pmid:`9585556`. :download:`PDF </pdf/9585556.pdf>`.

.. [9690470_17] Lai CC, Boguski M, Broek D, Powers S. **Influence of guanine nucleotides on complex formation between Ras and CDC25 proteins.** Mol Cell Biol. 1993 Mar;13(3):1345-52. :pmid:`8441380`.

.. [9690470_18] Mistou MY, Jacquet E, Poullet P, Rensland H, Gideon P, Schlichting I, Wittinghofer A, Parmeggiani A. **Mutations of Ha-ras p21 that define important regions for the molecular mechanism of the SDC25 C-domain, a guanine nucleotide dissociation stimulator.** EMBO J. 1992 Jul;11(7):2391-7. :pmid:`16286121`.

.. [9690470_21] Powers S, O'Neill K, Wigler M. **Dominant yeast and mammalian RAS mutants that interfere with the CDC25-dependent activation of wild-type RAS in Saccharomyces cerevisiae.** Mol Cell Biol. 1989 Feb;9(2):390-5. :pmid:`2651897`.

.. [9690470_22] Haney SA, Broach JR. **Cdc25p, the guanine nucleotide exchange factor for the Ras proteins of Saccharomyces cerevisiae, promotes exchange by stabilizing Ras in a nucleotide-free state.** J Biol Chem. 1994 Jun 17;269(24):16541-8. :pmid:`8206969`.

.. [9690470_23] Klebe C, Prinz H, Wittinghofer A, Goody RS. **The kinetic mechanism of Ran--nucleotide exchange catalyzed by RCC1. Biochemistry.** 1995 Oct 3;34(39):12543-52.:pmid:`7548002`.

.. [9585556_1] Rogge RD, Karlovich CA, Banerjee U. **Genetic dissection of a neurodevelopmental pathway: Son of sevenless functions downstream of the sevenless and EGF receptor tyrosine kinases.** Cell. 1991 Jan 11;64(1):39-48. :pmid:`1846090`.

.. [9585556_2] Bonfini L, Karlovich CA, Dasgupta C, Banerjee U. **The Son of sevenless gene product: a putative activator of Ras.** Science. 1992 Jan 31;255(5044):603-6. :pmid:`1736363`.

.. [9585556_3] Bowtell D, Fu P, Simon M, Senior P. **Identification of murine homologues of the Drosophila son of sevenless gene: potential activators of ras.** Proc Natl Acad Sci U S A. 1992 Jul 15;89(14):6511-5.  :pmid:`1631150`.

.. [9585556_4] Chardin P, Camonis JH, Gale NW, van Aelst L, Schlessinger J, Wigler MH, Bar-Sagi D. **Human Sos1: a guanine nucleotide exchange factor for Ras that binds to GRB2.** Science. 1993 May 28;260(5112):1338-43. :pmid:`8493579`.

.. [9585556_5] Martegani E, Vanoni M, Zippel R, Coccetti P, Brambilla R, Ferrari C, Sturani E, Alberghina L. **Cloning by functional complementation of a mouse cDNA encoding a homologue of CDC25, a Saccharomyces cerevisiae RAS activator.** EMBO J. 1992 Jun;11(6):2151-7. :pmid:`1376246`.

.. [9585556_6] Shou C, Farnsworth CL, Neel BG, Feig LA. **Molecular cloning of cDNAs encoding a guanine-nucleotide-releasing factor for Ras p21.** Nature. 1992 Jul 23;358(6384):351-4. :pmid:`1379346`.

.. [9585556_7] Wei W, Mosteller RD, Sanyal P, Gonzales E, McKinney D, Dasgupta C, Li P, Liu BX, Broek D. **Identification of a mammalian gene structurally and functionally related to the CDC25 gene of Saccharomyces cerevisiae.** Proc Natl Acad Sci U S A. 1992 Aug 1;89(15):7100-4. :pmid:`1379731`.

.. [9585556_8] Fam NP, Fan WT, Wang Z, Zhang LJ, Chen H, Moran MF. **Cloning and characterization of Ras-GRF2, a novel guanine nucleotide exchange factor for Ras.** Mol Cell Biol. 1997 Mar;17(3):1396-406. :pmid:`9032266`.

.. [11438727] Allin C, Ahmadian MR, Wittinghofer A, Gerwert K. **Monitoring the GAP catalyzed H-Ras GTPase reaction at atomic resolution in real time.** Proc Natl Acad Sci U S A. 2001 Jul 3;98(14):7754-9. :pmid:`11438727`.

.. [9247124] Wittinghofer A, Scheffzek K, Ahmadian MR. **The interaction of Ras with GTPase-activating proteins.** FEBS Lett. 1997 Jun 23;410(1):63-7. Review. :pmid:`9247124`.

.. [1569940] Gideon P, John J, Frech M, Lautwein A, Clark R, Scheffler JE, Wittinghofer A. **Mutational and kinetic analyses of the GTPase-activating protein (GAP)-p21 interaction: the C-terminal domain of GAP is not sufficient for full activity.** Mol Cell Biol. 1992 May;12(5):2050-6. :pmid:`1569940`. :download:`PDF </pdf/1569940.pdf>`.

.. [8262937] Eccleston JF, Moore KJ, Morgan L, Skinner RH, Lowe PN. **Kinetics of interaction between normal and proline 12 Ras and the GTPase-activating proteins, p120-GAP and neurofibromin. The significance of the intrinsic GTPase rate in determining the transforming ability of ras.** J Biol Chem. 1993 Dec 25;268(36):27012-9. :pmid:`8262937`.

.. [7548002] Klebe C, Prinz H, Wittinghofer A, Goody RS. **The kinetic mechanism of Ran--nucleotide exchange catalyzed by RCC1.** Biochemistry. 1995 Oct 3;34(39):12543-52. :pmid:`7548002`.

.. [2502546] John J, Schlichting I, Schiltz E, Rösch P, Wittinghofer A.  **C-terminal truncation of p21H preserves crucial kinetic and structural properties.** J Biol Chem. 1989 Aug 5;264(22):13086-92. :pmid:`2502546`. :download:`PDF </pdf/2502546.pdf>`.

.. [9219684] Scheffzek K, Ahmadian MR, Kabsch W, Wiesmüller L, Lautwein A, Schmitz F, Wittinghofer A. **The Ras-RasGAP complex: structural basis for GTPase activation and its loss in oncogenic Ras mutants.** Science. 1997 Jul 18;277(5324):333-8. :pmid:`9219684`. :download:`PDF </pdf/9219684.pdf>`.

