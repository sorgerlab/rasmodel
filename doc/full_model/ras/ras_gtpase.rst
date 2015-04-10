Ras proteins are GTPases
========================

Anatomy of Ras
--------------

For a diagram of the primary structure of Ras, see
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915522/figure/F1/).

.. image:: /images/18568040_ras_anatomy.jpg
    :height: 100px

::

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
                 'oncogenic': ['y', 'n']})

    # Guanine nucleotides
    Monomer('GTP', ['p'])
    Monomer('GDP', ['p'])

Anatomy of Ras regulation:

    [PMID18568040]_: The structural differences between the RasGDP and the
    RasGTP conformations reside mainly in two highly dynamic regions, termed
    switch i (residues 30–40) and switch ii (residues 60–76). Both regions are
    required for the interactions of Ras with upstream as well as downstream
    partners (see also FIG. 2a). The binding of GTP alters the conformation of
    switch i, primarily through the inward reorientation of the side chain of
    Thr35, thereby enabling its interactions with the GTP γ-phosphate as well
    as the Mg2+ ion. Similarly, the γ-phosphate induces significant changes in
    the orientation of the switch ii region through interactions it establishes
    with Gly60 (FIG. 4b).

A side note on alternative splice variants:

    [PMID3304147]_: The three functional ras genes code for highly related proteins
    generically known as p21 (71). The p21 coding sequences of each of these
    genes are equally distributed in four exons except for the K-ras-2 gene,
    which possesses two alternative fourth coding exons (exons IVA and IVB)
    that allow the synthesis of two isomorphic p21 proteins of 188 and 189
    residues that differ in their carboxy terminal domains (54, 60, 72).

Ras binds GTP and GDP
---------------------

Mechanism and rates associated with GTP/GDP binding by Ras.

    [PMID3304147]_: ras proteins, independently of their phylogenetic origin,
    have been shown to bind guanine nucleotides (GTP and GDP)
    ([PMID3304147_22]_ [PMID3304147_23]_ [PMID3304147_24]_ [PMID3304147_25]_)
    and possess intrinsic GTPase activity ([PMID3304147_25]_ [PMID3304147_26]_
    [PMID6147754]_ [PMID6148703]_ [PMID3304147_29]_)

The following statements were taken from a kinetic analysis of Ras and
nucleotide interactions. All rates were measured at 20C.

    [PMID9585556]_: the intrinsic dissociation rate of Ras for GTP (1 × 10-5
    s-1) is 2-fold lower than that for GDP (2 × 10-5 s-1)...

    [PMID9585556]_: Numerically, it was more convenient to use the
    corresponding differential equations with the program FACSIMILE and to
    calculate for 1000 s with the assumption of fast association rate constants
    (in all cases: 10^7 M-1 s-1).

    [PMID9585556]_: The equilibrium dissociation constant for Ras-3′mdGDP (KD1)
    had been determined independently as 9 pM from nucleotide association and
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
-----------------------

GTP hydrolysis by wild-type Ras is very slow in the absence of RasGAPs.

    [PMID1569940]_: It has been reported that the in vitro GTPase activity of
    wild-type p21, which proceeds at a rate of 0.028 min^-1 at 37°C
    ([PMID2502546]_), is accelerated 100- to 200-fold by GAP, as measured under
    nonsaturating conditions.

::

    # Convert 2.8e-2 min^-1 to units of s^-1
    wt_ras_hydrolysis_rate = 2.8e-2 * 60

    def ras_converts_gtp_to_gdp(ras, kcat):
        k = Parameter('k_{0}_gtpase'.format(ras.name), 1.)
        Rule('{0}_converts_GTP_GDP'.format(ras.name),
             ras(gtp=1) % GTP(p=1) >>
             ras(gtp=1) % GDP(p=1),
             k)

    ras_converts_gtp_to_gdp(HRAS, wt_ras_hydrolysis_rate)

Oncogenic Ras mutants have reduced GTP binding and GTPase activity
-------------------------------------------------------------------

[PMID18568040]_: In 1984, three groups reported that mutated Ras oncoproteins
differ functionally from their normal counterparts [PMID6147754]_
[PMID18568040_42]_ [PMID6148703]_. The oncogenic forms of Ras exhibited impaired
GTPase activity, which suggested that the hydrolysis of GTP somehow terminates
the activated state of the protein, which is consistent with the presumed
analogy to the behaviour of G proteins...Furthermore, the link between the
much-studied Gly-to-Val substitution of residue 12 of H-Ras and GTP hydrolysis
was made the following year by Frank McCormick’s group, which noted that
antibodies that are specific to that region blocked GTP binding [PMID18568040_44]_.

[PMID3304147]_: Early studies have predicted that replacement of Gly12 by any other
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

[PMID18568040]_: Other oncogenic mutations (such as Gln61leu in H-Ras) were
also shown to impair GTP hydrolysis [PMID18568040_45]_ and other oncogenic forms of
Ras were later determined to be impaired in GTP hydrolysis (for example, REF.
[PMID18568040_46]_).

[PMID3304147]_: Substitution of Gln61 by 17 different amino acid residues
invariably results in decreased GTPase activity ([PMID3304147_25]_, 117).

.. _FIG4a: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915522/figure/F4/
.. _FIG4b: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915522/figure/F4/

[PMID18568040]_: The overall Ras structure was shown to consist of a
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
    # from mutations in Q61 appear to be attributed to an affect on the
    # catalytic rate.

    # As an implementation detail, note that the mutant rate should be
    # constrained to be less than the wild type rate through the use of an
    # Expression incorporating a scaling parameter between [0, 1].

Autophosphorylation of Ras A59T
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[PMID3304147]_: In addition to GTP/GDP binding and GTPase activity, ras proteins
carrying an Ala59 -> Thr59 mutation exhibit an autophosphorylating activity of
an, as yet, unknown biological significance [PMID3304147_23]_. In all cases, Thr59
has been found to be the phosphate receptor site (106). No transphosphorylating
activity has been detected with any ras protein, including those carrying Thr59
mutations::

    # Add autophosphorylation of Ras A59T if it later turns out to be significant.


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

