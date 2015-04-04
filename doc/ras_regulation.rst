Ras activity and regulation
===========================

Anatomy of Ras
--------------

For a diagram of the primary structure of Ras, see
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915522/figure/F1/).

[3304147]_: The three functional ras genes code for highly related proteins
generically known as p21 (71). The p21 coding sequences of each of these genes
are equally distributed in four exons except for the K-ras-2 gene, which
possesses two alternative fourth coding exons (exons IVA and IVB) that allow
the synthesis of two isomorphic p21 proteins of 188 and 189 residues that
differ in their carboxy terminal domains (54, 60, 72)::

    from pysb.macros import bind_table

    for ras_name in ['KRAS', 'NRAS', 'HRAS']:
        Monomer(ras_name,
                ['gtp', 'p_loop', 'switch1', 'switch2', 'CAAX', 'oncogenic'],
                {'oncogenic': ['y', 'n']})
    Monomer('GTP', ['p'])
    Monomer('GDP', ['p'])

Ras proteins are GTPases
------------------------

[3304147]_: ras proteins, independently of their phylogenetic origin, have been
shown to bind guanine nucleotides (GTP and GDP) ([3304147_22]_ [3304147_23]_
[3304147_24]_ [3304147_25]_) and possess intrinsic GTPase activity
([3304147_25]_ [3304147_26]_ [3304147_27]_ [3304147_28]_ [3304147_29]_)::

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

Oncogenic Ras mutants have reduced GTP binding and GTPase activity
-------------------------------------------------------------------

[18568040]_: In 1984, three groups reported that mutated Ras oncoproteins
differ functionally from their normal counterparts [18568040_41]_
[18568040_42]_ [18568040_43]_. The oncogenic forms of Ras exhibited impaired
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
state.

.. _FIG4a: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915522/figure/F4/
.. _FIG4b: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915522/figure/F4/

::

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
-------------------------------

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

Post-translational modifications
--------------------------------

An initial study in this area, published in 1982, showed that the mature form
of viral H-Ras localized to the cell membrane47. Several months later it was
demonstrated that viral H-Ras is palmitoylated at the C terminus; the resulting
attached lipid moiety facilitated its association with the membrane48. The
functional connection between this lipid modification and Ras function was made
by Douglas Lowy’s group in 1984, which showed that lipid binding and membrane
association were actually required for the transforming activity of the viral
H-Ras oncoprotein49,50.

working with cellular H-Ras, Stuart Aaronson’s group proceeded to demonstrate that this C-terminal processing and membrane recruitment of Ras is a prerequisite to its biochemical activation51.

The molecular mechanisms of Ras lipid processing were laid out over the subsequent 5 years through a series of observations using yeast genetics, protein biochemistry and in vitro cellular systems52–57 (FIGS 2,3).3).

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

.. [3304147] Barbacid M. ras genes. Annu Rev Biochem. 1987;56:779-827. Review. :pmid:`3304147`.

.. [3304147_22] Scolnick EM, Papageorge AG, Shih TY. Guanine nucleotide-binding activity as an assay for src protein of rat-derived murine sarcoma viruses. Proceedings of the National Academy of Sciences of the United States of America. 1979;76(10):5355-5359. :pmid:`228288`.

.. [3304147_23] Shih TY, Papageorge AG, Stokes PE, Weeks MO, Scolnick EM. Guanine nucleotide-binding and autophosphorylating activities associated with the p21src protein of Harvey murine sarcoma virus. Nature. 1980 Oct 23;287(5784):686-91. :pmid:`6253810`.

.. [3304147_24] Tamanoi, F., Walsh, M., Kataoka, T., & Wigler, M. (1984). A product of yeast RAS2 gene is a guanine nucleotide binding protein. Proceedings of the National Academy of Sciences of the United States of America, 81(22), 6924–6928. :pmid:`6438624`.

.. [3304147_25] Temeles GL, Gibbs JB, D'Alonzo JS, Sigal IS, Scolnick EM. Yeast and mammalian ras proteins have conserved biochemical properties. Nature. 1985 Feb 21-27;313(6004):700-3. :pmid:`3919305`.

.. [3304147_26] Gibbs JB, Sigal IS, Poe M, Scolnick EM. Intrinsic GTPase activity distinguishes normal and oncogenic ras p21 molecules. Proc Natl Acad Sci U S A. 1984 Sep;81(18):5704-8. :pmid:`6148751`.

.. [3304147_27] McGrath JP, Capon DJ, Goeddel DV, Levinson AD. Comparative biochemical properties of normal and activated human ras p21 protein. Nature. 1984 Aug 23-29;310(5979):644-9. :pmid:`6147754`.

.. [3304147_28] Sweet RW, Yokoyama S, Kamata T, Feramisco JR, Rosenberg M, Gross M. The product of ras is a GTPase and the T24 oncogenic mutant is deficient in this activity. Nature. 1984 Sep 20-26;311(5983):273-5. :pmid:`6148703`.

.. [3304147_29] Manne V, Bekesi E, Kung HF. Ha-ras proteins exhibit GTPase activity: point mutations that activate Ha-ras gene products result in decreased GTPase activity. Proc Natl Acad Sci U S A. 1985 Jan;82(2):376-80. :pmid:`2982154`.

.. [18568040] Karnoub AE, Weinberg RA. Ras oncogenes: split personalities. Nature reviews Molecular cell biology. 2008;9(7):517-531. doi:10.1038/nrm2438.

.. [18568040_41] McGrath JP, Capon DJ, Goeddel DV, Levinson AD. Comparative biochemical properties of normal and activated human ras p21 protein. Nature. 1984;310:644–649. :pmid:`6147754`.

.. [18568040_42] Gibbs JB, Sigal IS, Poe M, Scolnick EM. Intrinsic GTPase activity distinguishes normal and oncogenic ras p21 molecules. Proc Natl Acad Sci USA. 1984;81:5704–5708. :pmid:`6148751`.

.. [18568040_43] Sweet RW, et al. The product of ras is a GTPase and the T24 oncogenic mutant is deficient in this activity. Nature. 1984;311:273–275. :pmid:`6148703`.

.. [18568040_44] Clark R, Wong G, Arnheim N, Nitecki D, McCormick F. Antibodies specific for amino acid 12 of the ras oncogene product inhibit GTP binding. Proc Natl Acad Sci USA. 1985;82:5280–5284.:pmid:`3927300`.

.. [18568040_45] Der CJ, Finkel T, Cooper GM. Biological and biochemical properties of human rasH genes mutated at codon 61. Cell. 1986;44:167–176. :pmid:`3510078`.

.. [18568040_46] Trahey M, McCormick F. A cytoplasmic protein stimulates normal N-ras p21 GTPase, but does not affect oncogenic mutants. Science.  1987;238:542–545. References 41–46 established that oncogenic mutation of ras affects its nucleotide cycle. :pmid:`2821624`.
