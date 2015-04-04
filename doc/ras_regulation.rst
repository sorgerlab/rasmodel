Karnoub and Weinberg: Ras oncogenes: split personalities
========================================================

* Title: Ras oncogenes: split personalities
* PMID: 18568040
* Link: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915522/

Activated mutants and the GTPase cycle
--------------------------------------

Structure of Ras (see
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915522/figure/F1/)::

    for ras_name in ['KRAS', 'NRAS', 'HRAS']:
        Monomer(ras_name,
                ['gtp', 'p_loop', 'switch1', 'switch2', 'CAAX', 'oncogenic'],
                {'oncogenic': ['y', 'n']})
    Monomer('GTP', ['p'])
    Monomer('GDP', ['p'])

:pmid:`18568040`: In 1984, three groups reported that mutated Ras oncoproteins
differ functionally from their normal counterparts [18568040_41]_
[18568040_42]_ [18568040_43]_. The oncogenic forms of Ras exhibited impaired
GTPase activity, which suggested that the hydrolysis of GTP somehow terminates
the activated state of the protein, which is consistent with the presumed
analogy to the behaviour of G proteins...Furthermore, the link between the
much-studied Gly-to-Val substitution of residue 12 of H-Ras and GTP hydrolysis
was made the following year by Frank McCormick’s group, which noted that
antibodies that are specific to that region blocked GTP binding [18568040_44]_::

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

:pmid:`18568040`: Other oncogenic mutations (such as Gln61leu in H-Ras) were
also shown to impair GTP hydrolysis [18568040_45]_ and other oncogenic forms of
Ras were later determined to be impaired in GTP hydrolysis (for example, REF.
[18568040_46]_).

:pmid:`18568040`: The overall Ras structure was shown to consist of a
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

    # Note that the mutant rate should be constrained to be less than
    # the wild type rate through the use of an Expression incorporating a
    # scaling parameter between [0, 1].
    Parameter('k_mut_gtpase', 0.1)

    # Mutant Ras has diminished GTPase activity:
    for ras in [KRAS, HRAS, NRAS]:
        ras_mut = ras(oncogenic='y')

        Rule('{0.name}_mut_converts_GTP_GDP'.format(ras),
             ras_mut(gtp=1) % GTP(p=1) >>
             ras_mut(gtp=1) % GDP(p=1),
             k_mut_gtpase)

However, the extent of such impairment did not always correlate with
transformation, suggesting that compromised intrinsic GTP hydrolysis was
necessary but not sufficient for aberrant Ras activation.

STL::

    # How to indicate necessity vs. sufficiency?
Anatomy of Ras regulation
-------------------------

:pmid:`18568040`: The structural differences between the RasGDP and the RasGTP
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

.. [18568040_41] McGrath JP, Capon DJ, Goeddel DV, Levinson AD. Comparative biochemical properties of normal and activated human ras p21 protein. Nature. 1984;310:644–649. :pmid:`6147754`.

.. [18568040_42] Gibbs JB, Sigal IS, Poe M, Scolnick EM. Intrinsic GTPase activity distinguishes normal and oncogenic ras p21 molecules. Proc Natl Acad Sci USA. 1984;81:5704–5708. :pmid:`6148751`.

.. [18568040_43] Sweet RW, et al. The product of ras is a GTPase and the T24 oncogenic mutant is deficient in this activity. Nature. 1984;311:273–275. :pmid:`6148703`.

.. [18568040_44] Clark R, Wong G, Arnheim N, Nitecki D, McCormick F. Antibodies specific for amino acid 12 of the ras oncogene product inhibit GTP binding. Proc Natl Acad Sci USA. 1985;82:5280–5284.:pmid:`3927300`.

.. [18568040_45] Der CJ, Finkel T, Cooper GM. Biological and biochemical properties of human rasH genes mutated at codon 61. Cell. 1986;44:167–176. :pmid:`3510078`.

.. [18568040_46] Trahey M, McCormick F. A cytoplasmic protein stimulates normal N-ras p21 GTPase, but does not affect oncogenic mutants. Science.  1987;238:542–545. References 41–46 established that oncogenic mutation of ras affects its nucleotide cycle. :pmid:`2821624`.
