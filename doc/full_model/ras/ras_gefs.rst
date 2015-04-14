.. _ras_gefs:

Ras regulation by GEFs
======================

Identities of Ras-specific GEFs
-------------------------------

There are several guanosine exchange factors (GEFs) that are specific to Ras
subfamily proteins.

    [PMID9585556]_: Several genes have been isolated from different organisms
    encoding proteins that have a GEF activity specific for Ras (for which we
    use the general name RasGEFs throughout this paper): SOS1 and SOS2
    ([PMID9585556_1]_ [PMID9585556_2]_ [PMID9585556_3]_ [PMID9585556_4]_) ;
    Cdc25Mm, also called RasGrf ([PMID9585556_5]_ [PMID9585556_6]_
    [PMID9585556_7]_); and mRas-GRF2 ([PMID9585556_8]_).

RasGEFs contain a specific domain responsible for activating Ras proteins.

    [PMID9585556]_: The RasGEFs are proteins of considerable length, 120 - 160
    kDa, and contain several regions which are generally accepted to represent
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
-----------------

Some key features of the mechanism:

1. Ras binds GEFs in the absence of nucleotides
2. GEF binding causes Ras to release GTP/GDP.
3. Rebinding of nucleotides causes Ras to release the GEF.
4. Reloading of Ras with GTP vs. GDP is not determined by the GEF, but rather
   by the relative cellular concentrations of the nucleotides.

    [PMID9690470]_: Biochemical studies of Ras exchange factors have shown that
    the complex of Ras with these proteins is stable in the absence of
    nucleotides and is dissociated by the rebinding of either GDP or GTP
    ([PMID9585556]_ [PMID9690470_17]_ [PMID9690470_18]_ [PMID9690470_21]_
    [PMID9690470_22]_). The principal role for the exchange factor is to
    facilitate nucleotide release, and it does not seem to control
    significantly the preferential rebinding of GTP over GDP ([PMID9585556]_,
    [PMID9690470_22]_, [PMID9690470_23]_).  Cellular concentrations of GTP are
    10-fold higher than GDP, which results in the loading of GTP onto Ras.

The fact that GTP/GDP can displace GEFs, while GEFs can also displace GTP/GDP,
leads to a paradox that is resolved by the fact that Ras undergoes a
conformational change that retains the necessary "state". The structural basis
of this conformational change is described as follows:

    [PMID9690470]_: As a nucleotide-exchange factor, Sos functions under two
    apparently conflicting imperatives. The interaction between Sos and Ras
    must be strong enough to dislodge the tightly bound nucleotide, but the Ras
    – Sos complex must also be poised for subsequent displacement by incoming
    nucleotides. The structure of the Ras – Sos complex shows that Ras and Sos
    meet these demands by forming a tight complex that is anchored at one end
    of the nucleotide- binding site, where phosphate and magnesium are normally
    bound. The interface between Sos and Ras is mainly hydrophilic, suggesting
    a ready unzippering through water-mediated displacements of the
    coordinating side chains. The main interacting elements of Sos avoid direct
    occlusion of the nucleotide-binding site, except the region where the
    terminal phosphate groups and the magnesium ion are bound. This feature
    allows incoming nucleotides to reverse the process by competing for the
    groups that ligate the phosphate and metal ion.

This conformational state change has been analyzed kinetically and identified
as the process of the nucleotide being first loosely and then tightly bound,
(see also :ref:`ras_gtpase`):

    [PMID9690470]_: Kinetic analysis of nucleotide association shows that the
    reaction proceeds by the formation of a ternary complex of a loosely bound
    nucleotide and Ras – Cdc25Mm followed by conversion to a form in which the
    nucleotide is tightly bound to Ras [PMID9585556]_. In light of the
    structure of the Ras–Sos complex, the first step can be interpreted as the
    interaction of the base and the ribose of the nucleotide with the part of
    the Ras binding site that is not occluded by Sos. The second step would
    involve a conformational change in the Switch 2 segment and release of
    Switch 1, resulting in the restructuring of a competent binding site for
    phosphate and magnesium, and the subsequent dissociation of Sos.

The kinetic analysis described in [PMID9585556]_ resulted in the following reaction scheme for the interactions between Ras, GTP/GDP, and GEFs:

.. image:: /images/9585556_rasgef_cycle.png
    :width: 600px

Note that the upper equilibria for Ras-nucleotide binding, K1a and K1b, were
implemented in the section :ref:`ras_gtpase`, along with corresponding rates.
Here we implement only the equilibria involving GEFs: K2, K3, K4a and K4b.

::

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

Rates of GEF activation
-----------------------

::

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

# = kf1a, i.e., on rate is insensitive to presence of GEF

::

    # Isomerization of Ras-RasGEF-GXP from loose to tight
    kf4b = 20.4         # s^-1
    kr4b = 3.9          # s^-1



The following study used purified HRAS and mouse RASGRF1:

    [PMID9690470]_: The mechanism of nucleotide release by the catalytic domain
    of murine Cdc25 (Cdc25Mm) has been investigated recently using
    fluorescently labelled nucleotides [PMID9585556]_.  The affinity of Cdc25Mm
    for nucleotide-free Ras (Kd = 4.6 nM) is found to be several orders of
    magnitude higher than that for nucleotide-bound Ras, and the maximal
    acceleration by Cdc25Mm of the rate of dissociation of nucleotide is more
    than 10^5.

    [PMID9585556]_: The best fit of our data resulted in similar quantum yields
    and a value of 4.6 nM for KD2 (NOTE: Kd between nucleotide-free H-Ras and
    RasGRF1). A variation in the value for KD2 of approximately 2-fold resulted
    in fits of comparable quality.

The activity of GEF (RASGRF1 in this case) does not depend on whether Ras
(HRAS) is loaded with GTP or GDP.

    [PMID9585556]_: However, since the intrinsic dissociation rate of Ras for
    GTP (1 × 10-5 s-1) is 2-fold lower than that for GDP (2 × 10-5 s-1), the
    stimulatory action of Cdc25Mm285 is practically independent of the nature
    of the bound nucleotide.

    [PMID9585556]_: Although we did not reach complete saturation at 600 μM
    Ras‚nucleotide, the data could be fitted to obtain a maximal rate of
    3′mdGDP release from Ras of 3.9 s-1 and an apparent Km value of 386 μM.
    Since the intrinsic dissociation rate of 3′mdGDP is 2 × 10-5 s-1 (Table 1),
    the acceleration of GDP dissociation from Ras by this GEF is approximately
    2 × 10^5-fold. An apparent Km of approximately 300 μM was obtained
    for the triphosphate-bound form of Ras, confirming that there is no
    pronounced specificity toward the nature of the Ras-bound nucleotide (data
    not shown).

.. warning:: GEF binding to GTP bound Ras?

    Can GEFs bind to Ras and cause ejection of nucleotide before the GTP/GDP
    conversion is complete? Moreover, if GEF binds to Ras-GTP, can the
    hydrolysis to GDP proceed while GEF is bound?

Instantiate the RasGEF cycle for HRAS and RASGRF1::

    #ras_gef_exchange_cycle(HRAS, RASGRF1, GTP, GDP)

[PMID9585556]_: Therefore, we tested the nucleotide specificity of the
interaction of Cdc25Mm285 (CdcMm285 is the fragment of CdcMm/RasGRF1 containing
the RasGEF domain) with Ras. Figure 1 shows the release of Ras-bound 3′mdGDP or
3′mdGTP (4 μM), in the presence of an excess of unlabeled nucleotide and in the
presence or absence of 1 μM Cdc25Mm285. The Cdc25Mm285-stimulated dissociation
rate of Ras-3′mdGDP is approximately twice that of Ras-3′mdGTP, with values of
0.0098 and 0.0046 s-1, respectively.  However, since the intrinsic dissociation
rate of Ras for GTP (1 × 10-5 s-1) is 2-fold lower than that for GDP (2 × 10-5
s-1), the stimulatory action of Cdc25Mm285 is practically independent of the
nature of the bound nucleotide. The difference in stimulated dissociation rates
is somewhat smaller than the results of Jacquet et al. (16) but is similar to
the results with the yeast proteins CDC25 and RAS2 obtained by Haney and Broach
(28).

[PMID9690470]_: The overall shape of the catalytic domain of Sos is that of an
oblong bowl (Fig. 2), with Ras bound at the centre of the bowl. The regions of
Ras that interact most closely with Sos include the phosphate-binding P-loop
(residues 10 – 17) and surrounding segments (including strand 􏰧1 and helix 􏰦1),
the Switch 1 region (defined here as residues 25–40) and the Switch 2 region
(defined here as residues 57 – 75). Additional interactions are seen with helix
3 (residues 95–105; Fig. 3a, b). The interface between Ras and Sos is primarily
hydrophilic and very extensive, with 3,600 A^2 of surface area buried in the
complex.

[PMID9690470]_: The most obvious effect of Sos binding to Ras is the opening of
the nucleotide binding site as a result of the displacement of Switch 1 of Ras
by the insertion of the helical hairpin formed by aH and aI of Sos (Fig. 5)

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
---------------------------------------

[PMID9585556]_: Three mammalian isoforms of Ras, H-, K-, and N-Ras, have been
identified which are highly conserved intheirprimarysequence.
Thesignificanceofhavingmore than one isoform is not understood at present,
although the isoforms may have different functions in different tissues, since
certain types of tumors have a preference for a particular activated Ras gene,
such as K-Ras for lung, colon and pancreas cancers and N-Ras for myeloid
leukemias (25). To see whether Cdc25Mm285 acts differently on the three
isoforms, we tested the GEF activity of Cdc25Mm285 on these proteins. As
summarized in Table 1, Cdc25Mm285 is active on all isoforms, being somewhat
more active on N-Ras, in accordance with the results of Leonardsen et al. (26).

.. raw:: html

    <script>
        window.setTimeout(function() {
        $('div.highlight-python pre > span.c:last-child').each(
            function () {
                if ($(this).text() == '#') {
                    $(this.nextSibling).detach();
                    $(this).detach();
                }
            }
        );
        }, 1000);
    </script>

