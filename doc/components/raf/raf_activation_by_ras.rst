Raf activation by Ras
=====================

In its basal state, RAF is present in a 'closed' conformation, wherein the N terminus of the RAF protein interacts with, and inhibits the C terminus. Ligand binding to receptor tyrosine kinase (RTK) results in activation of the RTK, leading to RAS actviation (RAS-GTP). It is not completely understood RAS activates RAF. GTP bound RAS binds to the RAS binding domain (RBD) of  RAF, and prevetns the inhibitory interaction between the C and N termini of RAF. This results in RAF assuming an 'open' conformation. A further step in the activation of RAF is dimerization, which stabilizes its 'open' conformation. 28% of Melanomas have RAS mutations, of which 94% are mutations in NRAS. Hence, we instantiate the model for the vemurfenib resistance scenario using NRAS.

::
   
    # BRAF dimerization
    Rule('BRAF_dimerization',
         BRAF(d=None, ras=None) + BRAF(d=None, ras=None, vem=None) <>
         BRAF(d=1, ras=None) % BRAF(d=1, ras=None, vem=None), kaf, kar)

    # KRAS binding BRAF monomers
    Rule('KRAS_binding_BRAF_monomers',
         BRAF(ras=None, d=None) + KRAS(raf=None, state='gtp') <>
         BRAF(ras=1, d=None) % KRAS(raf=1, state='gtp'), kdf, kdr)

    # KRAS binding BRAF dimers
    Rule('KRAS_binding_BRAF_dimers',
         BRAF(ras=None, d=1) % BRAF(ras=None, d=1) +
         KRAS(raf=None, state='gtp') + KRAS(raf=None, state='gtp') <>
         BRAF(ras=2, d=1) % BRAF(ras=3, d=1) %
         KRAS(raf=2, state='gtp') % KRAS(raf=3, state='gtp'), kbf, kbr)

    # KRAS:BRAF dimerization
    Rule('KRASBRAF_dimerization',
         BRAF(d=None, ras=ANY) + BRAF(d=None, ras=ANY, vem=None) <>
         BRAF(d=1, ras=ANY) % BRAF(d=1, ras=ANY, vem=None), kcf, kcr)
	 
    # Release KRAS:GDP from BRAF
    Rule('KRAS_GDP_dissoc_BRAF',
         KRAS(state='gdp', raf=1) % BRAF(ras=1) >>
         KRAS(state='gdp', raf=None) + BRAF(ras=None), koff)
	 

Vemurafenib inhibits RAF
========================
Vemurfenib is considered to be a BRAF selective inhibitor. However, in vitro experiments show that Vemurafenib can bind BRAF V600E, BRAF WT, CRAF WT, and ARAF with similar affinity. Moreover, the IC50 values computed for each isoforms in vitro are comparable. In contrast, in vivo experiments indicate that Vemurafenib is effective only against BRAF V600E mutants, and cause paradoxial activation in WT BRAF and WT CRAF cell lines. hence, there is a disconnect between conclusions from in vitro and in vivo observations. The mechanims for paradoixal activation are obscure. Answering these questions require a molecular-level analysis and undestanding of interaction beteen RAF isoforms and Vemurafenib.

Recent experiments have shown that Vemurafenib is iniffective in binding RAF dimers. ~30 fold higer concetration of Vemurafenib is required to inhibit constitiutuve BRAF dimers. Based on this knowledge the authors propose the following mechanism for the development of BRAF resistanse.

1. BRAFV600E mutants increase level of phosphorylated ERK.
2. Higher level of ERK translates to inhibition of RTK and SOS
3. Therefore, active RAS is present at low level in BRAFV600E mutants, and mutant BRAF can continue signaling independent of RAS.
4. Addition of Vemurafenib inhibits BRAFV600E, lowers ERKlevel and relieves feedback on RTK and SOS.
5. RAS is activated, which induces dimerization of BRAF.
6. Vemurafenib binds one monomer in the BRAF dimer effecitvely, but binds the second protomer in the dimer with significantly lower affinity.
7. As a result, BRAF dimers can activate ERK via the protomeer that is not bound to Vemurafenib, resulting in the partial restoration of ERK phosphorylation (imperfect adaptation)


::

    # BRAF:Vem dimerization to give 2(BRAF:Vem) g = a * f
    Rule('BRAF_Vem_dimerization',
         BRAF(d=None, ras=None, vem=ANY) + BRAF(d=None, ras=None, vem=ANY) <>
         BRAF(d=1, ras=None, vem=ANY) % BRAF(d=1, ras=None, vem=ANY), kgf, kgr)

    # KRAS:BRAF:Vem dimerization to give 2( KRAS:BRAF:Vem) h = c * a
    Rule('KRAS_BRAF_Vem_dimerization',
         BRAF(d=None, ras=ANY, vem=ANY) + BRAF(d=None, ras=ANY, vem=ANY) <>
         BRAF(d=1, ras=ANY, vem=ANY) % BRAF(d=1, ras=ANY, vem=ANY), khf, khr)

    # 1st Vemurafenib binds
    Rule('First_binding_Vemurafenib',
         BRAF(vem=None) % BRAF(vem=None) + Vem(raf=None) <>
         BRAF(vem=1) % BRAF(vem=None) % Vem(raf=1), kef, ker)

    # 2nd Vemurafenib binding
    Rule('Second_binding_vemurafenib',
         BRAF(vem=None) % BRAF(vem=ANY) + Vem(raf=None) <>
         BRAF(vem=1) % BRAF(vem=ANY) % Vem(raf=1), kff, kfr)

    # Vemurafenib binds BRAF monomer
    Rule('Vemurafenib_binds_BRAF_monomer',
         BRAF(vem=None, d=None) + Vem(raf=None) <>
         BRAF(vem=1, d=None) % Vem(raf=1), kef, ker)


References
----------

.. [PMID11237210] Avruch J, Khokhlatchev A, Kyriakis JM, Luo Z, Tzivion G, Vavvas D, Zhang XF.  **Ras activation of the Raf kinase: tyrosine kinase recruitment of the MAP kinase cascade.** Recent Prog Horm Res. 2001;56:127-55. :pmid:`11237210`. :download:`PDF </pdf/11237210.pdf>`

.. [PMID21862573] Hibino K, Shibata T, Yanagida T, Sako Y. **Activation kinetics of RAF protein in the ternary complex of RAF, RAS-GTP, and kinase on the plasma membrane of living cells: single-molecule imaging analysis.** J Biol Chem. 2011 Oct 21;286(42):36460-8. :doi:`10.1074/jbc.M111.262675.` :pmid:`21862573` :download:`PDF </pdf/21862573.pdf>`

.. [PMID11447113] Chong H, Lee J, Guan KL. **Positive and negative regulation of Raf kinase activity and function by phosphorylation.** EMBO J. 2001 Jul 16;20(14):3716-27. :pmid:`11447113` :download:`PDF </pdf/11447113.pdf>`

.. [PMID15664184] Dumaz N, Marais R. **Raf phosphorylation: one step forward and two steps back.** Mol Cell. 2005 Jan 21;17(2):164-6. :pmid:`15664184` :download:`PDF </pdf/15664184.pdf>`

.. [PMID15664191] Dougherty MK1, Müller J, Ritt DA, Zhou M, Zhou XZ, Copeland TD, Conrads TP, Veenstra TD, Lu KP, Morrison DK. **Regulation of Raf-1 by direct feedback phosphorylation.** Mol Cell. 2005 Jan 21;17(2):215-24. :pmid:`15664191` :download:`PDF </pdf/15664191.pdf>`

.. [lavoie] Lavoie H, Therrien M. **Regulation of RAF protein kinases in ERK signalling.** :doi:`10.1038/nrm3979` :download:`PDF </pdf/lavoie.pdf>`
	    
.. [2634358] Yao Z, Torres NM, Tao A, Gao Y, Luo L, Li Q, Stanchina E, Abdel-Wahab O, Solit DB, Poulikakos PI, Rosen N. **BRAF mutants evade ERK-dependent feedback by different mechanisms that determine their sensitivity to pharmacological inhibition.** Cancer Cell. 2015 Sept 1 15;28:270-83. :pmid:`26343582`

.. [2420239] Lito P, Rosen N, Solit DB. ** Tumor adaptation and resistance to RAF inhibitors.** Nature MEdicine. 2013 Nov; 19(11):1401-9. :pmid:`24202393`

.. [23153539] Lito P, Pratilas CA, Joseph EW, Tadi M, Halilovic E, Zubrowski M, Huan A, Wong WL, Callahan MK, Merghoun T, Wolchok JD, Stanchina E, Chandrarlapaty S, Paulikakos PI, Fagin JA, Rosen N, **Relief of profound feedback inhibition of mitogenic signaling by RAF inhibitors attenuates their activity in BRAFV600E melanomas.** Cancer Cell, 2012 Nov 12;22:668-82. pmid:`23153539`

.. [21107323] Nazarian R, Shi H, Wanf Q, Kon X, Koya RC, Lee H, Chen Z, Lee M-K, Attar N, Sazegar H, chodon T, Nelson SF, McArthur G, Sosman JA, Ribas A, Lo RS. ** Melanomas acquire resistance to B-RAF(V600E) inhibition by RTK or N-RAS upregulation.** Nature. 2010 DEc 16;468: 973-7. :pmid:`2110732`	      
