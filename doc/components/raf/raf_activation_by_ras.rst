Raf activation by Ras
=====================

In its basal state, RAF is present in a 'closed' conformation, wherein the N terminus of the RAF protein interacts with, and inhibits the C terminus. Ligand binding to receptor tyrosine kinase (RTK) results in activation of the RTK, leading to RAS actviation (RAS-GTP). It is not completely understood RAS activates RAF. GTP bound RAS binds to the RAS binding domain (RBD) of  RAF, and prevetns the inhibitory interaction between the C and N termini of RAF. This results in RAF assuming an 'open' conformation. A further step in the activation of RAF is dimerization, which stabilizes its 'open' conformation. 28% of Melanomas have RAS mutations, of which 94% are mutations in NRAS. Hence, we instantiate the model for the vemurfenib resistance scenario using NRAS.

::
   from pysb.macros import bind, _macro_rule
   
   def ras_activates_raf(ras, raf, klist):  
    
        kaf, kar, kdf, kdr, kbf, kbr, kcf, kcr, kf5, koff = klist

	# RAF dimerization
	bind(raf(ras=None), 'd', raf(ras=None, vem=None), 'd', [kaf, kar])
	
	# RAS binding RAF monomers
	bind(raf(d=None), 'ras', ras(gtp=ANY), 'raf', [kdf, kdr])

	# RAS binding RAF dimers
	Rule('RAS_binding_RAF_dimers',
	     raf(ras=None, d=1) % raf(ras=None, d=1) +
	     ras(raf=None, gtp=ANY) + ras(raf=None, gtp=ANY) <>
	     raf(ras=2, d=1) % ras(ras=3, d=1) %
	     ras(raf=2, gtp=ANY) % ras(raf=3, gtp=ANY), kbf, kbr)

	# RAS:RAF dimerization
	bind(raf(ras=ANY), 'd', raf(ras=ANY, vem=None), 'd', [kcf, kcr])
		
	# KRAS deactivates itself
	# Making this step reversible increased combinatorial complexity manifold     
	Rule('KRAS_inactivation',
             ras(gtp=1) % GTP(ras=1) >>
             ras(gtp=1) % GDP(ras=1),
             kf5)

	# Release RAS:GDP from RAF
	Rule('RAS_GDP_dissoc_RAF',
	     GDP(ras=2) % ras(gtp=2, raf=1) % raf(ras=1) >>
	     GDP(ras=2) % ras(gtp=2, raf=None) + raf(ras=None), koff)

	     
    # def kras_activates_braf(model):
    #     KRAS = model.components['KRAS']
    #     BRAF = model.components['BRAF']
    #     ras_activates_raf(KRAS, BRAF, ...)
    

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

    def vemurafenib_binds_raf(model, raf, klist):

	Vemurafenib = model.monomers["Vemurafenib"]
        kgf, kgr, khf, khr, kef, ker, kff, kfr = klist
	
	# RAF:Vem dimerization to give 2(RAF:Vem) g = a * f
	bind(raf(ras=None, vem=ANY), 'd', raf(ras=None, vem=ANY), 'd', [kgf, kfr], ('kf', 'kr')]
	
	# RAS:RAF:Vem dimerization to give 2(RAS:RAF:Vem) h = c * a
        bind(raf(ras=ANY, vem=ANY), 'd', raf(ras=ANY, vem=ANY) 'd', [khf, khr], ('kf', 'kr'))
	
	# 1st Vemurafenib binds
	_macro_rule('First_binding_Vemurafenib',
	     raf(vem=None) % raf(vem=None) + Vem(raf=None) <>
	     raf(vem=1) % raf(vem=None) % Vem(raf=1), [kef, ker], ('kf', 'kr'))

	# 2nd Vemurafenib binding
	_macro_rule('Second_binding_vemurafenib',
	     raf(vem=None) % raf(vem=ANY) + Vem(raf=None) <>
	     raf(vem=1) % raf(vem=ANY) % Vem(raf=1), [kff, kfr], ('kf', 'kr'))

	# Vemurafenib binds RAF monomer
	bind(raf(d=None), 'vem', Vem(), 'raf', [kef, ker], ('kf', 'kr')


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
