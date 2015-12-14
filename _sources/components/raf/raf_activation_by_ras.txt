Raf activation by Ras
=====================

In its basal state, RAF is present in a 'closed' conformation, wherein the N terminus of the RAF protein interacts with, and inhibits the C terminus. Ligand binding to receptor tyrosine kinase (RTK) results in activation of the RTK, leading to RAS actviation (RAS-GTP). It is not completely understood how RAS activates RAF. GTP bound RAS binds to the RAS binding domain (RBD) of  RAF, and prevents the inhibitory interaction between the C and N termini of RAF. This results in RAF assuming an 'open' conformation. A further step in the activation of RAF is dimerization, which stabilizes its 'open' conformation.

::

   from pysb import Monomer, Rule, Parameter, Annotation, ANY
   from pysb.macros import bind, _macro_rule, catalyze_state
   from pysb.util import alias_model_components


   def RAF_monomers():
       Monomer('RAF', ['ras', 'd', 'vem', 'erk'])
       Monomer('Vem', ['raf'])

       # IC values
       # --------
       Parameter('RAF_0', 1e5)
       Parameter('Vem_0', 1000)

       alias_model_components()

       # Initial conditions
       # ------------------
       Initial(RAF(d=None, ras=None, erk=None, vem=None), RAF_0)
       Initial(Vem(raf=None), Vem_0)


      def RAF_dynamics():
       # Parameters
       # -----------
       Parameter('kaf', 1e-6)
       Parameter('kar', 1)
       Parameter('kbf', 0.5)
       Parameter('kbr', 1e-11)
       Parameter('kcf', 1)
       Parameter('kcr', 0.0001)
       Parameter('kdf', 1)
       Parameter('kdr', 0.1)
       Parameter('kef', 1e-2)
       Parameter('ker', 0.1)
       Parameter('kff', 1e-5)
       Parameter('kfr', 1)
       Parameter('kgf', 1e-11)
       Parameter('kgr', 1)
       Parameter('khf', 1e-2)  # 100)
       Parameter('khr', 1)  # 1)
       Parameter('koff', 1)

       alias_model_components()

       # Rules
       # -----

       # RAF that is not bound to RAS binds RAF
       Rule('RAF_dimerization',
	    RAF(d=None) + RAF(d=None, ras=None) <>
	    RAF(d=1) % RAF(d=1, ras=None), kaf, kar)

       # RAS binds RAF
       Rule('RAS_binding_RAF_monomers',
	    RAF(ras=None) + RAS(raf=None, sos1=None, gtp=ANY) <>
	    RAF(ras=1) % RAS(raf=1, sos1=None, gtp=ANY), kdf, kdr)


       # RAF, bound to RAS, binds RAF bound to RAS
       Rule('RASRAF_dimerization',
	    RAF(d=None, ras=ANY) + RAF(d=None, ras=ANY) <>
	    RAF(d=1, ras=ANY) % RAF(d=1, ras=ANY), kcf, kcr)

      
Vemurafenib inhibits RAF
========================

Vemurfenib is considered to be a BRAF selective inhibitor. However, in vitro experiments show that Vemurafenib can bind BRAF V600E, BRAF WT, CRAF WT, and ARAF with similar affinity. Moreover, the IC50 values computed for each isoforms in vitro are comparable. In contrast, in vivo experiments indicate that Vemurafenib is effective only against BRAF V600E mutants, and cause paradoxial activation in WT BRAF and WT CRAF cell lines. hence, there is a disconnect between conclusions from in vitro and in vivo observations. The mechanims for paradoixal activation are obscure. Answering these questions require a molecular-level analysis and undestanding of interaction beteen RAF isoforms and Vemurafenib.

Recent experiments have shown that Vemurafenib is iniffective in binding RAF dimers. ~30 fold higer concetration of Vemurafenib is required to inhibit constitiutuve BRAF dimers. Based on this knowledge the authors propose the following mechanism for the development of BRAF resistanse.

1. BRAFV600E mutants increase level of phosphorylated ERK.
2. Higher level of ERK translates to inhibition of RTK and SOS
3. Therefore, active RAS is present at low level in BRAFV600E mutants, and mutant BRAF can continue signaling independent of RAS.
4. Addition of Vemurafenib inhibits BRAFV600E, lowers ERK level and relieves feedback on RTK and SOS.
5. RAS is activated, which induces dimerization of BRAF.
6. Vemurafenib binds one protomer in the BRAF dimer effectively, but binds the second protomer in the dimer with significantly lower affinity.
7. As a result, BRAF dimers can activate ERK via the promoter that is not bound to Vemurafenib, resulting in the partial restoration of ERK phosphorylation (imperfect adaptation)

::

   # The RAF-RAF complex binds Vemurafenib (1st Vemurafenib binds)
   Rule('First_binding_Vemurafenib',
	RAF(vem=None) % RAF(vem=None) + Vem(raf=None) <>
	RAF(vem=1) % RAF(vem=None) % Vem(raf=1), kef, ker)

   # The RAF-RAF complex binds Vemurafenib (2nd Vemurafenib binding
   Rule('Second_binding_vemurafenib',
	RAF(vem=None) % RAF(vem=ANY) + Vem(raf=None) <>
	RAF(vem=1) % RAF(vem=ANY) % Vem(raf=1), kff, kfr)

   # RAF that is not bound to RAF binds Vemurafenib
   Rule('Vemurafenib_binds_RAF_monomer',
	RAF(vem=None, d=None) + Vem(raf=None) <>
	RAF(vem=1, d=None) % Vem(raf=1), kef, ker)

   # Release RAS:GDP from RAF
   Rule('RAS_GDP_dissoc_RAF',
	RAS(gtp=None, raf=1) % RAF(ras=1) >>
	RAS(gtp=None, raf=None) + RAF(ras=None), koff)

   def observables():    
       # Observables
       # ----------
       Observable('RAF_WT_active',
		  RAF(d=ANY, vem=None)) 

       Observable('RAF_V600E_active',
		  RAF(vem=None))
       Observable('active_KRAS', RAS(gtp=ANY))
       Observable('active_SOS1', SOS1(state='up'))
       Observable('ERK_P', ERK(state='p'))    
       Observable('MEK_P', MEK(state='p'))


     

References
----------

.. [PMID11237210] Avruch J, Khokhlatchev A, Kyriakis JM, Luo Z, Tzivion G, Vavvas D, Zhang XF.  **Ras activation of the Raf kinase: tyrosine kinase recruitment of the MAP kinase cascade.** Recent Prog Horm Res. 2001;56:127-55. :pmid:`11237210`. :download:`PDF </pdf/11237210.pdf>`

.. [PMID21862573] Hibino K, Shibata T, Yanagida T, Sako Y. **Activation kinetics of RAF protein in the ternary complex of RAF, RAS-GTP, and kinase on the plasma membrane of living cells: single-molecule imaging analysis.** J Biol Chem. 2011 Oct 21;286(42):36460-8. :doi:`10.1074/jbc.M111.262675.` :pmid:`21862573` :download:`PDF </pdf/21862573.pdf>`

.. [PMID11447113] Chong H, Lee J, Guan KL. **Positive and negative regulation of Raf kinase activity and function by phosphorylation.** EMBO J. 2001 Jul 16;20(14):3716-27. :pmid:`11447113` :download:`PDF </pdf/11447113.pdf>`

.. [PMID15664184] Dumaz N, Marais R. **Raf phosphorylation: one step forward and two steps back.** Mol Cell. 2005 Jan 21;17(2):164-6. :pmid:`15664184` :download:`PDF </pdf/15664184.pdf>`

.. [PMID15664191] Dougherty MK1, Müller J, Ritt DA, Zhou M, Zhou XZ, Copeland TD, Conrads TP, Veenstra TD, Lu KP, Morrison DK. **Regulation of Raf-1 by direct feedback phosphorylation.** Mol Cell. 2005 Jan 21;17(2):215-24. :pmid:`15664191` :download:`PDF </pdf/15664191.pdf>`

.. [Lavoie] Lavoie H, Therrien M. **Regulation of RAF protein kinases in ERK signalling.** :doi:`10.1038/nrm3979` :download:`PDF </pdf/lavoie.pdf>`

.. [PMID2634358] Yao Z, Torres NM, Tao A, Gao Y, Luo L, Li Q, Stanchina E, Abdel-Wahab O, Solit DB, Poulikakos PI, Rosen N. **BRAF mutants evade ERK-dependent feedback by different mechanisms that determine their sensitivity to pharmacological inhibition.** Cancer Cell. 2015 Sept 1 15;28:270-83. :pmid:`26343582`

.. [PMID2420239] Lito P, Rosen N, Solit DB. ** Tumor adaptation and resistance to RAF inhibitors.** Nature MEdicine. 2013 Nov; 19(11):1401-9. :pmid:`24202393`

.. [PMID23153539] Lito P, Pratilas CA, Joseph EW, Tadi M, Halilovic E, Zubrowski M, Huan A, Wong WL, Callahan MK, Merghoun T, Wolchok JD, Stanchina E, Chandrarlapaty S, Paulikakos PI, Fagin JA, Rosen N, **Relief of profound feedback inhibition of mitogenic signaling by RAF inhibitors attenuates their activity in BRAFV600E melanomas.** Cancer Cell, 2012 Nov 12;22:668-82. pmid:`23153539`

.. [PMID21107323] Nazarian R, Shi H, Wanf Q, Kon X, Koya RC, Lee H, Chen Z, Lee M-K, Attar N, Sazegar H, chodon T, Nelson SF, McArthur G, Sosman JA, Ribas A, Lo RS. ** Melanomas acquire resistance to B-RAF(V600E) inhibition by RTK or N-RAS upregulation.** Nature. 2010 DEc 16;468: 973-7. :pmid:`2110732`
