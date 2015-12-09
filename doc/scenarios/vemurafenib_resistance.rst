.. component::
   :module: rasmodel.scenarios.vemurafenib_resistance
	    
Vemurafenib resistance
======================

::
   
    from pysb import *

    Model()

    import rasmodel.components.egfr as egfr
    egfr.egfr_minimal()

    import rasmodel.components.ras as ras
    ras.ras_minimal()

    import rasmodel.components.raf as raf
    raf.raf_monomers()
    raf.ras_activates_raf()
    raf.vemurafenib_monomers()
    raf.vemurafenib_binds_raf()
    raf.mek_phosphorylation()

    import rasmodel.components.mapk as mapk
    mapk.erk_dynamics()

  
