Model reduction
===============

SOS can bind Grb2 only when it is unphosphorylated
::
   
   r = model.rules['GRB2_EGFR_SOS1_bind']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   SOS1(grb2=None, ras=None, state='up')
   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   SOS1(grb2=1, ras=None, state='up')

   r = model.rules['GRB2_EGFR_SOS1_dissociate']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[1] = SOS1(grb2=1, ras=None, state='up')
   r.product_pattern.complex_patterns[1].monomer_patterns[0] = SOS1(grb2=None, ras=None, state='up') 
   
SOS1 binds RAS when it is not bound to its substrate RAF or bound to GTP
::

   r = model.rules['SOS1_GRB2_RAS_bind']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   RAS(sos1=None, gtp=None, raf=None)
   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   RAS(sos1=1, gtp=None, raf=None)

   r.reactant_pattern.complex_patterns[0].monomer_patterns[0] =\
   SOS1(grb2=ANY, ras=None, state='up')
   r.product_pattern.complex_patterns[0].monomer_patterns[0] =\
   SOS1(grb2=ANY, ras=1, state='up')

   
SOS1 and RAS dissociate only when SOS1 is bound to Grb2
::
   r = model.rules['SOS1_GRB2_RAS_dissociate']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[0] =\
   SOS1(ras=1, grb2=ANY, state='up')
   r.product_pattern.complex_patterns[0].monomer_patterns[0] =\
   SOS1(ras=None, grb2=ANY, state='up')

   r.reactant_pattern.complex_patterns[0].monomer_patterns[1] =\
   RAS(sos1=1, gtp=None, raf=None)
   r.product_pattern.complex_patterns[1].monomer_patterns[0] =\
   RAS(sos1=None, gtp=None, raf=None)

Dissociation of RAS:GTP and SOS1
::

   r = model.rules
   r.add(Rule('RASGTP_SOS1_dissociate', RAS(gtp=ANY, sos1=1, raf=None) % SOS1(ras=1, grb2=ANY) >>\
   RAS(gtp=ANY, sos1=None, raf=None) + SOS1(ras=None, grb2=ANY), Parameter('ke4', 50)))

EGFR dimers dissociate only when bound to EGF
::

   r = model.rules['EGFR_EGF_EGFR_EGF_dissociate']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[0] =\
   EGFR(egfr=1, egf=ANY)
   r.reactant_pattern.complex_patterns[0].monomer_patterns[1] =\
   EGFR(egfr=1, egf=ANY)

   r.product_pattern.complex_patterns[0].monomer_patterns[0] =\
   EGFR(egfr=None, egf=ANY)
   r.product_pattern.complex_patterns[1].monomer_patterns[0] =\
   EGFR(egfr=None, egf=ANY)

EGFR binds Grb2 only when when EGFR is in the active dimer form.
::
   r = model.rules['GRB2_EGFR_bind']
   r.reactant_pattern.complex_patterns[1].monomer_patterns[0] =\
   EGFR(Y='p', grb2=None, egfr=ANY)

   r.product_pattern.complex_patterns[0].monomer_patterns[1] =\
   EGFR(Y='p', grb2=1, egfr=ANY)

Grb2 and EGFR dissociate only wwhen EGFR is in the active dimer form.
::
   r = model.rules['GRB2_EGFR_dissociate']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[1] =\
   EGFR(grb2=1, Y='p', egfr=ANY)
   r.product_pattern.complex_patterns[1].monomer_patterns[0] =\
   EGFR(grb2=None, Y='p', egfr=ANY)

Grb2 and SOS1 dissociate when Grb2 is bound to bound to EGFR
::
   r = model.rules['GRB2_EGFR_SOS1_dissociate']
   r.reactant_pattern.complex_patterns[0].monomer_patterns[0] =\
   GRB2(sos1=1, egfr=ANY)
   r.product_pattern.complex_patterns[0].monomer_patterns[0] =\
   GRB2(sos1=None, egfr=ANY)
   
