# Added these in for building of the model file
CODEDIR = src
DOCDIR = doc
DOCMODEL = $(DOCDIR)/doc_model
OUTPUTDIR = output
MODEL = $(OUTPUTDIR)/ras_model
PYSB_MODEL = $(MODEL).py
BNG_MODEL = $(MODEL).bngl
KAPPA_MODEL = $(MODEL).ka
SBML_MODEL = $(MODEL).sbml
RXN_NET = $(MODEL)_rxns

all: model bngl kappa sbml rxn_net contact_map influence_map simulation doc

doc: model
	cd doc; make html

clean:
	cd $(DOCDIR); make clean
	cd $(OUTPUTDIR); rm -f *.pdf *.png *.dot *.py *.ka *.bngl

model: $(PYSB_MODEL)
bngl: $(BNG_MODEL)
kappa: $(KAPPA_MODEL)
sbml: $(SBML_MODEL)

$(PYSB_MODEL): $(CODEDIR)/extract_model.py \
       $(DOCMODEL)/ras/ras_gtpase.rst \
       $(DOCMODEL)/ras/ras_gefs.rst \
       $(DOCMODEL)/ras/ras_gaps.rst \
       $(DOCMODEL)/raf/raf_activation_by_ras.rst \
       $(DOCMODEL)/initial_conditions.rst \
       $(DOCMODEL)/observables.rst
	mkdir -p $(OUTPUTDIR)
	python $(CODEDIR)/extract_model.py $(DOCMODEL)/ras/ras_gtpase.rst > $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCMODEL)/ras/ras_gefs.rst >> $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCMODEL)/ras/ras_gaps.rst >> $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCMODEL)/raf/raf_activation_by_ras.rst >> $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCMODEL)/initial_conditions.rst >> $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCMODEL)/observables.rst >> $(PYSB_MODEL)

%.bngl: %.py
	python -m pysb.export $< bngl > $(BNG_MODEL)

%.ka: %.py
	python -m pysb.export $< kappa > $(KAPPA_MODEL)

%.sbml: %.py
	python -m pysb.export $< sbml > $(SBML_MODEL)

rxn_net: $(PYSB_MODEL)
	python -m pysb.tools.render_reactions $(PYSB_MODEL) > $(RXN_NET).dot
	dot $(RXN_NET).dot -T pdf -o $(RXN_NET).pdf
	dot $(RXN_NET).dot -T png -o $(RXN_NET).png

contact_map: $(PYSB_MODEL) $(CODEDIR)/contact_map.py
	python $(CODEDIR)/contact_map.py $(PYSB_MODEL) $(MODEL)
	dot $(MODEL)_cm.dot -T pdf -o $(MODEL)_cm.pdf
	dot $(MODEL)_cm.dot -T png -o $(MODEL)_cm.png

influence_map: $(PYSB_MODEL) $(CODEDIR)/influence_map.py
	python $(CODEDIR)/influence_map.py $(PYSB_MODEL) $(MODEL)_im
	dot $(MODEL)_im_fixed.dot -T pdf -o $(MODEL)_im.pdf
	dot $(MODEL)_im_fixed.dot -T png -o $(MODEL)_im.png

simulation: $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/results/simulation.rst > $(OUTPUTDIR)/simulation.py
	python $(OUTPUTDIR)/simulation.py

