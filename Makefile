# Added these in for building of the model file
CODEDIR = src
DOCDIR = doc
OUTPUTDIR = output
MODEL = $(OUTPUTDIR)/ras_model
PYSB_MODEL = $(MODEL).py
BNG_MODEL = $(MODEL).bngl
KAPPA_MODEL = $(MODEL).ka
RXN_NET = $(MODEL)_rxns

all: model bngl kappa doc

doc: model
	cd doc; make html

clean:
	cd $(DOCDIR); make clean
	rm -f $(PYSB_MODEL)
	rm -f $(BNG_MODEL)

model: $(PYSB_MODEL)
bngl: $(BNG_MODEL)
kappa: $(KAPPA_MODEL)

$(PYSB_MODEL): $(CODEDIR)/extract_model.py \
       $(DOCDIR)/ras/ras_gtpase.rst \
       $(DOCDIR)/ras/ras_gefs.rst \
       $(DOCDIR)/ras/ras_gaps.rst \
       $(DOCDIR)/raf/raf_activation_by_ras.rst \
       $(DOCDIR)/initial_conditions.rst
	mkdir -p $(OUTPUTDIR)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/ras/ras_gtpase.rst > $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/ras/ras_gefs.rst >> $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/ras/ras_gaps.rst >> $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/raf/raf_activation_by_ras.rst >> $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/initial_conditions.rst >> $(PYSB_MODEL)

%.bngl: %.py
	python -m pysb.export $< bngl > $(BNG_MODEL)

%.ka: %.py
	python -m pysb.export $< kappa > $(KAPPA_MODEL)

rxn_net: $(PYSB_MODEL)
	python -m pysb.tools.render_reactions $(PYSB_MODEL) > $(RXN_NET).dot
	dot $(RXN_NET).dot -T pdf -o $(RXN_NET).pdf
	dot $(RXN_NET).dot -T png -o $(RXN_NET).png

contact_map: $(PYSB_MODEL) $(CODEDIR)/contact_map.py
	python $(CODEDIR)/contact_map.py $(PYSB_MODEL) $(MODEL)
	dot $(MODEL)_cm.dot -T pdf -o $(MODEL)_cm.pdf
	dot $(MODEL)_cm.dot -T png -o $(MODEL)_cm.png
