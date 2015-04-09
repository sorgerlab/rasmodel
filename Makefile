# Added these in for building of the model file
CODEDIR = src
DOCDIR = doc
OUTPUTDIR = output
MODEL = $(OUTPUTDIR)/ras_model
PYSB_MODEL = $(MODEL).py
BNG_MODEL = $(MODEL).bngl

doc: model
	cd doc; make html

clean:
	cd $(DOCDIR); make clean
	rm -f $(PYSB_MODEL)
	rm -f $(BNG_MODEL)

model: $(PYSB_MODEL)
bngl: $(BNG_MODEL)

$(PYSB_MODEL): $(CODEDIR)/extract_model.py \
       $(DOCDIR)/ras/ras_gtpase.rst \
       $(DOCDIR)/ras/ras_gefs.rst \
       $(DOCDIR)/ras/ras_gaps.rst \
       $(DOCDIR)/raf/raf_activation_by_ras.rst
	mkdir -p $(OUTPUTDIR)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/ras/ras_gtpase.rst > $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/ras/ras_gefs.rst >> $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/ras/ras_gaps.rst >> $(PYSB_MODEL)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/raf/raf_activation_by_ras.rst >> $(PYSB_MODEL)

%.bngl: %.py
	python -m pysb.export $< bngl > $(BNG_MODEL)
