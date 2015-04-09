CODEDIR = src
DOCDIR = doc
OUTPUTDIR = output
OUTPUTFILE = $(OUTPUTDIR)/ras_model.py

all: doc model

.PHONY: clean doc model

clean:
	cd $(CODEDIR); make clean
	rm -f $(OUTPUTFILE)

doc: model
	cd $(DOCDIR); make html

model: $(CODEDIR)/extract_model.py \
       $(DOCDIR)/ras/ras_gtpase.rst \
       $(DOCDIR)/ras/ras_gefs.rst \
       $(DOCDIR)/ras/ras_gaps.rst \
       $(DOCDIR)/raf/raf_activation_by_ras.rst
	mkdir -p $(OUTPUTDIR)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/ras/ras_gtpase.rst > $(OUTPUTFILE)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/ras/ras_gefs.rst >> $(OUTPUTFILE)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/ras/ras_gaps.rst >> $(OUTPUTFILE)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/raf/raf_activation_by_ras.rst >> $(OUTPUTFILE)

