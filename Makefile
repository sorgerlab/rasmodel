CODEDIR = src
DOCDIR = doc
OUTPUTDIR = output
OUTPUTFILE = $(OUTPUTDIR)/extracted.py

all: doc model

.PHONY: clean doc model

clean:
	cd $(CODEDIR); make clean
	rm -f $(OUTPUTFILE)

doc:
	cd $(DOCDIR); make html

model: $(DOCDIR)/ras_regulation.rst $(CODEDIR)/extract_model.py
	mkdir -p $(OUTPUTDIR)
	python $(CODEDIR)/extract_model.py $(DOCDIR)/ras_regulation.rst > $(OUTPUTFILE)

