DOCDIR = doc

.PHONY: all tangle html clean

all: html tangle

tangle:
	cd $(DOCDIR); make tangle
	cd ..
	for f in rasmodel/{components,scenarios,experiments}/*.py; do \
	    [ `basename "$$f"` != __init__.py ] && rm -f "$$f"{,c} || true; \
	done
	cp -R $(DOCDIR)/_build/tangle/ .

html:
	cd $(DOCDIR); make html

clean:
	cd $(DOCDIR); make clean
