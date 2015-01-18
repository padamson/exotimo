# configure the next two lines
SPEC=ExotiMO
SEAMLESSHOME=~/research/seamless
SPECHOME=~/research/exotimo/spec


DATE=$(shell date +'%H%M.%d.%m.%Y')

all:
	touch $(SPEC).ind
	pdflatex -file-line-error -halt-on-error $(SPEC)
	makeindex $(SPEC)
	bibtex $(SPEC).aux
	pdflatex -file-line-error -halt-on-error $(SPEC)
	pdflatex -file-line-error -halt-on-error $(SPEC)

usedvi:
	touch $(SPEC).ind
	latex $(SPEC)
	makeindex $(SPEC)
	bibtex $(SPEC).aux
	latex $(SPEC)
	latex $(SPEC)
	dvips -P cmz -t letter -o $(SPEC).ps $(SPEC).dvi
	ps2pdf $(SPEC).ps $(SPEC).pdf

publish:
	cp $(SPEC).pdf release/$(SPEC)_release_$(DATE).pdf

draft:
	cp $(SPEC).pdf release/$(SPEC)_draft_$(DATE).pdf

clean:
	rm -f *.aux *.blg *.bbl *.dvi *.log

clobber: clean
	rm -f $(SPEC).ps $(SPEC).pdf $(SPEC).ilg $(SPEC).out $(SPEC).toc $(SPEC).idx $(SPEC).ind

TESTOPTS=-no-chpl-home-warn -compopts '-M $(SPECHOME)/src -M $(SPECHOME)/helper'
test: FORCE
	cd test && start_test $(TESTOPTS)

SPECHELPER_DIR= ./helper
helpers: FORCE
	rm -rf $(SPECHELPER_DIR)
	$(SEAMLESSHOME)/util/extract_helpers -o $(SPECHELPER_DIR) *.tex

SPECSOURCE_DIR = ./src
sources: FORCE
	rm -rf $(SPECSOURCE_DIR)
	$(SEAMLESSHOME)/util/extract_sources -o $(SPECSOURCE_DIR) *.tex

SPECTEST_DIR = ./test
tests: FORCE
	rm -rf $(SPECTEST_DIR)
	$(SEAMLESSHOME)/util/extract_tests -o $(SPECTEST_DIR) *.tex

speclist: FORCE
	sed -n '/\\item\\spec/p' *.tex

reqlist: FORCE
	grep '\\req' $(SPECHOME)/Requirements.tex

FORCE:
