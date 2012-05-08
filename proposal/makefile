# Only add file roots, not file extensions
LATEX_FILES = proposal

# Modify these commands to add command-line flags and options
LATEX = latex -interaction=nonstopmode
BIBTEX = bibtex
DVIPS = dvips -t letter
PS2PDF = ps2pdf

all: $(addsuffix .pdf, $(LATEX_FILES))
	-

%.pdf: %.tex
	$(LATEX) $*.tex
	-$(BIBTEX) $*
	$(LATEX) $*.tex
	$(LATEX) $*.tex
	$(DVIPS) $*.dvi
	$(PS2PDF) $*.ps

tidy:
	-rm $(addsuffix .aux, $(LATEX_FILES)) \
        $(addsuffix .bbl, $(LATEX_FILES)) \
        $(addsuffix .blg, $(LATEX_FILES)) \
        $(addsuffix .toc, $(LATEX_FILES)) \
        $(addsuffix .log, $(LATEX_FILES)) \
        $(addsuffix .ps, $(LATEX_FILES)) \
        $(addsuffix .dvi, $(LATEX_FILES)) \

clean: tidy
	-rm $(addsuffix .pdf, $(LATEX_FILES))