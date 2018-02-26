TEXSOURCE = Building.tex Manual.tex

all:	$(TEXSOURCE:.tex=.dvi)
pdf:	$(TEXSOURCE:.tex=.pdf)
ps:	$(TEXSOURCE:.tex=.ps)
html:   $(TEXSOURCE:.tex=.html)

force_build:

%.pdf:	$(TEXSOURCE)
	pdflatex $*.tex
	pdflatex $*.tex

%.dvi:	force_build
	latex $*.tex
	latex $*.tex

%.html: $(TEXSOURCE)
	latex2html $*.tex
