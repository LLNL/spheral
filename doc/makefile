TEXSOURCE = Building.tex

all:	$(TEXSOURCE:.tex=.dvi)
pdf:	$(TEXSOURCE:.tex=.pdf)
ps:	$(TEXSOURCE:.tex=.ps)
html:   $(TEXSOURCE:.tex=.html)

force_build:

%.dvi:	force_build
	latex $*.tex
	latex $*.tex

%.pdf:	$(TEXSOURCE)
	pdflatex $*.tex
	pdflatex $*.tex

%.html: $(TEXSOURCE)
	latex2html $*.tex
