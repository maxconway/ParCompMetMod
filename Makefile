geobacter.pdf : ParCompMetMod.pdf figures/*

%.pdf : %.tex
	pdflatex $*.tex
	bibtex *.aux
	pdflatex $*.tex
	pdflatex $*.tex

clean : 
	rm -f ParCompMetMod.pdf ParCompMetMod.aux ParCompMetMod.bbl ParCompMetMod.blg ParCompMetMod.log ParCompMetMod.out ParCompMetMod.rtf *~

%.rtf : %.tex
	latex2rtf $*.tex

wordcount : ParCompMetMod.tex
	texcount ParCompMetMod.tex
