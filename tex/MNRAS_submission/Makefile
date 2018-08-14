
TEXFILE=Ross_etal_J1100_MNRAS
default: all

# just run it


bib: 
	bibtex ${TEXFILE}

all :
	pdflatex  ${TEXFILE}
	pdflatex  ${TEXFILE}
	bibtex ${TEXFILE}
	bibtex ${TEXFILE}
	bibtex ${TEXFILE}
	pdflatex  ${TEXFILE}
	pdflatex  ${TEXFILE}
	pdflatex  ${TEXFILE}
#	dvipdf ${TEXFILE}
	rm -f  ${TEXFILE}.aux
	rm -f  ${TEXFILE}.bbl
	rm -f  ${TEXFILE}.blg
	rm -f  ${TEXFILE}.dvi
	rm -f  ${TEXFILE}.log
	rm -f  ${TEXFILE}.out
	rm -f  ${TEXFILE}.toc
	rm -f  ${TEXFILE}Notes.bib
#	dvips ERQs_PaperI_rv2.dvi 
#	ps2pdf ERQs_PaperI_rv2.ps
#	rm -f ERQs_PaperI_rv2.ps


#	pdflatex  ${TEXFILE_TWO}
#	pdflatex  ${TEXFILE_TWO}
#	bibtex ${TEXFILE_TWO}
#	bibtex ${TEXFILE_TWO}
#	pdflatex  ${TEXFILE_TWO}
#	pdflatex  ${TEXFILE_TWO}
#	rm -f  ${TEXFILE_TWO}.aux
#	rm -f  ${TEXFILE_TWO}.bbl
#	rm -f  ${TEXFILE_TWO}.blg
#	rm -f  ${TEXFILE_TWO}.dvi
#	rm -f  ${TEXFILE_TWO}.log
#	rm -f  ${TEXFILE_TWO}.out
#	rm -f  ${TEXFILE_TWO}.toc
#	rm -f  ${TEXFILE_TWO}Notes.bib

	open ${TEXFILE}.pdf
#	open ${TEXFILE_TWO}.pdf


clean:
	rm -f \
	${TEXFILE}.dvi \
	${TEXFILE}.ps \
	${TEXFILE}.out \
	${TEXFILE}.aux \
	${TEXFILE}.bbl \
	${TEXFILE}.blg \
	${TEXFILE}.toc \
	${TEXFILE}.log 
#	${TEXFILE}.pdf 


sup:
	pdflatex  ${TEXFILE_TWO}
	pdflatex  ${TEXFILE_TWO}
	bibtex ${TEXFILE_TWO}
	bibtex ${TEXFILE_TWO}
	pdflatex  ${TEXFILE_TWO}
	pdflatex  ${TEXFILE_TWO}
	rm -f  ${TEXFILE_TWO}.aux
	rm -f  ${TEXFILE_TWO}.bbl
	rm -f  ${TEXFILE_TWO}.blg
	rm -f  ${TEXFILE_TWO}.dvi
	rm -f  ${TEXFILE_TWO}.log
	rm -f  ${TEXFILE_TWO}.out
	rm -f  ${TEXFILE_TWO}.toc
	rm -f  ${TEXFILE_TWO}Notes.bib
	open ${TEXFILE_TWO}.pdf

