   set(LATEX_OUTPUT_PATH ..)
   set(LATEX_DEFAULT_BUILD "pdf")
   ADD_LATEX_DOCUMENT( UppASDmanual.tex
   INPUTS structure.tex
   #BIBFILES mybib.bib
   IMAGE_DIRS Pictures)
