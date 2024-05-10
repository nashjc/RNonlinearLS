#!/bin/bash
## make article
Rscript -e 'rmarkdown::render("2023-15-copyeditresponse.Rmd", "rjtools::rjournal_web_article", output_options=list(), clean=FALSE)'
echo "Files:"
ls
read tmp
pdflatex RJwrapper.tex
pdflatex RJwrapper.tex
