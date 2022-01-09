
my.knit = knitr::knit("document.Rnw")
## document.tex is the latex file that will be compiled by the two following command:

system(paste0("pdflatex ", "document.tex")) 
system(paste0("pdflatex ", "document.tex")) 