check:
    R -e 'devtools::check()'
render:
    pkgload::load_all(".", quiet = TRUE); rmarkdown::render("README.Rmd", quiet = TRUE)
doc:
    R -e 'devtools::document()'
