script = 
`Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
if (!require("salso")) {
    install.packages("salso", lib = Sys.getenv("R_LIBS_USER"), repos = "http://cran.us.r-project.org")
}'`
run(script)