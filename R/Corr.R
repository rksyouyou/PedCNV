##' TODO
##'
##' TODO
##' @title TODO
##' @param ped TODO
##' @param type TODO
##' @return TODO
##' @author Meiling
##' @export
Corr <- function (ped, type=c("empirical","theoretical","ibs","hybrid")){
    td = tempdir()
    tf = tempfile("CORR", tmpdir = td)
    wisardPath = system.file("bin", "wisard", package="PedCNV")
    type = match.arg(type)
    switch(type,
           theoretical = system(paste(wisardPath, "--in", ped, "--pddt --cormat --out", tf)),
           empirical = system(paste(wisardPath, "--in", ped, "--cormat --out", tf)),
           ibs = system(paste(wisardPath, "--in", ped, "--ibs --cormat --out", tf)),
           hybrid=system(paste(wisardPath, "--in", ped, "--hybrid --cormat --out", tf))
           )
    dat = read.table(paste0(tf, ".pdt"))
    return(dat)
}
