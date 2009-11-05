`vcov.cusp` <-
function (object, ...) 
{
    so <- summary.cusp(object, corr = FALSE, logist = FALSE, 
        ...)
    so$dispersion * so$cov.unscaled
}

