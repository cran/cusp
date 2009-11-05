`pcusp` <-
function (y, alpha, beta) 
Vectorize(function(x) integrate(dcusp, -Inf, x, alpha = alpha, 
    beta = beta)$value)(y)

