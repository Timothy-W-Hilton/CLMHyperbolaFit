library(DEoptim)

#==============================================================
##' calculates sum of squares of X.
##'
##' if there are no valid data points, return an arbitrary (and very
##' high) SSE. This is a little bit of a hack - it really ought to
##' return NA in that situation, but it can"t return NA, because the
##' DEoptim implementation of DE requires that the optimized function
##' return a non-NA scalar.
##' @title calculate sum of squares
##' @param x values for which to calculate sum of squares
##' @return sum of squares of x
##' @author Timothy W. Hilton
getSSE <- function(x) {


  if (length(which(!is.na(x))) == 0) return(1e20)
  else return(sum( x^2, na.rm=TRUE))
}

WB_hyperbola <- function(x, theta_1, theta_2, x_0, beta_0, delta)
{
    ## Purpose: implement the rotated hyperbola of Watts & Bacon (1974)
    ##
    ##
    ## from Watts and Bacon (1974):
    ## (i) the dependent variable y is a single valued function of the
    ##     independent variable x
    ## (ii) the left asymptote has slope theta_1
    ## (iii) the right asymptote has slope theta_2
    ## (iv) the asymptotes intersect at the point (x_0, beta_0),
    ## (v) the radius of curvature at x = x_0, is proproportional to a
    ##     quantity delta
    ##
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Timothy W. Hilton, Date: 18 Nov 2016, 13:39

    beta_1 <- (theta_1 + theta_2) / 2.0
    beta_2 <- (theta_2 - theta_1) / 2.0
    y <- beta_0 + beta_1*(x - x_0) + beta_2 * sqrt((x - x_0)^2 + (delta^2)/4)
    return(y)
}

WB_hyperbola_SSE <- function(pars, pcp, npp_obs) {
    ## Purpose: evaluate sum of squared errors (SSE) for WB_hyperbola()
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##    pars: parmeter vector.  in order, [x, theta_1, theta_2, x_0, beta_0]
    ## ----------------------------------------------------------------------
    ## Author: Timothy W. Hilton, Date: 18 Nov 2016, 13:58

    npp_est <- WB_hyperbola(pcp, pars[[1]], pars[[2]], pars[[3]], pars[[4]], pars[[5]])
    res <- npp_obs - npp_est
    sse <- getSSE(res)
    ## cat('SSE:', sse, '\n')
    return(sse)
}

generate_pseudodata <- function(pcp, lower, upper, n) {
    ## generate random hyperbola parameters within specified bounds
    pars <- mapply(function(LB, UB, n){runif(n, LB, UB)},
                   lower, upper, MoreArgs = list(n=n))
    pars <- as.data.frame(t(pars))
    ## generate hyperbola from the random parameters
    data <- lapply(pars, function(pars, pcp){
        WB_hyperbola(pcp,
                     pars[[1]], pars[[2]], pars[[3]],
                     pars[[4]], pars[[5]])
    },
                   pcp=pcp)
    ## add some noise to the psuedodata
    data <- lapply(data, function(x) return(x + rnorm(length(x), sd = max(abs(x)) / 20.0)))
    return(list(pars=pars, data=data))
}

fit_WB_hyperbola<- function(pcp, npp, upper, lower, ctl)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Timothy W. Hilton, Date: 21 Nov 2016, 12:30
    r <- DEoptim(fn=WB_hyperbola_SSE, lower=lower, upper=upper, control=ctl, pcp, npp)
    return(r)
}

plot_pd_fit <- function(pcp, npp_pd, pd_pars, fit)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Timothy W. Hilton, Date: 21 Nov 2016, 12:38
    pars <- fit[['optim']][['bestmem']]
    npp_mod <- WB_hyperbola(pcp, pars[[1]], pars[[2]], pars[[3]], pars[[4]], pars[[5]])
    par(mar=c(12.1, 4.1, 4.1, 2.1), xpd=TRUE)
    plot(pcp, npp_pd, ylab='annual NPP', xlab='annual pcp (mm)')
    points(pcp, npp_mod, type='l', lwd=3, col='#1b9e77')
    points(pars[[3]], pars[[4]], pch=13, col='#1b9e77', cex=3.0, lwd=5.0)
    ##
    legend(x='bottomleft',
           legend=c('pseudodata', 'best fit', 'inflection'),
           col=c('black', '#1b9e77', '#1b9e77'),
           lty=c(0, 1, 0), pch=c(1, NA, 13),
           inset=c(0.0, -0.6))
    ##
    df_pars <- data.frame(fit=fit[['optim']][['bestmem']], pseudo=pd_pars)
    df_pars[['par']] <- c('slope L', 'slope R', 'inflection (x)', 'inflection (y)', 'delta')
    df_pars <- df_pars[c('par', 'fit', 'pseudo')]
    leg_text <- apply(format(df_pars, digits=3, scientific=TRUE), 1,
                      function(x) paste(x, collapse=','))
    legend(x='bottomright',
           inset=c(0.0, -0.6),
           legend=leg_text)
}


ndata <- 400
pcp <- seq(0, 1000, length.out=ndata)
theta_1 <- 0.4
theta_2 <- 0.3
x_0 <- 500
beta_0 <- 450
delta <- 10
npp_obs <- WB_hyperbola(pcp, theta_1, theta_2, x_0, beta_0, delta) + rnorm(ndata, sd=5)

ctl <- DEoptim.control(itermax=1e3, trace=500, strategy=1)
lower <- c(-100, -100, 0, 0, 1e-10)
upper <- c(100, 100, 1000, 3000, 30)
## r <- DEoptim(fn=WB_hyperbola_SSE, lower=lower, upper=upper, control=ctl, pcp, npp_obs)
## pars <- r[['optim']][['bestmem']]
## npp_mod <- WB_hyperbola(pcp, pars[[1]], pars[[2]], pars[[3]], pars[[4]], pars[[5]])

n_pseudo <- 3
pd <- generate_pseudodata(pcp, lower, upper, n_pseudo)
fits <- lapply(pd[['data']], fit_WB_hyperbola, pcp=pcp, lower=lower, upper=upper, ctl=ctl)
pdf(file='pseudodata.pdf')
for (i in seq(1, n_pseudo)){
    plot_pd_fit(pcp=pcp, npp_pd=pd[['data']][[i]], pd_pars=pd[['pars']][[i]], fit=fits[[i]])
}
dev.off()

## pdf(file='pseudodata.pdf')
## for (i in seq(1, n_pseudo)){
##     plot(pcp, pd[['data']][[i]])
## }
## dev.off()



## SSE_correct <- WB_hyperbola_SSE(list(theta_1, theta_2, x_0, beta_0, delta),
##                                 pcp, npp_obs)
## plot(pcp, npp_obs, xlim=c(-100, 3000), ylim=c(-100, 3000))
## points(pcp, npp_mod, pch='*', col='red')
