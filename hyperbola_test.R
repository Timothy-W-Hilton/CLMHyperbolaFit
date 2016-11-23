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
    ## doesn't work for n=1 for some reason, need to figure out why
    ## generate random hyperbola parameters within specified bounds
    parvals <- mapply(function(LB, UB, n){runif(n, LB, UB)},
                   lower, upper, MoreArgs = list(n=n))
    parvals <- as.data.frame(t(parvals))
    ## generate hyperbola from the random parameters
    data <- lapply(parvals, function(p, pcp){
        WB_hyperbola(pcp, p[[1]], p[[2]], p[[3]], p[[4]], p[[5]])},
                   pcp=pcp)
    ## add some noise to the psuedodata
    data <- lapply(data, function(x) return(x + rnorm(length(x), sd = max(abs(x)) / 20.0)))
    return(list(pars=parvals, data=data))
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

fit.AIC <- function(nobs, npars, sse) {
  aic <- nobs * log(sse) + 2*npars
  return(aic)
}

fit_line <- function(npp, pcp) {
    linfit <- lm(npp~pcp, data=data.frame(pcp=pcp, npp=npp))
    sse <- sum(linfit[['residuals']]^2)
    return(linfit)
}

compare_lm_hyperbola<- function(linfit, hyfit)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Timothy W. Hilton, Date: 21 Nov 2016, 15:11
    aic_lin <- fit.AIC(nobs=length(linfit[['residuals']]),
                       npars=length(coef(linfit)),
                       sse=sum(linfit[['residuals']]^2))
    aic_hy <- fit.AIC(nobs=length(linfit[['residuals']]),
                      npars=length(hyfit[['optim']][['bestmem']]),
                      sse=hyfit[['optim']][['bestval']])
    return(list(aic_lin=aic_lin, aic_hy=aic_hy))
}


plot_pd_fit <- function(pcp, npp_pd, pd_pars, fit, linfit, aic)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Timothy W. Hilton, Date: 21 Nov 2016, 12:38
    my_lty <- c(NA, 1, 2, 1, 2)
    if (aic[['aic_lin']] < aic[['aic_hy']]) {
        my_lty <- c(NA, 2, 1, 1, 2)
    }

    pars <- fit[['optim']][['bestmem']]
    npp_mod <- WB_hyperbola(pcp, pars[[1]], pars[[2]], pars[[3]], pars[[4]], pars[[5]])
    par(mar=c(12.1, 4.1, 4.1, 2.1), xpd=TRUE)
    plot(pcp, npp_pd, ylab='annual NPP', xlab='annual pcp (mm)')
    points(pcp, npp_mod, type='l', lwd=5, col='#1b9e77', lty=my_lty[2])
    points(pars[[3]], pars[[4]], pch=13, col='#1b9e77', cex=3.0, lwd=5.0)
    lines(pcp, predict(linfit), col='#d95f02', lwd=5, lty=my_lty[3])
    ##

    legend(x='bottomleft',
           legend=c('pseudodata', 'best hyperbola', 'best line', 'AIC chose', 'AIC rejected'),
           col=c('black', '#1b9e77', '#d95f02', 'black', 'black'),
           lty=my_lty, pch=c(1, 13, NA, NA, NA), lwd=c(NA, 2, 2, 2, 2),
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

pseudodata_main <- function()
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Timothy W. Hilton, Date: 23 Nov 2016, 09:13

    ndata <- 400
    pcp <- seq(0, 1000, length.out=ndata)

    ctl <- DEoptim.control(itermax=1e3, trace=500, strategy=1)
    lower <- c(-100, -100, 0, 0, 1e-10)
    upper <- c(100, 100, 1000, 3000, 30)

    n_pseudo <- 100
    pd <- generate_pseudodata(pcp, lower, upper, n_pseudo)
    hyfits <- lapply(pd[['data']], fit_WB_hyperbola, pcp=pcp, lower=lower, upper=upper, ctl=ctl)
    linfits <- lapply(pd[['data']], fit_line, pcp=pcp)
    fit_comparison <- as.data.frame(t(mapply(compare_lm_hyperbola, linfits, hyfits)))
    fit_comparison[['aic_lin']] <- as.numeric(fit_comparison[['aic_lin']])
    fit_comparison[['aic_hy']] <- as.numeric(fit_comparison[['aic_hy']])

    pdf(file='pseudodata.pdf')
    for (i in seq(1, n_pseudo)){
        plot_pd_fit(pcp=pcp, npp_pd=pd[['data']][[i]],
                    pd_pars=pd[['pars']][[i]],
                    fit=hyfits[[i]],
                    linfit=linfits[[i]],
                    aic=fit_comparison[i, ])
    }
    dev.off()
}
