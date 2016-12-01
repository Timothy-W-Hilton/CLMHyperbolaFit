library(DEoptim)
library(chron)
library(Hmisc)  ## for monthDays
library(plyr)
library(lattice)
library(RColorBrewer)
library(latticeExtra)
library(ncdf4)

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

## ================================================================================
##                                data pre-processing
## ================================================================================

parse_monthly_data <- function(fpath='./monthly_vals.txt')
{
    ## Purpose: parse monthly_vals.txt to a data frame
    ##
    ## also calculates timestamps
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Timothy W. Hilton, Date: 23 Nov 2016, 09:14
    df <- read.csv(fpath)
    df[['X']] <- NULL   ## this column contains row number; delete it
    datestrs <- strsplit(as.character(df[['date']]), ' ')
    datestrs <- unlist(strsplit(as.character(df[['date']]), ' '))
    timestrs <- datestrs[seq(2, length(datestrs), 2)]
    datestrs <- datestrs[seq(1, length(datestrs), 2)]
    df[['date']] <- chron(datestrs, timestrs, format=c('y-m-d', 'h:m:s'))
    df[['ndays']] <- monthDays(df[['date']])
    return(df)
}

get_annual_pcp_npp <- function(df) {

    df[['var']] <- revalue(df[['var']], c("FPSN"="NPP"))
    df[['loc']] <- revalue(df[['loc']],
                           c("Sierra Foothill Research Extension Center"=
                                 "Sierra Foothill",
                             "Loma Ridge Global Change Experiment"=
                                 "Loma Ridge"))
    df <- df[df[['var']] %in% c('NPP', 'RAIN'), ]
    df[['year']] <- years(df[['date']])
    S_PER_DAY <- 60 * 60 * 24  ## seconds per day
    df[['monthsum']] <- df[['value']] * df[['ndays']] * S_PER_DAY
    npp_idx <- df[['var']] == 'NPP'
    MOL_PER_UMOL <- 1e-6
    C_G_PER_MOL <- 12
    df[npp_idx, 'monthsum'] <- df[npp_idx, 'monthsum'] * MOL_PER_UMOL * C_G_PER_MOL
    annsum <- aggregate(x=df[['monthsum']], by=df[, c('case', 'loc', 'var', 'year')], FUN=sum)
    annsum <- reshape(annsum, timevar = "var",
                      idvar = c("case", "loc", "year"),
                      direction = "wide")
    annsum <- rename(annsum, c('x.NPP'='NPP', 'x.RAIN'='RAIN'))
    return(annsum)
}

## ================================================================================
##                                fit curves to site data
## ================================================================================

fit_curves_to_site <- function(ad) {
    site_name <- as.character(unique(ad[['loc']]))
    if (length(site_name) != 1) {
        stop('annual data (ad) must contain one and only one site')
    }
    ctl <- DEoptim.control(itermax=1e3, trace=500, strategy=1)
    lower <- c(-100, -100, 0, 0, 1e-10)
    upper <- c(100, 100, 1000, 3000, 30)
    hyfit <- fit_WB_hyperbola(pcp=ad[['RAIN']], npp=ad[['NPP']], lower=lower, upper=upper, ctl=ctl)
    linfit <- fit_line(npp=ad[['NPP']], pcp=ad[['RAIN']])
    fit_comparison <- compare_lm_hyperbola(linfit, hyfit)
    ## fit_comparison[['aic_lin']] <- as.numeric(fit_comparison[['aic_lin']])
    ## fit_comparison[['aic_hy']] <- as.numeric(fit_comparison[['aic_hy']])
    ## return(list(hyfit=hyfit, linfit=linfit, comparison=fit_comparison))
    hypars <- hyfit[['optim']][['bestmem']]
    result <- data.frame(loc=site_name,
                         theta_1=hypars[[1]],
                         theta_2=hypars[[2]],
                         x_0=hypars[[3]], beta_0=hypars[[4]],
                         delta=hypars[[5]],
                         m = coef(linfit)[[1]], b=coef(linfit)[[2]],
                         AIC.hy=fit_comparison[['aic_hy']],
                         AIC.lin=fit_comparison[['aic_lin']])
    class(result) <- append(class(result), 'hyperbolafit')
    return(result)
}

## ================================================================================
##                                California site fits
## ================================================================================

predict_hy <- function(df) {
    df[['best_hy']] <- WB_hyperbola(df[['RAIN']],
                                    df[['theta_1']],
                                    df[['theta_2']],
                                    df[['x_0']],
                                    df[['beta_0']],
                                    df[['delta']])
    df[['best_lin']] <- df[['RAIN']] * df[['m']] + df[['b']]
    return(df)
}

fit_california_sites <- function() {
    md <- parse_monthly_data()  ## monthly data
    ad <- get_annual_pcp_npp(md)  ## annual data
    fits <- ddply(ad, "loc", fit_curves_to_site)
    ad_fits <- merge(fits, ad)
    ad_mod <- ddply(ad_fits, "loc", predict_hy)  ## annual data, modeled
    ad_mod <- ddply(ad_mod, "loc", function(x) x[order(x[['RAIN']]), ])
    ad_mod <- within(ad_mod, best_hy[AIC.hy > AIC.lin] <- NA)

    pal <- brewer.pal(n=3, name='Dark2')
    plt <- xyplot(NPP~RAIN|loc, groups=case, data=ad_mod,
                  xlab=expression(Rain~(mm~yr^{-1})),
                  ylab=expression(NPP~(g~C~m^{-2}~yr^{-1})),
                  col=pal[1:2], pch=c(24, 25),
                  key=list(text=list(levels(ad[['case']])), space='top',
                           points=list(pch=c(24, 25)), col=pal[1:2],
                           ## lines=list(col=pal[1:2]),
                           columns=nlevels(ad[['case']])))
    plt <- plt + xyplot(best_hy~RAIN|loc,
                        groups=case,
                        data=ad_mod,
                        type=c('s', 'l'),
                        col=pal[1:2],
                        panel=function(x, y, ...) {
                            panel.lines(x, y, col.line='black')})
    return(plt)
}

## ================================================================================
##                                read CLM netCDF
## ================================================================================

fit_global <- function() {
    nc <- nc_open('NPP_RAIN_annual.nc')
    nppctl <- ncvar_get(nc, 'NPPctl')
    nppide <- ncvar_get(nc, 'NPPide')
    rainctl <- ncvar_get(nc, 'RAINctl')
    rainide <- ncvar_get(nc, 'RAINide')
    nc_close(nc)

    fits <- vector("list", prod(dim(nppctl)[1:2]))
    n <- 1
    n_land_point <- 0
    cat('begin ', date(), '\n')
    for (i in seq(1, dim(nppctl)[1])) {
        for (j in seq(1, dim(nppctl)[2])) {
            this_ad <- data.frame(loc=paste(i, j, sep='_'),
                                  NPP=c(nppctl[i, j,], nppide[i, j,]),
                                  RAIN=c(rainctl[i, j,], rainide[i, j,]))
            land_point <- (!(all(is.na(this_ad[['NPP']]))) &
                           any(this_ad[['NPP']] > 0.0))
            if (land_point) {
                fits[[n]] <- fit_curves_to_site(this_ad)
                n_land_point <- n_land_point + 1
            }
            ## if (n_land_point > 5) {
            ##     break
            ## }
            n <- n + 1
            if ((n %% 100) == 0) {
                cat('n: ', n, '\n')
            }
        }
    }
    fitsdf <- do.call('rbind', fits)
    cat('done ', date(), '\n')
    return(fitsdf)
}

fillvals <- function(df, valscol, nrows, ncols) {
    result <- matrix(data=NA, nrow=nrowsCLM, ncol=ncolsCLM)
    result[ cbind(df$i, df$j) ] <- df[[valscol]]
    return(result)
}

fitsdf_ncdf <- function(fitsdf, fname_nc) {
    fitsdf[['hy_best']] <- fitsdf[['AIC.hy']] < fitsdf[['AIC.lin']]
    fitsdf[['AICrat']] <- fitsdf[['AIC.hy']] / fitsdf[['AIC.lin']]
    coords <- strsplit(as.character(fitsdf[['loc']]), '_')
    fitsdf[['i']] <- sapply(coords, function(x) as.numeric(x[[1]]))
    fitsdf[['j']] <- sapply(coords, function(x) as.numeric(x[[2]]))
    nrowsCLM <- max(fitsdf[['i']])
    ncolsCLM <- max(fitsdf[['j']])
    fieldnames <- c("theta_1", "theta_2", "x_0", "beta_0", "delta", "m",
                    "b", "AIC.hy", "AIC.lin", 'hy_best', 'AICrat')
    fields <- vector(mode='list', length=length(fieldnames))
    names(fields) <- fieldnames
    ncvars <- vector(mode='list', length=length(fieldnames))
    names(ncvars) <- fieldnames
    latdim <- ncdim_def(name='lat', units='index', vals=seq(1, nrowsCLM))
    londim <- ncdim_def(name='lon', units='index', vals=seq(1, ncolsCLM))
    for (this_name in names(fields)) {
        ncvars[[this_name]] <- ncvar_def(name=this_name,
                                         units="parameter values",
                                         dim=list(latdim, londim),
                                         missval=NA,
                                         longname="",
                                         prec="double")
    }
    ncnew <- nc_create(filename=fname_nc, vars=ncvars)
    for (this_name in names(fields)) {
        fields[[this_name]] <- fillvals(fitsdf, this_name, nrowsCLM, ncolsCLM)
        ncvar_put(nc=ncnew,
                  varid=ncvars[[this_name]],
                  vals=fields[[this_name]])
    }
    nc_close(ncnew)
    return(fields)
}

## save(fitsdf, file='fitsdf.RData')
load('fitsdf.RData')
