# Rscript for aseInferenceSingle
# Adopted from QuASAR

aseInferenceSingle <- function (gts, eps.vect, priors, ref.mat, alt.mat, min.cov, sample.name, annos) {
    ##################################################################
    ## inference
    ##################################################################
    coverage <- (ref.mat + alt.mat)
    coverage.floor <- min.cov
    coverage.ind <- (coverage > coverage.floor)
    ref <- ref.mat[coverage.ind]
    alt <- alt.mat[coverage.ind]
    phi <- priors[coverage.ind]
    eps <- eps.vect
    het <- gts[coverage.ind, 2]
    het.ind <- (het > 0.99)
    numb.hets <- sum(het.ind)
    annotations <- annos[coverage.ind, ][het.ind, ]

    cat("============================================================\n", sep = "")
    cat("========== Processing Sample: ", sample.name, " ==========\n", sep = "")
    cat("========== ", numb.hets, " heterozygotes with [P(het)>.99] ==========", "\n", sep = "")
    
    ##################################################################
    ## rho ~ calculate rho using eps
    ## rho0 ~ simple rho estimate
    ## lrt ~ likelihood ratio test using the simple rho estimate
    ## pval ~ pval of the llk ratio test using the chi-square approximation
    rho <- comp.rho(ref, alt, eps)
    rho0 <- plogis(log(ref) - log(alt))
    lrt <- lrtEpsRhoBinom(ref, alt, eps, rho)
    pval <- (1 - pchisq(2 * lrt, df = 1))


    ##################################################################
    ## D ~ a grid of possible dispersion values for the Beta-binomial model
    ## aux ~ loglikelihood of the beta bionmial model across D values
    ## with null \rho value
    ## Dmax ~ disperison which maximizes the llk
    D <- exp((0:500)/50)
    aux <- sapply(D, function(D) {
        sum(logLikBetaBinomialRhoEps(0.5, eps, D, ref[het.ind], alt[het.ind]))
    })
    Dmax <- D[which.max(aux)]

    cat("========== Dispersion estimate: ", round(Dmax, 3), " ==========\n", sep = "")

    ##################################################################
    ## Find MLE for \rho using the the Beta-Binomial model
    ## Dmax2 ~ the dispersian parameter estimed from the llk in the
    ##   previous step
    ## aux2 ~ optimization of the Beta-biomial model in terms of
    ##            logit(\rho)
    ## rho3 ~ vector of rho estimates from expit(logit(\rho))
    ## lrt3 ~ Recaluclate Het LRT using the beta-bionomial
    ## pval3 ~ pval of thr llk ratio test using the chi-square approximation
    ## betas.beta.binom ~ logit(\rho) or beta value for heterozygotes
    ## betas.se ~ standard error of the beta value
    ## betas.z ~ Z scores of the heterozygote beta values
    ## betas.pval ~ pvalues for the above z-scores
    Dmax2 <- Dmax
    aux2 <- t(sapply(1:sum(het.ind), function(ii) {
                auxLogis <- optim(0, fn = logLikBetaBinomial2, 
                                    gr = gLogLikBetaBinomial,
                                    D = Dmax2,
                                    R = ref[het.ind][ii],
                                    A = alt[het.ind][ii],
                                    method = "L-BFGS-B", hessian = TRUE)
                c(auxLogis$par, 1/(auxLogis$hessian)^0.5)
                }))

    rho3 <- plogis(aux2[, 1])
    betas.beta.binom <- aux2[, 1]
    lrt3 <- logLikBetaBinomialRhoEps(rho3, eps, Dmax2, ref[het.ind], alt[het.ind]) -
        logLikBetaBinomialRhoEps(0.5, eps, Dmax, ref[het.ind], alt[het.ind])
    pval3 <- (1 - pchisq(2 * lrt3, df = 1))
    betas.se <- abs(betas.beta.binom/qnorm(pval3/2))
    betas.se[which(betas.se == "NaN")] <- aux2[, 2][which(betas.se == "NaN")]

    betas.z <- betas.beta.binom/betas.se
    betas.pval <- 2 * pnorm(-abs(betas.z))

    ##################################################################
    ## rho2 ~ reassign rho calculated with a simple estimate from epsilon
    ## rho2[het.ind] ~ heterozygotes are re-assigned rho values
    ##         calucalted from grid optimization above
    ## aux ~ choose the null model whith largest probability
    ## lrt2 ~ llkRT with Beta-bionmial and the possibility that
    ##        that the heterozygote is of a different genotype
    ## pval2 ~ pval of thr llk ratio test using the chi-square approximation
    ##         (includes uncertainty in genotyping)
    ## qv2 ~ qvalue object calucalted from the llkRT
    ## qv2.qvals ~ qvalues from the previous calculation
    rho2 <- rho
    rho2[het.ind] <- rho3
    aux <- pmax(logLikBetaBinomialRhoEps(0, eps, Dmax, ref, alt), 
                logLikBetaBinomialRhoEps(1, eps, Dmax, ref, alt), 
                logLikBetaBinomialRhoEps(0.5, eps, Dmax, ref, alt))
    lrt2 <- logLikBetaBinomialRhoEps(rho2,eps,Dmax2,ref,alt) - aux;
    pval2 <- (1 - pchisq(2*lrt2, df=1))

    rsID <- annotations
    betas <- betas.beta.binom
    inference.data <- data.frame(annotations$rsID, annotations$chr,
        annotations$pos0, ref[het.ind], alt[het.ind], rho2[het.ind], betas, betas.se, pval2[het.ind])
    #colnames(inference.data)<-c("rsID", "chr", "pos0", paste(sample.name, c("beta", "se", "pval"), sep="_"))
    colnames(inference.data)<-c("rsID", "chr", "pos0", "R", "A", "rho", "beta", "se", "pval")

    return(inference.data)
}
