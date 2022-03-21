## Copyright (C) 2019 Robersy Sanchez <https://genomaths.com/>
##
## Author: Robersy Sanchez
#
## This file is part of the R package "usefr".
##
## 'usefr' is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

#' @rdname tableBoots
#'
#' @title tableBoots
#' @description Parametric Bootstrap of \eqn{n x m} Contingency independence
#' test. It is assumed that the provided contingency table summarizes the
#' results where two or more groups are compared and the outcome is a
#' categorical variable, i.e., disease vs. no disease, pass vs. fail, treated
#' patient vs. control group, etc.
#'
#' The goodness of fit statistic is the root-mean-square statistic (RMST)
#' or Hellinger divergence, as proposed by Perkins et al. [1, 2]. Hellinger
#' divergence (HD) is computed as proposed in [3].
#'
#' @param x A numerical matrix corresponding to cross tabulation \eqn{n x m}
#' table (contingency table).
#' @param stat Statistic to be used in the testing: 'rmst','hdiv', or 'all'.
#' @param num.permut Number of permutations.
#' @param out.stat logical(1). Whether to return the values of the statistics
#' used: the bootstrap mean and the original value estimated.
#'
#' @details For goodness-of-fit the following null hypothesis is tested
#'     \eqn{H_\theta: p = p(\theta)}
#'     To conduct a single simulation, we perform the following three-step
#'      procedure [1,2]:
#' \enumerate{
#'     \item To generate m i.i.d. draws according to the model distribution
#'           \eqn{p(\theta)}, where \eqn{\theta'} is the estimate calculated
#'           from the experimental data,
#'     \item To estimate the parameter \eqn{\theta} from the data generated in
#'           Step 1, obtaining a new estimate \eqn{\theta}est.
#'     \item To calculate the statistic under consideration (HD,
#'           RMST), using the data generated in Step 1 and taking the model
#'           distribution to be \eqn{\theta}est, where \eqn{\theta}est is the
#'           estimate calculated in Step 2 from the data generated in Step 1.
#' }
#'     After conducting many such simulations, the confidence level for
#'     rejecting the null hypothesis is the fraction of the statistics
#'     calculated in step 3 that are less than the statistic calculated from
#'     the empirical data. The significance level \eqn{\alpha} is the same as a
#'     confidence level of \eqn{1-\alpha}.
#'
#' @return A p-value probability
#' @references
#' \enumerate{
#'     \item Perkins W, Tygert M, Ward R. Chi^2 and Classical Exact Tests
#'           Often Wildly Misreport Significance; the Remedy Lies in Computers
#'           [Internet]. Uploaded to ArXiv. 2011. Report No.:
#'           arXiv:1108.4126v2.
#'     \item Perkins, W., Tygert, M. & Ward, R. Computing the confidence
#'           levels or a root-mean square test of goodness-of-fit. 217,
#'           9072-9084 (2011).
#'     \item Basu, A., Mandal, A. & Pardo, L. Hypothesis testing for two
#'           discrete populations based on the Hellinger distance. Stat.
#'           Probab. Lett. 80, 206-214 (2010).
#' }
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @examples
#' set.seed(123)
#' ## The numbers of methylated and unmethylted "read-counts" targeting DNA
#' ## cytosine sites are the typical raw data resulting from bisulfate
#' ## sequencing experiments (BiSeq). Methylation analysis requires for the
#' ## application of some statistical tests to identify differentially
#' ## methylated sites/ positions (DMSs or DMPs). BiSeq experiments on humans
#' ## are relatively expensive (still in 2020) and some researcher is willing
#' ## to accept a site coverage (total) of 4 read-counts as threshold for a
#' ## site to be included in the analysis. Let's suppose that at a given site
#' ## we have the following contingency table:
#'
#' site_res <- matrix(c(3, 1, 1, 3), nrow = 2,
#'                    dimnames = list(status = c("ctrl", "treat"),
#'                                    treat = c("meth", "unmeth")))
#' site_res
#'
#' ## The application of classical test variants yield the following results
#' fisher.test(site_res)$p.value
#' fisher.test(site_res, alternative = "greater")$p.value
#' fisher.test(site_res, simulate.p.value = TRUE)$p.value
#' chisq.test(site_res)$p.value
#' chisq.test(site_res, simulate.p.value = TRUE)$p.value
#'
#' ## That is, these test did not detect differences between the read-counts
#' ## for control and the treated patients. 'tableBoots' function suggests
#' ## that the site is borderline.
#'
#' tableBoots( site_res, stat = 'all', num.permut = 1e3 )
#'
#' ## A slight change in the read counts identify a meaningful difference
#' ## not detected by the classicla tests
#' site_res <- matrix(c(3, 0, 1, 4), nrow = 2,
#'                    dimnames = list(status = c("ctrl", "treat"),
#'                                    treat = c("meth", "unmeth")))
#' site_res
#'
#' set.seed(23)
#' tableBoots( site_res, stat = 'all', num.permut = 1e3 )
#' fisher.test(site_res)$p.value
#' fisher.test(site_res, alternative = "greater")$p.value
#' fisher.test(site_res, simulate.p.value = TRUE)$p.value
#' chisq.test(site_res)$p.value
#' chisq.test(site_res, simulate.p.value = TRUE)$p.value
#'
#' ## Now, let's suppose that we want to test whether methylation status from
#' ## two DNA sites the same differentially methylated region (DMR)
#' ## are independent, i.e. whether the treatment affected the methylation
#' ## status of the two sites. In this case, the counts grouped into four
#' ## categories.
#' set.seed(1)
#' site_res <- matrix(c(3, 1, 1, 3, 3, 0, 1, 4), nrow = 2, byrow = FALSE,
#'                    dimnames = list(status = c("ctrl", "treat"),
#'                                    treat = c("meth.1", "unmeth.1",
#'                                              "meth.2", "unmeth.2")))
#' site_res
#'
#' ## Chi-squared from the R package 'stats'
#' chisq.test(site_res)$p.bvalue
#' chisq.test(site_res, simulate.p.value = TRUE, B = 2e3)$p.value
#'
#' tableBoots( site_res, stat = 'all', num.permut = 999 )
#'
#' ## Results above are in border. If we include, third site,
#' ## then sinces would different.
#' site_res <- matrix(c(3, 1, 1, 3, 3, 0, 1, 4,  4, 0, 1, 4),
#'                    nrow = 2, byrow = FALSE,
#'                    dimnames = list(status = c("ctrl", "treat"),
#'                                    treat = c("meth.1", "unmeth.1",
#'                                              "meth.2", "unmeth.2",
#'                                              "meth.3", "unmeth.3")))
#' site_res
#'
#' ## That is, we have not reason to believe that the observed methylation
#' ## levels in the treatment are not independent
#' chisq.test(site_res)$p.value
#' chisq.test(site_res, simulate.p.value = TRUE, B = 2e3)$p.value
#' tableBoots( site_res, stat = 'all', num.permut = 999 )
#'
tableBoots <- function(
                        x,
                        stat = c("rmst", "hd", "chisq", "all"),
                        out.stat = FALSE,
                        num.permut = 100) {

    # obsf_ij is the observed cell count in the ith row and jth column of
    # the table expf_ij is the expected cell count in the ith row and jth
    # column of the table, computed as
    # expf_ij = row i total x col j / total-grand total

    m0 <- rowSums(x)
    n0 <- colSums(x)
    if ((N0 <- sum(x)) == 0)
        stop("\n*** at least one entry of 'x' must be positive")

    ## The expected number of counts
    expf <- as.vector(outer(m0, n0)) / N0 #
    prob <- expf / N0
    d <- dim(x)

    ## Function to randomly generate a nxm table based on the expected
    ## number of counts estimated from the observed nxm contingency
    ## table
    y <- rep(0, prod(d))  ## all initial counts equal to zero

    ## Randomly select a cell in the table with probability equal to
    ## the expected counts divided by the total number of counts (N) in
    ## the table & increment the value in this cell by one. Repeat this
    ## procedure N times

    tb <- replicate(num.permut, {
                        r <- table(sample(
                                        x = seq_len(prod(d)),
                                        size = N0,
                                        replace = TRUE,
                                        prob = prob))
                        ## updates initial counts
                        y[as.numeric(names(r))] <- r
                        y
                    },
                    simplify = FALSE
    )

    x <- as.vector(x)

    st <- lapply(tb, function(y) {
        ctb <- matrix(y, nrow = d[1])
        m <- rowSums(ctb)
        n <- colSums(ctb)
        N <- sum(ctb)
        freq <- as.vector(outer(m, n)) / N # Expected frequencies
        ## Compute the specified statistic for the randomly generated table
        st <- switch(stat,
                     rmst = sum((y - freq)^2, na.rm = TRUE)/4,
                     hd = HD(y, freq),
                     chisq = sum((y - freq)^2/freq, na.rm = TRUE),
                     all = c(rmst = sum((y - freq)^2, na.rm = TRUE)/4,
                             hdiv = HD(y, freq),
                             chisq = sum((y - freq)^2/freq, na.rm = TRUE))
        )

        return(st)
    })

    if (stat == "all") {
        st <- do.call(rbind, st)
        st0 <- c(
                    rmst = sum((x - expf)^2, na.rm = TRUE)/4,
                    hd = HD(x, expf),
                    chisq = sum((x - expf)^2/expf, na.rm = TRUE)
                )
        res <- c(
                rmst.p.value = (sum(st[, 1] > st0[1], na.rm = TRUE) + 1)/
                    (num.permut + 1),
                hdiv.p.value = (sum(st[, 2] > st0[2], na.rm = TRUE) + 1)/
                    (num.permut + 1),
                chisq.p.value = (sum(st[, 3] > st0[3], na.rm = TRUE) + 1)/
                    (num.permut + 1)
                )
    } else {
        st <- unlist(st)
        st0 <- switch(stat,
                      rmst = sum((x - expf)^2, na.rm = TRUE)/4,
                      hd = HD(x, expf),
                      chisq = sum((x - expf)^2/expf, na.rm = TRUE))
        if (out.stat) {
            p.value <- (sum(st > st0, na.rm = TRUE) + 1)/(num.permut + 1)
            boot.stat <- mean(c(st,st0), na.rm = TRUE)
            res <- data.frame(stat = st0, boot.stat = boot.stat,
                            p.value = p.value)
        }
        else
            res <- (sum(st > st0, na.rm = TRUE) + 1)/(num.permut + 1)
    }
    return(res)
}


# ========================= Auxiliary function ============================ #

# ------------------ Hellinger divergence  statistic------------------------ #

HD <- function(x, y) {

    ## Function to compute the Hellinger divergence from an observed
    ## table
    v <- c(x, y)
    n1 <- sum(x, na.rm = TRUE)
    n2 <- sum(y, na.rm = TRUE)
    n <- cbind(n1, n2)
    ## pseudo-counts added if at least one of the cell counts is zero
    if (sum(v == 0) > 0) {
        p1 <- (x + 1)/(n1 + length(x))
        p2 <- (y + 1)/(n2 + length(x))
    } else {
        p1 <- x/n1
        p2 <- y/n2
    }
    p <- cbind(p1, p2)

    w <- (2 * n[1] * n[2]) / (n[1] + n[2])
    sum_hd <- sapply(seq_along(x),
                     function(k) (sqrt(p[k, 1]) - sqrt(p[k, 2]))^2)

    return(w * sum(sum_hd, na.rm = TRUE))
}
