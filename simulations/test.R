compute_factors <- function(y, no.de) {
    standard <- exp(colMeans(log1p(y)))
    clrm1 <- standard - 1
    best <- colSums(no.de)
    list(
        standard = standard / mean(standard),
        clrm1 = clrm1 / mean(clrm1),
        best = best / mean(best)
    )
}

plot_factors <- function(truth, factors, main) {
    plot(truth, factors, xlab="true bias", ylab="size factor", main=main, col="grey40", log="xy")
    abline(a=0, b=1, col="black", lty=2)
    mad <- median(abs(factors/truth - 1))
    legend("topleft", legend=paste0("MAD = ", signif(mad, 2)), bty="n")
}

unlink("results", recursive=TRUE)
dir.create("results", showWarnings=FALSE)

#########################################################
# Background-only, no differential abundance between any cells.

local({
    set.seed(7836490) # random-button-mash seed

    tag.means <- 1:50 / 20
    cell.biases <- 2^rnorm(1000)
    cell.biases <- cell.biases / mean(cell.biases)
    mu <- outer(tag.means, cell.biases)
    y <- matrix(rpois(length(mu), lambda=mu), ncol=length(cell.biases))

    factors <- compute_factors(y, y)

    png("results/bg_only.png", res=150, width=8, height=3, units="in")
    par(mfrow=c(1,3))
    plot_factors(cell.biases, factors$standard, main="Standard")
    plot_factors(cell.biases, factors$clrm1, main="CLRm1")
    plot_factors(cell.biases, factors$best, main="Optimal")
    dev.off()
})

#########################################################
# Each cell has one randomly chosen tag that is increased by 100-fold.

local({
    set.seed(310734893) # random-button-mash seed

    tag.means <- 1:50 / 20
    cell.biases <- 2^rnorm(1000)
    cell.biases <- cell.biases / mean(cell.biases)
    mu <- outer(tag.means, cell.biases)
    y <- matrix(rpois(length(mu), lambda=mu), ncol=length(cell.biases))

    # Injecting one highly abundant tag.
    y2 <- y
    for (i in seq_along(cell.biases)) {
        j <- sample(length(tag.means), 1L)
        y2[j,i] <- rpois(1, tag.means[j] * 100)
    }

    factors <- compute_factors(y2, y)

    png("results/one_tag.png", res=150, width=8, height=3, units="in")
    par(mfrow=c(1,3))
    plot_factors(cell.biases, factors$standard, main="Standard")
    plot_factors(cell.biases, factors$clrm1, main="CLRm1")
    plot_factors(cell.biases, factors$best, main="Optimal")
    dev.off()
})

#########################################################
# Each cell has zero to two tags that are increased by 100-fold.

local({
    set.seed(2367419) # random-button-mash seed

    tag.means <- 1:50 / 20
    cell.biases <- 2^rnorm(1000)
    cell.biases <- cell.biases / mean(cell.biases)
    mu <- outer(tag.means, cell.biases)
    y <- matrix(rpois(length(mu), lambda=mu), ncol=length(cell.biases))

    # Injecting 0-2 highly abundant tags.
    y2 <- y
    for (i in seq_along(cell.biases)) {
        chosen <- sample(0:2, 1L)
        j <- sample(length(tag.means), chosen)
        y2[j,i] <- rpois(chosen, tag.means[j] * 100)
    }

    factors <- compute_factors(y2, y)

    png("results/multi2_tag.png", res=150, width=8, height=3, units="in")
    par(mfrow=c(1,3))
    plot_factors(cell.biases, factors$standard, main="Standard")
    plot_factors(cell.biases, factors$clrm1, main="CLRm1")
    plot_factors(cell.biases, factors$best, main="Optimal")
    dev.off()
})

#########################################################
# Each cell has zero to 10 tags that are increased by 100-fold.

local({
    set.seed(2367419) # random-button-mash seed

    tag.means <- 1:50 / 20
    cell.biases <- 2^rnorm(1000)
    cell.biases <- cell.biases / mean(cell.biases)
    mu <- outer(tag.means, cell.biases)
    y <- matrix(rpois(length(mu), lambda=mu), ncol=length(cell.biases))

    # Injecting 0-10 highly abundant tags.
    y2 <- y
    for (i in seq_along(cell.biases)) {
        chosen <- sample(0:10, 1L)
        j <- sample(length(tag.means), chosen)
        y2[j,i] <- rpois(chosen, tag.means[j] * 100)
    }

    factors <- compute_factors(y2, y)

    png("results/multi10_tag.png", res=150, width=8, height=3, units="in")
    par(mfrow=c(1,3))
    plot_factors(cell.biases, factors$standard, main="Standard")
    plot_factors(cell.biases, factors$clrm1, main="CLRm1")
    plot_factors(cell.biases, factors$best, main="Optimal")
    dev.off()
})


#########################################################
# Highly variable background.

local({
    set.seed(2367419) # random-button-mash seed

    tag.means <- rep(c(0, 100), 20)
    cell.biases <- 2^rnorm(1000)
    cell.biases <- cell.biases / mean(cell.biases)
    mu <- outer(tag.means, cell.biases)
    y <- matrix(rpois(length(mu), lambda=mu), ncol=length(cell.biases))

    factors <- compute_factors(y, y)

    png("results/bg_variable.png", res=150, width=8, height=3, units="in")
    par(mfrow=c(1,3))
    plot_factors(cell.biases, factors$standard, main="Standard")
    plot_factors(cell.biases, factors$clrm1, main="CLRm1")
    plot_factors(cell.biases, factors$best, main="Optimal")
    dev.off()
})
