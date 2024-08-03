# Normalization of ADT data

## Background

Most workflows for scaling normalization of ADT data use the geometric mean as the size factor, based on the CLR method introduced by Stoeckius et al. (2017).
This is a simple and pragmatic solution to the problem of composition biases introduced by a minority of high-abundance tags.

Consider a cell $i$ with $n$ tags where the count for tag $t$ is $`y_{it}`$.
Assume we have another cell $j$ with identical counts except for one tag $t'$, where $`y_{jt'} = fy_{it'}`$ for $f >> 1$.
With library size scaling, the size factors will have a fold-difference of $1 + (f-1)s_{it'}$, where $s_{it'}$ is the ratio of $y_{it'}$ to the sum of $y_{it}$ for all $t$;
this represents the composition bias introduced by the differential abundance in $l$.
If we assume all counts are equal in cell $i$, the above expression simplifies to $1 + (f-1)n^{-1}$.
In contrast, the composition bias is only $\sqrt[n]{f}$ with geometric mean.

An obvious issue with the geometric mean is that it is equal to zero when one or more values are zero.
As such, we usually add a pseudo-count - typically 1 - to ensure that some information is preserved from the non-zero counts.
(Alternatively, we could directly replace zeros with a value of 1.)
This workaround introduces its own bias in the form of a fold-change from the expected value of the tag with the zero count and its pseudo-count-based replacement,
effectively overestimating the size factor. 

## Improving performance at low counts

The geometric mean for cell $i$ can be expressed as $`f(z_i) = \exp(z_i)`$ where $`z_i`$ is the mean of $`g(y_{it}) = \log(y_{it})`$. 
The addition of a pseudo-count is typically performed by replacing $`g(y_{it}) = \log(y_{it} + 1)`$,
which is the implementation used by [**Seurat**](https://github.com/satijalab/seurat/blob/1549dcb3075eaeac01c925c4b4bb73c73450fc50/R/preprocessing5.R#L345)
and [**muon**](https://github.com/scverse/muon/blob/94917d23291f329a19b3c282276c960d414319ad/muon/_prot/preproc.py#L229).
However, attentive readers will notice that $f()$ is no longer the inverse of our new $g()$.
Perhaps we should consider defining $`f(z_i) = \exp(z_i) - 1`$ for the sake of symmetry.

Surprisingly enough, this works quite well.

INSERT GRAPH HERE.

Let's try to rationalize this behavior by considering two cells $k$ and $k'$ that only differ in their counts by some scaling factor $a$. 

- In the case where $`y_{kt}`$ is equal to some constant $`c_k`$ is equal for all $t$, the size factor simplifies to $`(c_k + 1) - 1`$ for $k$ and $`(ac_k + 1) - 1`$ for $k'$.
  The ratio in the size factors will be equal to $a$, meaning that we correctly eliminate the scaling difference between the two cells.
- A generalization of the previous point involves approximating the geometric mean with the arithmetic mean.
  This approximation is satisfactory if the variance in $`y_{kt}`$ is low relative to the mean (see Equation 31 and related discussion of Rodin, 2014).
  Doing so yields a size factor of $`n^{-1}\sum_t(y_{kt} + 1) - 1`$ for cell $k$ and $`n{-1}\sum_t(ay_{kt} + 1) - 1`$ for cell $k'$
  again simplifying down to a relative difference of $a$.
- If all $`y_{kt}`$ and $`ay_{kt}`$ are much greater than 1, the addition or subtraction of the pseudo-count can be ignored entirely.
  The size factors for the two cells cancel out perfectly, leaving us with $a$.
- In the rare case that all $`y_{kt}`$ are much less than 1, we can approximate $`\prod_t (1 + y_{kt}) \approx 1 + \sum_t y_{kt}`$.
  We can further approximate $`\sqrt[n]{1 + z} \approx 1 + zn^{-1}`$ zhen $z$ is close to zero.
  This allows us to obtain a size factor of $`(1 + n^{-1}\sum_t(y_{kt})) - 1`$ for $k$ and $`(1 + an^{-1}\sum_t(y_{kt})) - 1`$ for $k'$,
  which again cancels out to $a$ between the two cells.

This analysis suggests that our approach will deteriorate when $`y_{kt}`$ is highly variable with at least one small/zero value.
Our hope is that this does not happen too frequently in real data, as there should not be large fluctuations in the ambient concentrations of different tags.
Of course, differentially abundant tags will also introduce variation in $`y_{kt}`$ but some loss of accuracy is to be expected from composition bias.

## References

Stoeckius M, Hafemeister C, Stephenson W, et al. (2017).
Simultaneous epitope and transcriptome measurement in single cells.
_Nature Methods_ 14, 865-868.

Rodin (2014).
Variance and the Inequality of Arithmetic and Geometric Means.
_arXiv_ doi:10.48550/arXiv.1409.0162.
