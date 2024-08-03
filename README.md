# Normalization of ADT data

## Background

Most workflows for scaling normalization of ADT data use the geometric mean as the size factor, based on the CLR method used by Stoeckius et al. (2017).
This is a simple and pragmatic solution to the problem of composition biases introduced by a minority of high-abundance tags.

Consider a cell $i$ with $n$ tags where the count for tag $t$ is $`y_{it}`$.
Assume we have another cell $j$ with the same counts as $i$ except for one tag $t'$, where $`y_{jt'} = fy_{it'}`$ for $f \gg 1$.
If we use the total count as the size factor for each cell (i.e., $`\sum_t y_{it}`$),
the ratio of the size factors between $i$ and $j$ is a linear function of $f$;
this represents the composition bias introduced by the differential abundance in $l$.
For comparison purposes, let's consider the case where all $`y_{it}`$ are equal, such that the composition bias simplifies to $`1 + (f-1)n^{-1}`$.
In contrast, if we use the geometric mean (i.e., $`\sqrt[n]{\prod_t y_{it}}`$), the composition bias is $\sqrt[n]{f}$,
which is always smaller than the bias of the total count-derived size factors when $f > 1$.

An obvious issue with the geometric mean is that it is equal to zero when one or more values are zero.
As such, we usually add a pseudo-count - typically 1 - to ensure that some information is preserved from the non-zero counts.
(Alternatively, we could directly replace zeros with a value of 1.)
This workaround introduces its own bias in the form of a fold-change from the expected value of the tag with the zero count and its pseudo-count-based replacement,
effectively overestimating the size factor. 

## Improving performance at low counts

The "standard" CLR formula for the size factor for cell $i$ is $`\sqrt[n]{\prod_t (y_{it} + 1)}`$.
This is typically implemented as $`f(z_i) = \exp(z_i)`$ where $`z_i`$ is the mean of $`g(y_{it}) = \log(y_{it} + 1)`$,
as used by [**Seurat**](https://github.com/satijalab/seurat/blob/1549dcb3075eaeac01c925c4b4bb73c73450fc50/R/preprocessing5.R#L345)
and [**muon**](https://github.com/scverse/muon/blob/94917d23291f329a19b3c282276c960d414319ad/muon/_prot/preproc.py#L229).
However, attentive readers will notice that the addition of a pseudo-count means that $f()$ is not the inverse of $g()$.
Perhaps we should consider defining $`f(z_i) = \exp(z_i) - 1`$ for the sake of symmetry.
The size factor would then be $`\sqrt[n]{\prod_t (y_{it} + 1)} - 1`$, which we call the "CLRm1" size factor.

Despite its rather _ad hoc_ derivation, the CLRm1 approach works quite well.

INSERT GRAPH HERE.

Let's try to rationalize this behavior by considering two cells $k$ and $k'$ that only differ in their counts by some scaling factor $a$. 
The ideal normalization method would produce a size factor for $k'$ that is $a$-fold larger than that of $k$,
thus eliminating the scaling difference between the two cells.

- In the case where $`y_{kt}`$ is equal to some constant $`c_k`$ for all $t$, the CLRm1 size factor simplifies to $`(c_k + 1) - 1`$ for $k$ and $`(ac_k + 1) - 1`$ for $k'$.
  The ratio in the size factors will be equal to $a$.
- A generalization of the previous point involves approximating the geometric mean with the arithmetic mean.
  This approximation is satisfactory if the variance in $`y_{kt}`$ is low relative to the mean (see Equation 31 and related discussion in Rodin, 2014).
  Doing so simplifies the CLRm1 factor to $`n^{-1}\sum_t(y_{kt} + 1) - 1`$ for cell $k$ and $`n^{-1}\sum_t(ay_{kt} + 1) - 1`$ for cell $k'$, again yielding a ratio of $a$.
- If all $`y_{kt}`$ and $`ay_{kt}`$ are much greater than 1, the addition or subtraction of the pseudo-count can be ignored entirely.
  The size factors for the two cells cancel out, leaving us with $a$.
- In the rare case that all $`y_{kt}`$ are much less than 1, we can approximate $`\prod_t (1 + y_{kt}) \approx 1 + \sum_t y_{kt}`$.
  We can further approximate $`\sqrt[n]{1 + z} \approx 1 + zn^{-1}`$ when $z$ is close to zero.
  This allows us to obtain a size factor of $`(1 + n^{-1}\sum_t y_{kt}) - 1`$ for $k$ and $`(1 + an^{-1}\sum_t y_{kt}) - 1`$ for $k'$,
  which again cancels out to $a$.

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
