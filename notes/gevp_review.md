---
title: Review of `gevp` sorting
author: Martin Ueding (<ueding@hiskp.uni-bonn.de>)
date: 2019-02-21
colorlinks: true
urlcolor: blue
...

We have experienced a problem with the sorting of the eigenvalues in the GEVP
for $t < t_0$. In [hadron PR
157](https://github.com/HISKP-LQCD/hadron/pull/157) I was asked to review lines
93--121 of `R/gevp.R`. I will review this in commit
`2863b591666bb6399bd9a09bbdf22bc99a3ffebb`.

# Analytic problem statement

To quote from `man/gevp.Rd`:

> The generalised eigenvalue problem
> $$ C(t) v(t,t_0) = C(t_0) \lambda(t,t_0) v(t,t_0) $$
> is solved by performing a Cholesky decomposition of $C(t_0)=L^t L$ and
> transforming the GEVP into a standard eigenvalue problem for all values of
> $t$. The matrices $C$ are symmetrised for all $t$. So we solve for $\lambda$
> $$ (L^t)^{-1} C(t) L^{-1} w = \lambda w $$
> with $w = L v$ or the wanted $v = L^{-1} w$.
> 
> The amplitudes can be computed from
> $$ A_i^{(n)}(t) = \frac{\sum_{j}C_{ij}(t) v_j^{(n)}(t,t_0)}
> {\sqrt{(v^{(n)}, Cv^{(n)})(\exp(-mt)\pm \exp(-m(t-t)))}} $$
> and this is what the code returns up to the factor
> $$ 1/\sqrt{\exp(-mt)\pm \exp(-m(t-t))} $$
> The states are sorted by their eigenvalues when "values" is chosen. If
> "vectors" is chosen, we take
> $$ \max\big( \sum_i \langle v(t_0,i), v(t, j)\rangle\big) $$
> with $v$ the eigenvectors.

We want to consistently sort the eigenvalues for $t > t_0$.

# Actual implementation

## Sorting time

The `gevp` function has the parameter `sort.type`, which we will only use with
`"vectors"` here. And there is a logical `sort.t0`, whose action we will
investigate in detail.

Depending on `sort.t0`, an integer `t.sort` will be computed which holds the
time slice with which the scalar product is formed. That is not entirely true
though, as there is a fallback to $t_0 + 2$ later.

We have $t_\text{sort} = t_0 + 2$ by default. But if the user specified
`sort.t0 = FALSE`, the sorting does not happen with respect the fixed value of
$t_0 + 2$ but rather depends on the time slice $t$ that we are working with:
$$
t_\text{sort}(t) = \begin{cases}
t & \text{for } t > t_0 + 1 \\
t + 2 & \text{for } t < t_0 - 1 \\
t_0 + 2 & \text{otherwise} \,.
\end{cases}
$$
Just to be explicit: The special values are $t \in \{ t_0-1, t_0, t_0+1 \}$.
Applying these three cases to a range of time slices $t$, we get the following
$t_\text{sort}$:

| Case | $t$ | $t_\text{sort}$ |
| ---: | :--- | :--- |
| 2 | $t_0-3$ | $t_0 - 1$ |
| 2 | $t_0-2$ | $t_0$ |
| 3 | $t_0-1$ | $t_0 + 2$ |
| 3 | $t_0+1$ | $t_0 + 2$ |
| 1 | $t_0+2$ | $t_0 + 2$ |
| 1 | $t_0+3$ | $t_0 + 3$ |

This looks a bit strange, but perhaps it will make sense later on.

## Eigenvectors

The eigenvectors depend on $t$ and $t_0$. $t_0$ is arbitrary but fixed, so the
remaining dependency is on $t$. Each of them has length `matrix.size` because
it tells us how the various correlators in the matrix have to be combined to
give a pure state. There are `matrix.size` many eigenvectors on each time
slice.

The variable `evectors` holds them and is initialized with `NA` at first in
line 58:

```r
evectors <- array(NA, dim=c(Thalf + 1, matrix.size, matrix.size))
```

The loop over $t$ goes from $t_0 + 1$ to $T/2$ and then from $t_0 - 1$ to 0. In
line 95 the eigenvalues on time slice $t$ are computed:

```r
variational.solve <- eigen(t(invL) %*% cM %*% invL,
  symmetric = TRUE, only.values = FALSE, EISPACK = FALSE)
```

The `variational.solve` is a `list` that contains `values` and `vectors`.

Then in line 106 the scalar product between the eigenvectors on the current
time slice and the ones on $t_\text{sort}$ is computed:

```r
X <- abs(t(variational.solve$vectors) %*% evectors[t.sort, , ])
```

For the first iteration ($t = t_0 + 1$) the array `evectors` is completely filled
with `NA`. Therefore sorting by eigenvalues is performed for this case.

The `variational.solve$vectors` contain the eigenvectors $w$, but the
array `evectors` that is assigned slices in line 141 contains the eigenvectors
$v = L^{-1} w$. The documentation of `gevp` promises to form the scalar
product $\langle v(t_0, i), v(t, j) \rangle$, which mismatches the
implementation. I cannot see whether that causes a problem, at this point.

In line 141 the eigenvectors $w$ are sorted and multiplied with $L^{-1}$ to
give the eigenvalues $v$. These are then stored in `evectors` at slice $t + 1$.
The following is an overview over the slices that are created with the
iterations:

| $t$ | $t_\text{sort}$ | generated `evectors` |
| --- | --- | --- |
| $t_0 + 1$ | — | $t_0 + 2$ |
| $t_0 + 2$ | $t_0 + 2$ | $t_0 + 3$ |
| $t_0 + 3$ | $t_0 + 3$ | $t_0 + 4$ |
| … | … | … |
| $T/2$ | $T/2$ | $T/2 + 1$ |
| $t_0 - 1$ | $t_0 + 2$ | $t_0$ |
| $t_0 - 2$ | $t_0$ | $t_0 - 1$ |
| $t_0 - 3$ | $t_0 - 1$ | $t_0 - 2$ |
| … | … | … |
| 0 | 2 | 1 |

We can see that the the `evectors` slice that is needed for the next iterations
is always generated before it is used via $t_\text{sort}$. Also the first index
of `evectors` has dimension $T/2 + 1$, therefore the indices can be from the
interval $[1, T/2+1]$, which fits perfectly.

## Ordering

The eigenvalues on $t_0 + 1$ have been sorted simply by their values. This also
sorts the eigenvectors correspondingly. Now we want to sort the other
eigenvalues by matching the eigenvectors to the ones at $t_\text{sort}$. The
rational is that the mixing of the different operators has not changed much,
but the values of the correlators at the various time slices might have changed
a lot. The eigenvectors are therefore a better proxy for identifying the same
state across time slices.

As hinted earlier, the ordering for all $t$ except $t_0 + 1$ shall be computed
with the scalar products $X_{ij} = |\langle v(t, i), v(t_\text{sort},
j)\rangle |$:

```r
X <- abs(t(variational.solve$vectors) %*% evectors[t.sort, , ])
```

This is a square matrix of length `matrix.size`. We want to get an ordering (a
permutation of $1, \ldots, N$) for every $i$. First we have to decide whether
it shall be ascending or descending ordering. This is done in
lines 100f:[^ternary]

```r
decreasing <- TRUE
if (t < t0) decreasing <- FALSE
```

[^ternary]:

    I'd prefer to write this assignement as `decreasing <- t >= t0`, or if it
    has to be as like this using R's equivalent of a ternary operator:

    ```r
    decreasing <- if (t < t0) FALSE else TRUE
    ```

    This way the variable is assigned its value once and is constant. But that
    is a different story.

The important bit is that we do a decreasing sort for the time slices $t \ge
t_0$. But actually this value is ignored, at least most of the time. We'll see.

The first attempt to find an ordering is in line 108:

```r
idx <- apply(X, 1, order, decreasing = TRUE)
```

For each $i$ it gives the values $\{ \sigma(j) \}$ such that the overlaps
$X_{i,\sigma(j)}$ for fixed $i$ are sorted from largest to smallest. Also note
that the parameter `decreasing` is always set to `TRUE`, no matter what the
variable `decreasing` contains.

From all these orderings we just take the first one per $i$:

```r
sortindex <- idx[1, ]
```

This means that for each eigenvector at slice $t$ we pick the one eigenvector
at $t_\text{sort}$ which has the greatest overlap. This definition does not
exclude ambiguities, multiple vectors could have greatest overlap with a common
reference vector. For instance a very large eigenvector at $t_\text{sort}$
might have more overlap with an eigenvector at $t$ which points into a
different direction but is much smaller. Let us take these eigenvectors:
$$
v(t_\text{sort}, 1) = (100, 2)^\mathrm T
\,,\quad
v(t_\text{sort}, 2) = (0, 1)^\mathrm T
\,,\qquad
v(t, 1) = (0, 1)^\mathrm T
\,.
$$

From the looks of it, we would want to associate $v(t, 1)$ with
$v(t_\text{sort}, 2)$ because they are identical. The scalar product is
$X_{1,2} = 1$. But $v(t_\text{sort}, 1)$ has a much larger overlap, $X_{1,1} =
2$. This mean that we would associate it to the wrong eigenvector just because
the reference one is very large. Perhaps it would be better to consider instead
something corresponding to the angle in eigenspace:
$$
X_{ij} = \frac{|\langle v(t, i), v(t_\text{sort}, j)\rangle |}
{\|v(t, i)\| \, \|v(t_\text{sort}, j)\|} \,.
$$
With that definition my pathological example would be resolved better.

Still there could be ambiguities. Unless we are in `sort.t0` mode (which always
uses $t_\text{sort} = t_0 + 2$) already, we just try using that as a sorting
time:

```r
idx <- apply(abs(t(variational.solve$vectors) %*% evectors[t0 + 2, , ]),
             1, order, decreasing = TRUE)
sortindex <- idx[1, ]
```

There are a couple of time slices that already use $t_0 + 2$ as the sorting
time (see previous table), so that might just be redundant. We also always use
decreasing sort here.

If that fails again, we just sort by eigenvalues on that time slice, but this
time honoring the `decreasing` setting.

```r
sortindex <- order(variational.solve$values, decreasing = decreasing)
```

So does the direction of the sorting make sense? For plain eigenvalue sorting
the `decreasing` parameter is honored and it makes sense that the order of
eigenvalues is reversed when passing $t_0$. This arises because the eigenvalues
at $t_0$ are all normalized to unity, therefore the smallest mass has the
smallest slope and therefore the highest eigenvalues for $t > t_0$ and the
smallest eigenvalues for $t < t_0$.

For the eigenvectors one has to think a bit more thorough. We expect that the
correlators $C_{ij}(t)$ in the matrix are linear combinations of the compositions
of the spectrum $\lambda_i(t)$. In general we have
([0808.1017](https://arxiv.org/abs/0808.1017v1))
$$ C_{ij}(t) = \sum_{n = 1}^\infty \exp(-E_n t) \psi_{ni} \psi_{nj} \,, $$
where the $\psi$ are overlaps of the operators $i$ and $j$ with the state $n$.
We of course want the spectrum $E_n$. I would think that the mixing that
determines how the states $|n\rangle$ couple to the operators is similar on all
time slices. The choice of $t_0$ is done after they have been computed, so
there is no physical explanation for the eigenvectors to depend on $t_0$ and
therefore the directions in eigenspace for all eigenvectors should be the same.

If that is true, then using `decreasing = TRUE` always is exactly the right
thing to do.

# Conclusion

There are two things that I find suspicious:

-   The scalar products are formed between $w(t)$ and $v(t_\text{sort})$ and
    not between $v(t)$ and $v(t_\text{sort})$. Depending on $L^{-1}$ this could
    make them point into completely different directions.

-   The matching is done via the scalar product, which might lead to strange
    results when the magnitude of reference eigenvectors is vastly different.
    Perhaps using a normalized expression corresponding to an angle would be
    better.

The actual sorting direction seems consistent between $t > t_0$ and $t < t_0$
for both the eigenvalue and eigenvector method.
