<hr>
<style>
body{
text-align: justify
}
</style>

# Introduction

This project studies the paper titled \`\`Linear Hypothesis Testing in
Linear Models With High-Dimensional Responses” by Li and Li (2022)
(Runze Li 2022) where the authors talk about the high-dimensional
(*p*≫*n*) scenarios where Hotelling *T*<sup>2</sup> cannot be used to do
one-sample and two-sample testing because of the singularity of sample
covariance matrix. Although, there have several papers on solving this
problem and its subsequent problems, the paper of our interest develops
general linear hypothesis testing frameworks using projection matrices,
which can be applied to one-sample, two-sample and multiple sample means
testing in high-dimensions. The main goal of this project is to
understand, simplify the concepts and lucidly explain this new
projection test and its theoretical properties through a blog. As we
proceed further in this project, we both shall elucidate the theorems in
the main paper as far as possible.

One sample and two sample tests have gained a lot of focus in the last
two decades. The Hotelling *T*<sup>2</sup> test for the one-sample and
two-sample mean problem is the uniformly most powerful test invariant
under affine transformations for fixed or low-dimensional data. The
singularity of the sample covariance matrix renders the Hotelling
*T*<sup>2</sup> test invalid when the dimension *p* of observations
exceeds the sample size *n*. Numerous tests for high-dimensional one and
two-sample mean testing have been proposed to address the singularity.
When examining the effects of dimensionality on two-sample mean testing,
Bai and Saranadasa (1996) (Bai and Saranadasa 1996) promoted
sum-of-squares-type statistics for the mean difference and disregarded
correlation matrices. Their test circumvents the issue brought on by the
sample covariance matrices’ rank deficiency. There are also various
other tests that have been created to address any potential issues that
may arise in this area. To create the methods, some of them used
supremum-type statistics, U-statistics, and so on. The authors of this
paper offer a projection test for a specific setting. The correlation
information is used by the projection test to increase the test’s power.
They build the projection test by first projecting the data into low
dimension and then performing the test on the projected data. They
demonstrated that there are low-dimensional projection matrices that
leverage covariance information to maximize the power of projection
tests. They further demonstrate that the dimension of the optimal
projections is not more than *m*, and they derive the *m* optimal
projection directions to maximize power. They adopt a sample-splitting
approach to perform projection tests, in which they use the first
portion of the sample for estimate of projection matrices and the second
part for testing. They also suggest U-projection tests to increase the
test’s power based on the sample-splitting technique by creating a
U-type test statistic from all samples. They have also established the
theoretical properties of the U-projection test.

# Notations

Let us consider the classical linear model setup:
*Y* = *X**B* + *E*
where *Y* is an *n* × *p* response matrix, *X* is an *n* × *d* design
matrix whose dimension *d* is fixed, *B* is a *d* × *p* regression
coefficient matrix, and *E* is an *n* × *p* error matrix with mean zero
and covariance matrix *I*<sub>*n*</sub> ⊗ *Σ* for a positive definite
matrix *Σ*. We denote
vec(*A*) = (*a*<sub>1</sub><sup>⊤</sup>,⋯,*a*<sub>*n*</sub><sup>⊤</sup>)<sup>⊤</sup>
for
*A* = (*a*<sub>1</sub>,*a*<sub>2</sub>,⋯,*a*<sub>*n*</sub>)<sup>⊤</sup>
from here onward. We are interested in testing the following hypotheses:
$$\begin{aligned}
    H_0&: A_0B = 0 \\
    H_1 & : A_0B\ne 0 \notag
\end{aligned}$$
where *A*<sub>0</sub> is an *m* × *d* known constant matrix with both
*m* and *d* fixed and *A*<sub>0</sub> is of full row rank which implies
*m* ≤ *d*. Under usual assumptions of dimensions, i.e., *n* \< *p*, the
likelihood ratio test (LRT) and its asymptotic convergence to
chi-squared distribution holds good. However, in high dimensions when
*p* \>  \> *n*, the usual properties of LRT fails. In this blog, we
elucidate a methods of testing hypotheses in such scenarios using
projection matrices that project the high-dimension problem to low
dimensions such that a tractable solution is achieved using the
conventional LRT.

# Optimal Projection Direction

Firstly, let us redefine the problem. Let *P* be a *p* × *r* full column
rank projection matrix with *r* ≪ *p* and *r* \< *n*, consider the
*P*-projected model:

*Y*<sup>⋆</sup> = *X**B*<sup>⋆</sup> + *E*<sup>⋆</sup>

where *Y*<sup>⋆</sup> = *Y**P*, *B*<sup>⋆</sup> = *B**P*, and
*E*<sup>⋆</sup> = *E**P*, which is the *n* × *r* projected error matrix
with mean zero and cov (vec(*E*<sup>⋆</sup>)) = *I*<sub>*n*</sub>⊗
(*P*<sup>⊤</sup>*Σ**P*)

The corresponding projected hypothesis becomes

*H*<sub>0*P*</sub> : *A*<sub>0</sub>*B*<sup>⋆</sup> = 0   versus   *H*<sub>1*P*</sub> : *A*<sub>0</sub>*B*<sup>⋆</sup> ≠ 0

and the projected test statistic is

$$\Lambda=\frac{\left\|G\_{P}\right\|}{\left\|G\_{P}+H\_{P}\right\|}$$

where *G*<sub>*P*</sub> is the residual sum of squares under
*H*<sub>1*P*</sub>, *G*<sub>*P*</sub> + *H*<sub>*P*</sub> is the
residual sum of squares under *H*<sub>0*P*</sub>.
*G*<sub>*P*</sub> = *P*<sup>⊤</sup>*Y*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>1</sub></sub>)*Y**P*,
*P*<sub>*H*<sub>1</sub></sub> = *X*(*X*<sup>⊤</sup>*X*)<sup>−1</sup>*X*<sup>⊤</sup>
and
*G*<sub>*P*</sub> + *H*<sub>*P*</sub> = *P*<sup>⊤</sup>*Y*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>0</sub></sub>)*Y**P*,
*P*<sub>*H*<sub>0</sub></sub> = *X**A*<sub>1</sub>(*A*<sub>1</sub><sup>⊤</sup>*X*<sup>⊤</sup>*X**A*<sub>1</sub>)<sup>−1</sup>*A*<sub>1</sub><sup>⊤</sup>*X*<sup>⊤</sup>,
and *A*<sub>1</sub> is a *d* × (*d*−*m*) matrix defined by
*A*<sub>1</sub> = (*A*<sub>0</sub><sup>⊤</sup>)<sup>⊥</sup>, which means
that (*A*<sub>0</sub><sup>⊤</sup>,*A*<sub>1</sub>) forms an orthogonal
matrix. The test is the LRT under normality assumption and is equivalent
to many useful test in various cases, like Hotelling *T*<sup>2</sup> in
the one-sample mean testing.

The one-sample and two-sample mean tests can be written as special
instances of the mentioned hypotheses test if *X* and *A*<sub>0</sub>
are appropriately specified. As a result, the approach suggested in this
report applies directly to the one-sample, two-sample mean and even
multiple-sample mean problems. The projected null hypothesis
*H*<sub>0*P*</sub> is rejected for small *Λ*, and *H*<sub>0</sub> is
rejected if *H*<sub>0*P*</sub> is rejected. In general, *H*<sub>0</sub>
and *H*<sub>0*P*</sub> are not equivalent. In this article, we will show
that there exists an optimal projection direction *P* which makes
*H*<sub>0*P*</sub> equivalent to *H*<sub>0</sub> and also maximizes the
power of the test.

The problem how to construct an optimal direction can be divided into
two sub-problems:

-   to find the dimension *r* of the optimal projection direction *P*,
    and

-   to find the optimal *P* of a particular dimension.

These issues are addressed in Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>.
We assume the following conditions on the multivariate linear model (1)
to derive the asymptotically optimal projection direction matrix.

**Condition 1**:
∃*C*<sub>1</sub> \> 0 ∋ ∥*X**X*<sup>⊤</sup>∥<sub>∞</sub> \< *C*<sub>1</sub>,
where
∥*A*∥<sub>∞</sub> = max<sub>*i*, *j*</sub>\|*a*<sub>*i*, *j*</sub>\| for
*A* = (*a*<sub>*i*, *j*</sub>)<sub>*i*, *j*</sub>.

**Condition 2**: The limit
lim<sub>*n* → ∞</sub>*n*<sup>−1</sup>(*X*<sup>⊤</sup>*X*) = *M*<sub>*X*</sub>
exists, where *M*<sub>*X*</sub> is a nonsingular *d* × *d* matrix.

**Condition 3**:
∃*C*<sub>2</sub> \> 0 ∋ 𝔼(*E*<sub>*i**j*</sub><sup>4</sup>) \< *C*<sub>2</sub>,
*i* = 1, …, *n*, *j* = 1, …, *p* for elements *E*<sub>*i*, *j*</sub> in
the error matrix.

Conditions 1-3 are quite mild and are used to guarantee asymptotic
normality of least squares estimate in linear models. For example,
Conditions 1 and 2 hold in the one-sample mean testing problem
automatically, and they also hold in the two-sample mean testing problem
if $\frac{n\_{1}}{n\_{1}+n\_{2}} \rightarrow \kappa \in(0,1)$, where
*n*<sub>1</sub> and *n*<sub>2</sub> are sample sizes of the first and
the second samples, respectively.

**Theorem 1** (Optimal Projection Direction). *Suppose that the error
matrix *E* follows a matrix normal distribution. Under Conditions 1-3,
the following statements are valid.*

*(A1) Let
*W* = *Σ*<sup>−1/2</sup>*B*<sup>⊤</sup>*X*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>0</sub></sub>)*X**B**Σ*<sup>−1/2</sup>,
where
*P*<sub>*H*<sub>0</sub></sub> = *X**A*<sub>1</sub>(*A*<sub>1</sub><sup>⊤</sup>*X*<sup>⊤</sup>*X**A*<sub>1</sub>)<sup>−1</sup>*A*<sub>1</sub><sup>⊤</sup>*X*<sup>⊤</sup>,
and *A*<sub>1</sub> is a *d*× (*d*−*m*) matrix which is defined by
*A*<sub>1</sub> = (*A*<sub>0</sub><sup>⊤</sup>)<sup>⊥</sup>. *W* is a
*p* × *p* matrix of rank *m*. Suppose *W* admits the eigenvalue
decomposition $W=\sum\_{i=1}^{m} \lambda\_{0 i} R\_{i} R\_{i}^\top$,
where 0 \< *λ*<sub>0*m*</sub> ≤ ⋯ ≤ *λ*<sub>01</sub> are *m* nonzero
eigenvalues, and *R*<sub>*i*</sub>, *i* = 1, …, *m*, are the
corresponding orthogonal eigenvectors. Then for any given *r* ≤ *m*,
*P*<sub>0</sub> = *Σ*<sup>−1/2</sup>(*R*<sub>1</sub>,…,*R*<sub>*r*</sub>)
is the projection direction matrix which is optimal among all *p* × *r*
matrices for the hypothesis testing problem (2).*

*(A2) The optimal dimension of projection matrix is less than or equal
to *m* for the hypothesis testing problem (2).*

*(A3) Set
*r* = *m*, *Σ*<sup>−1</sup>*B*<sup>⊤</sup>*A*<sub>0</sub><sup>⊤</sup> is
the projection direction matrix which is optimal among all *p* × *m*
matrices for the hypothesis testing problem (2).*

*(B) Without normality assumption, the statements in (A1)-(A3) are still
valid in asymptotic sense by changing optimal projection to
asymptotically optimal projection which maximizes the asymptotic local
power. *

*Proof.* To prove the theorem, we use the unexplained variability or
variability due to random error components. For the ordinary linear
model, orthogonal split of sum of squares is given by
*Y*<sup>2</sup> = *Y*<sup>*T*</sup>*P*<sub>*H*</sub>*Y* + *Y*<sup>*T*</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*</sub>)*Y*
where the later term in the summation is accountable for the sum of
squares due to errors.

We first prove Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>(A).
For linear model and under *H*<sub>1</sub>, the estimated *Y* is the
projection of *Y* onto the column space of *X*, denoted by
*P*<sub>*H*<sub>1</sub></sub>*Y*; and under *H*<sub>0</sub>, the
statement that *A*<sub>0</sub>*B* = 0 is equivalent to that *B* is of
the form *A*<sub>1</sub><sup>⊤</sup>*D*, where
(*A*<sub>0</sub>,*A*<sub>1</sub>) forms an orthogonal basis and *D* is a
(*d*−*m*) × *p* matrix. Hence the estimated *Y* under *H*<sub>0</sub> is
the projection of *Y* onto the column space of
*X**A*<sub>1</sub><sup>⊤</sup>, denoted by
*P*<sub>*H*<sub>0</sub></sub>*Y*. So
*G* = (*Y*−*P*<sub>*H*<sub>1</sub></sub>*Y*)<sup>⊤</sup>(*Y*−*P*<sub>*H*<sub>1</sub></sub>*Y*)=
*Y*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>1</sub></sub>)*Y*, *G* + *H* = (*Y*−*P*<sub>*H*<sub>0</sub></sub>*Y*)<sup>⊤</sup>(*Y*−*P*<sub>*H*<sub>0</sub></sub>*Y*) = *Y*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>0</sub></sub>)*Y*, *H* = *Y*<sup>⊤</sup>(*P*<sub>*H*<sub>1</sub></sub>−*P*<sub>*H*<sub>0</sub></sub>)*Y*.
Under normality assumption, LRT statistic for *H*<sub>0</sub> defined in
(1) is $\frac{\|G\|}{\|G+H\|}$.

Moreover, in the projection test, *Y* is replaced by *Y**P*. The LRT
statistic in the projection test can be written as
$\frac{\left\|G\_{P}\right\|}{\left\|G\_{P}+H\_{P}\right\|}$, where
*G*<sub>*P*</sub> is the residual sum of squares under
*H*<sub>1*P*</sub>, and *G*<sub>*P*</sub> + *H*<sub>*P*</sub> is the
residual sum of squares under *H*<sub>0*P*</sub>. Without loss of
generality, we assume that *Σ*<sup>1/2</sup>*P* is a projection matrix
*Q*, i.e., *Q* is symmetric and idempotent, so
*P* = *Σ*<sup>−1/2</sup>*Q* and

*G*<sub>*P*</sub> = *P*<sup>⊤</sup>*Y*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>1</sub></sub>)*Y**P* = *Q*<sup>⊤</sup>*Σ*<sup>−1/2</sup>*Y*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>1</sub></sub>)*Y**Σ*<sup>−1/2</sup>*Q*

Since
(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>1</sub></sub>)*Y* ∼ *N*(0,(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>1</sub></sub>)⊗*Σ*),
we have
*Σ*<sup>−1/2</sup>*Y*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>1</sub></sub>)*Y* ∼ *W*<sub>*p*</sub>(*n*−*d*,*I*<sub>*p*</sub>).
Furthermore, since *Q* is a *p* × *k* projection matrix,
*G*<sub>*P*</sub> = *Q*<sup>⊤</sup>*Σ*<sup>−1/2</sup>*Y*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>1</sub></sub>)*Y**Σ*<sup>−1/2</sup>*Q*∼
*W*<sub>*k*</sub>(*n*−*d*,*I*<sub>*k*</sub>). We have
*H*<sub>*P*</sub> = *P*<sup>⊤</sup>*Y*<sup>⊤</sup>(*P*<sub>*H*<sub>1</sub></sub>−*P*<sub>*H*<sub>0</sub></sub>)*Y**P* = *Q*<sup>⊤</sup>*Σ*<sup>−1/2</sup>*Y*<sup>⊤</sup>(*P*<sub>*H*<sub>1</sub></sub>−*P*<sub>*H*<sub>0</sub></sub>)*Y**Σ*<sup>−1/2</sup>*Q*.
Since
(*P*<sub>*H*<sub>1</sub></sub>−*P*<sub>*H*<sub>0</sub></sub>)*Y* ∼ *N*((*I*<sub>*n*</sub>−*P*<sub>*H*<sub>0</sub></sub>)*X**B*,(*P*<sub>*H*<sub>1</sub></sub>−*P*<sub>*H*<sub>0</sub></sub>)⊗*Σ*),
we have

*Σ*<sup>1/2</sup>*Y*<sup>⊤</sup>(*P*<sub>*I*<sub>1</sub></sub>−*P*<sub>*I*<sub>0</sub></sub>)*Y**Σ*<sup>−1/2</sup> ∼ *W*<sub>*p*</sub>(*m*,*I*<sub>*p*</sub>,*Σ*<sup>−1/2</sup>*B*<sup>⊤</sup>*X*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*I*<sub>0</sub></sub>)*X**B**Σ*<sup>−1/2</sup>)

Hence we have

$$\begin{aligned}
H\_{P} & =Q^\top \Sigma^{-1 / 2} Y^\top\left(P\_{H\_{1}}-P\_{H\_{0}}\right) Y \Sigma^{-1 / 2} Q \\
& \sim W\_{k}\left(m, I\_{k}, Q^\top \Sigma^{-1 / 2} B^\top X^\top\left(I\_{n}-P\_{H\_{0}}\right) X B \Sigma^{-1 / 2} Q\right) .
\end{aligned}$$

By Cochran’s Theorem, we know that *G*<sub>*P*</sub> is independent to
*H*<sub>*P*</sub> given *X* and *P*. Proof of Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>(A1).
Here we prove a more general Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>(A1)
for any *r* ≤ *p* by defining *R*<sub>*i*</sub>, *i* = *m* + 1, ⋯, *p*
as orthogonal eigenvectors corresponding to eigenvalue 0 . We first show
that *P*<sub>0</sub> = *Σ*<sup>−1/2</sup>*Q*<sub>0</sub> is the optimal
*p* × *r* projection direction matrix under normality assumption, where
*Q*<sub>0</sub> = (*R*<sub>1</sub>,⋯,*R*<sub>*r*</sub>) and
*R*<sub>*i*</sub> is defined as in Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a> .

Given *r* ≤ *m* fixed, it is easy to see that
$\lambda\_{P}=\frac{\left\|G\_{P}\right\|}{\left\|G\_{P}+H\_{P}\right\|}$
will be of the same distribution under *H*<sub>0</sub> no matter what
*Q* is, so the critical value will be the same for any *p* × *r*
projection matrix *Q*. So to compare the power of *p* × *r* projection
tests, we only need to compare
$P\left(\frac{\left\|G\_{P}\right\|}{\left\|G\_{P}+H\_{P}\right\|}\<c\right)$
for any positive *c*.

For any *p* × *r* projection matrix *Q*, it is easy to know that
*λ*<sub>*i*</sub> ≤ *λ*<sub>0*i*</sub>, where
*λ*<sub>*r*</sub> ≤ *λ*<sub>*r* − 1</sub>≤ ⋯ ≤ *λ*<sub>1</sub> are the
nonzero eigenvalues of
*Q*<sup>⊤</sup>*Σ*<sup>−1/2</sup>*B*<sup>⊤</sup>*X*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>0</sub></sub>)*X**B**Σ*<sup>−1/2</sup>*Q*,
and *λ*<sub>0*r*</sub>≤ *λ*<sub>0(*r*−1)</sub> ≤ ⋯ ≤ *λ*<sub>01</sub>
are the *r* largest eigenvalues of
*Σ*<sup>−1/2</sup>*B*<sup>⊤</sup>*X*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>0</sub></sub>)*X**B**Σ*<sup>−1/2</sup>.
Let
*Q*<sup>⊤</sup>*Σ*<sup>−1/2</sup>*B*<sup>⊤</sup>*X*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>0</sub></sub>)*X**B**Σ*<sup>−1/2</sup>*Q* = *U*<sup>⊤</sup>diag (*λ*<sub>1</sub>,…,*λ*<sub>*r*</sub>)*U*,
where the right-hand side is the eigenvalue decomposition of the left
and *U* is an *r* × *r* orthogonal matrix. We have

$$\begin{aligned}
        Q^\top \Sigma^{-1 / 2} B^\top X^\top\left(I\_{n}-P\_{H\_{0}}\right) X B \Sigma^{-1 / 2} Q & =U^\top \operatorname{diag}\left(\lambda\_{1}, \ldots, \lambda\_{r}\right) U \\
        & \leq U^\top \operatorname{diag}\left(\lambda\_{01}, \ldots, \lambda\_{0 r}\right) U \\
        & =U^\top Q\_{0}^\top \Sigma^{-1 / 2} B^\top X^\top\left(I\_{n}-P\_{H\_{0}}\right) X B \Sigma^{-1 / 2} Q\_{0} U .
    \end{aligned}$$

So

$$\begin{aligned}
        M= & U^\top Q\_{0}^\top \Sigma^{-1 / 2} B^\top X^\top\left(I\_{n}-P\_{H\_{0}}\right) X B \Sigma^{-1 / 2} Q\_{0} U \\
        & -Q^\top \Sigma^{-1 / 2} B^\top X^\top\left(I\_{n}-P\_{H\_{0}}\right) X B \Sigma^{-1 / 2} Q
    \end{aligned}$$

is a semi-definite matrix. Generate
*K* ∼ *W*<sub>*k*</sub>(0,*I*<sub>*r*</sub>,*M*), which is a *r* × *r*
random matrix of non-central Wishart distribution with zero degree of
freedom and non-central parameter *M* independent to (*X*,*Y*). Note
that *K* is semi-definite almost surely. From the property of
non-central Wishart distribution from Muirhead (2009) (Muirhead 2009),
we know that

$$H\_{\Sigma^{-1 / 2} Q}+K \stackrel{d}{=} H\_{\Sigma^{-1 / 2} Q\_{0} U} \stackrel{d}{=} H\_{\Sigma^{-1 / 2} Q\_{0}}$$

Since
$G\_{\Sigma^{-1 / 2} Q} \stackrel{d}{=} G\_{\Sigma^{-1 / 2} Q\_{0}}, G\_{\Sigma^{-1 / 2} Q}$
is independent to *H*<sub>*Σ*<sup>−1/2</sup>*Q*</sub> + *K* and
*G*<sub>*Σ*<sup>−1/2</sup>*Q*<sub>0</sub></sub> is independent to
*H*<sub>*Σ*<sup>−1/2</sup>*Q*<sub>0</sub></sub>, we have

$$\frac{\left\|G\_{\Sigma^{-1 / 2} Q}\right\|}{\left\|G\_{\Sigma^{-1 / 2} Q}+H\_{\Sigma^{-1 / 2} Q}+K\right\|} \stackrel{d}{=} \frac{\left\|G\_{\Sigma^{-1 / 2} Q\_{0}}\right\|}{\left\|G\_{\Sigma^{-1 / 2} Q\_{0}}+H\_{\Sigma^{-1 / 2} Q\_{0}}\right\|}$$

Thus, we have

$$\begin{aligned}
    P\left(\lambda\_{\Sigma^{-1 / 2} Q}\<c\right) & =P\left(\frac{\mid G\_{\Sigma^{-1 / 2} Q}}{\left\|G\_{\Sigma^{-1 / 2} Q}+H\_{\Sigma^{-1 / 2} Q}\right\|}\<c\right) \\
        & \leq P\left(\frac{\mid G\_{\Sigma^{-1 / 2} Q}}{\left\|G\_{\Sigma^{-1 / 2} Q}+H\_{\Sigma^{-1 / 2} Q}+K\right\|}\<c\right) \\
        & =P\left(\frac{\mid G\_{\Sigma^{-1 / 2} Q\_{0}}}{\left\|G\_{\Sigma^{-1 / 2} Q}+H\_{\Sigma^{-1 / 2} Q\_{0}}\right\|}\<c\right)=P\left(\lambda\_{\Sigma^{-1 / 2} Q\_{0}}\<c\right) .
\end{aligned}$$

Thus, *Σ*<sup>−1/2</sup>*Q*<sub>0</sub> is the optimal *p* × *r*
projection matrix.

Proof of Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>(A2)
For any *p* × *r* projection matrix *Q* and *r* \> *m* and from Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>(A2),
*P*<sub>0</sub> = *Σ*<sup>−1/2</sup>*Q*<sub>0</sub> is the optimal
*p* × *r* projection direction matrix, where *Q*<sub>0</sub>=
(*R*<sub>1</sub>,⋯,*R*<sub>*r*</sub>) and *R*<sub>*i*</sub> is defined
as in Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a> .
Note that
$P\_{0}=\left(\begin{array}{ll}P\_{01} & P\_{02}\end{array}\right)$,
where
*P*<sub>01</sub> = *Σ*<sup>−1/2</sup>(*R*<sub>1</sub>,⋯,*R*<sub>*m*</sub>)
is the optimal *p* × *m* projection direction matrix and
*P*<sub>02</sub>=
*Σ*<sup>−1/2</sup>(*R*<sub>*m* + 1</sub>,⋯,*R*<sub>*r*</sub>).

The LRT statistic for the *P*<sub>0</sub> projection test can be
decomposed correspondingly as follows:
$$\dfrac{\|G\_{P_0}\|}{\|G\_{P_0}+H\_{P_0}\|}=\dfrac{\|G\_{22}\|}{\|G\_{22}+H\_{22}\|}\dfrac{\left\|\begin{pmatrix}
    G\_{11} & G\_{12}\\
   G\_{21} & G\_{22}
\end{pmatrix}\right\|/\|G\_{22}\|}{\left\|\left(\begin{array}{cc}
    G\_{11}+H\_{11} & G\_{12}+H\_{12} \\
    G\_{21}+H\_{21} & G\_{22}+H\_{22}
    \end{array}\right)\right\| /\left\|G\_{22}+H\_{22}\right\|}$$
where *G*<sub>11</sub>, *G*<sub>11</sub> + *H*<sub>11</sub> are of
dimension *m* × *m*, and
*G*<sub>22</sub>, *G*<sub>22</sub> + *H*<sub>22</sub> are of dimension
(*r*−*m*)× (*r*−*m*).*G*<sub>22</sub> and
*G*<sub>22</sub> + *H*<sub>22</sub> can be seen as the residual sums of
squares under *H*<sub>1</sub> and *H*<sub>0</sub> for *P*<sub>02</sub>
projection test. It is easy to show that
(*G*<sub>22</sub>,*H*<sub>22</sub>) has the same distribution under
*H*<sub>1</sub> and *H*<sub>0</sub>. In fact, under *H*<sub>1</sub> and
*H*<sub>0</sub>, *G*<sub>22</sub> ∼ *W*<sub>*r* − *m*</sub>(*n*−*d*,*I*<sub>*r* − *m*</sub>), *H*<sub>22</sub> ∼ *W*<sub>*r* − *m*</sub>(*m*,*I*<sub>*m*</sub>)
and they are independent to each other. So
\|*G*<sub>22</sub>\|/\|*G*<sub>22</sub>+*H*<sub>22</sub>\| has the same
distribution under *H*<sub>1</sub> and *H*<sub>0</sub>. Moreover,
$\left\|\left(\begin{array}{ll}G\_{11} & G\_{12} \\ G\_{21} & G\_{22}\end{array}\right)\right\| /\left\|G\_{22}\right\|=\left\|G\_{11}-G\_{12} G\_{22}^{-1} G\_{21}\right\|$
can be seen as the generalized residual sum of squares of
*Y**P*<sub>01</sub> on *Y**P*<sub>02</sub> and *X*.
*W*<sub>1</sub> = *G*<sub>11</sub> − *G*<sub>12</sub>*G*<sub>22</sub><sup>−1</sup>*G*<sub>21</sub> ∼ *W*<sub>*m*</sub>(*n*−(*r*−*m*)−*d*,*I*<sub>*m*</sub>).

$$\left\|\left(\begin{array}{cc}
    G\_{11}+H\_{11} & G\_{12}+H\_{12} \\
    G\_{21}+H\_{21} & G\_{22}+H\_{22}
    \end{array}\right)\right\| /\left\|G\_{22}+H\_{22}\right\|=\mid G\_{11}+H\_{11}-\left(G\_{12}+H\_{12}\right)\left(G\_{22}+H\_{22}\right)^{-1}\left(G\_{21}+H\_{21}\right)$$

can be seen as the generalized residual sum of squares of
*Y**P*<sub>01</sub> on *Y**P*<sub>02</sub> and
*X**A*<sub>1</sub><sup>⊤</sup>. *W*<sub>0</sub>=
*G*<sub>11</sub> + *H*<sub>11</sub> − (*G*<sub>12</sub>+*H*<sub>12</sub>)(*G*<sub>22</sub>+*H*<sub>22</sub>)<sup>−1</sup>(*G*<sub>21</sub>+*H*<sub>21</sub>) = *W*<sub>1</sub> + *W*<sub>2</sub>,
where
*W*<sub>2</sub> ∼ *W*<sub>*m*</sub>(*m*,*I*<sub>*m*</sub>,*P*<sub>01</sub><sup>⊤</sup>*B*<sup>⊤</sup>*X*<sup>⊤</sup>(*I*<sub>*n*</sub>−
*P*<sub>*H*<sub>0</sub></sub>)*X**B**Σ*<sup>−1/2</sup>) and
*W*<sub>2</sub> is independent to *W*<sub>1</sub>.

In sum,

$$\frac{\left\|\left(\begin{array}{ll}
    G\_{11} & G\_{12} \\
    G\_{21} & G\_{22}
    \end{array}\right)\right\| /\left\|G\_{22}\right\|}{\left\|\left(\begin{array}{ll}
    G\_{11}+H\_{11} & G\_{12}+H\_{12} \\
    G\_{21}+H\_{21} & G\_{22}+H\_{22}
    \end{array}\right)\right\| /\left\|G\_{22}+H\_{22}\right\|} \stackrel{d}{=} \frac{\left\|W\_{1}\right\|}{\left\|W\_{1}+W\_{2}\right\|},$$

where
*W*<sub>1</sub> ∼ *W*<sub>*m*</sub>(*n*−(*r*−*m*)−*d*,*I*<sub>*m*</sub>), *W*<sub>2</sub> ∼ *W*<sub>*m*</sub>(*m*,*I*<sub>*m*</sub>,*P*<sub>01</sub><sup>⊤</sup>*B*<sup>⊤</sup>*X*<sup>⊤</sup>(*I*<sub>*n*</sub>−*P*<sub>*H*<sub>0</sub></sub>)*X**B**P*<sub>01</sub>),
and *W*<sub>1</sub> is dependent to *W*<sub>2</sub>.

Note that the distribution of
$\frac{\left\|W\_{1}\right\|}{\left\|W\_{1}+W\_{2}\right\|}$ does not
depend on *Y**P*<sub>02</sub>. So it is independent to
$\frac{\left\|G\_{22}\right\|}{\left\|G\_{22}+H\_{22}\right\|}$. The
proof here is similar to the proof in Lemma 8.4.3, Lemma 8.4.4 in
Anderson (2003) (Anderson 2003), which decomposes the test statistic
into the product of a series of independent Beta distributed variables.

Since $\frac{\left\|G\_{22}\right\|}{\left\|G\_{22}+H\_{22}\right\|}$
has the same distribution under *H*<sub>0</sub> and *H*<sub>1</sub>, the
test statistic

$$\frac{\left\|\left(\begin{array}{cc}
    G\_{11} & G\_{12} \\
    G\_{21} & G\_{22}
    \end{array}\right)\right\| /\left\|G\_{22}\right\|}{\left\|\left(\begin{array}{cc}
    G\_{11}+H\_{11} & G\_{12}+H\_{12} \\
    G\_{21}+H\_{21} & G\_{22}+H\_{22}
    \end{array}\right)\right\| /\left\|G\_{22}+H\_{22}\right\|}$$

is better than

$$\frac{\left\|G\_{22}\right\|}{\left\|G\_{22}+H\_{22}\right\|} \frac{\left\|\left(\begin{array}{ll}
    G\_{11} & G\_{12} \\
    G\_{21} & G\_{22}
    \end{array}\right)\right\| /\left\|G\_{22}\right\|}{\left\|\left(\begin{array}{ll}
    G\_{11}+H\_{11} & G\_{12}+H\_{12} \\
    G\_{21}+H\_{21} & G\_{22}+H\_{22}
    \end{array}\right)\right\| /\left\|G\_{22}+H\_{22}\right\|} .$$

Moreover,
$\frac{\left\|\left(\begin{array}{ll}G\_{11} & G\_{12} \\ G\_{21} & G\_{22}\end{array}\right)\right\| /\left\|G\_{22}\right\|}{\left\|\left(\begin{array}{ll}G\_{11}+H\_{11} & G\_{12}+H\_{12} \\ G\_{21}+H\_{21} & G\_{22}+H\_{22}\end{array}\right)\right\| /\left\|G\_{22}+H\_{22}\right\|}$
has the same distribution with *Λ*<sub>*P*<sub>01</sub></sub> only with
sample size reduced from *n* to *n* − (*r*−*m*).

Thus, *P*<sub>01</sub> projection test is powerful than
$\left(\begin{array}{ll}P\_{01} & P\_{02}\end{array}\right)$ projection
test. So the optimal *m* dimension projection test is better than the
optimal *r* dimension projection test with *r* \> *m*, which means that
the dimension for the optimal projection direction is less than or equal
to *m* under normality assumption.

Proof of Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>(A3).
Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>(A3)
is just a special case of Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>(A1)
and can be proved using the same techniques.

Proof of Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>
(B). It is sufficient to prove that the test statistic
$\Lambda=\frac{\left\|G\_{P}\right\|}{\left\|G\_{P}+H\_{P}\right\|}$
follows the same asymptotic distribution in both Gaussian and
non-Gaussian multivariate linear models, thus the asymptotic results in
the Gaussian case also hold in the non-Gaussian cases.

From Conditions 1, 2, and 3 of Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>,
we have $\frac{G\_{P}}{n} \stackrel{p}{\rightarrow} P^\top \Sigma P$
with convergence rate *O*<sub>*p*</sub>(*n*<sup>−1/2</sup>). So with the
delta method, we have

$$\begin{aligned}
        -\log (\Lambda) & =-\log \frac{\left\|G\_{P}\right\|}{\left\|G\_{P}+H\_{P}\right\|}=\log \left(\left\|I\_{k}+H\_{P} G\_{P}^{-1}\right\|\right) \\
        & =\log \operatorname{det}\left(I\_{k}+H\_{P}\left(P^\top \Sigma P\right)^{-1} / n\right)\left(1+O\_{p}\left(n^{-1 / 2}\right)\right)
    \end{aligned}$$

Hence the theorem holds if *H*<sub>*P*</sub> follows the same asymptotic
distribution under both the Gaussian and non-Gaussian cases.

For any projection direction matrix *P* of dimension *p* × *k*, we have
*B̂*<sub>*I*<sub>1</sub></sub><sup>\*</sup> = *B̂*<sup>\*</sup>=
(*X*<sup>⊤</sup>*X*)<sup>−1</sup>*X*<sup>⊤</sup>*Y**P*, and
*B̂*<sub>*H*<sub>0</sub></sub><sup>\*</sup> = *A*<sub>1</sub>(*A*<sub>1</sub><sup>⊤</sup>*X*<sup>⊤</sup>*X**A*<sub>1</sub>)<sup>−1</sup>*A*<sub>1</sub><sup>⊤</sup>*X*<sup>⊤</sup>*Y**P*,
where *A*<sub>1</sub> = (*A*<sub>0</sub><sup>⊤</sup>)<sup>⊥</sup>.
Hence,

$$d\_{B}=\hat{B\_{H\_{1}}^\*}-\hat{B\_{H\_{0}}^\*}=\left(\left(X^\top X\right)^{-1}-A\_{1}\left(A\_{1}^\top X^\top X A\_{1}\right)^{-1} A\_{1}^\top\right) X^\top Y P$$

Then from Conditions 1, 2 and 3 of Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>,
we have $d\_{B}=\hat{B\_{H\_{1}}^\*}-\hat{B\_{H\_{0}}^\*}$ has a
limiting normal distribution with

$$\begin{aligned}
        \mathbb{E} d\_{B} & =\left(\left(X^\top X\right)^{-1}-A\_{1}\left(A\_{1} X^\top X A\_{1}\right)^{-1} A\_{1}^\top\right) X^\top X B P \\
        & =\left(I\_{d}-A\_{1}\left(A\_{1}^\top X^\top X A\_{1}\right)^{-1} A\_{1}^\top X^\top X\right) B P
    \end{aligned}$$

and
$\sqrt{n} \operatorname{vec}\left(d\_{B}-\mathbb{E} d\_{B}\right) \stackrel{d}{\rightarrow} N\_{d k}(0, V)$,
where

$$\begin{aligned}
        V & =\lim \_{n \rightarrow \infty}\left(\left(X^\top X / n\right)^{-1}-A\_{1}\left(A\_{1}^\top X^\top X A\_{1} / n\right)^{-1} A\_{1}^\top\right) \otimes\left(P^\top \Sigma P\right) \\
        & =\left(M\_{X}^{-1}-A\_{1}\left(A\_{1}^\top M\_{X} A\_{1}\right)^{-1} A\_{1}^\top\right) \otimes\left(P^\top \Sigma P\right) .
    \end{aligned}$$

Hence we have
$H\_{P}=d\_{B}^\top X^\top X d\_{B} \stackrel{d}{\rightarrow} W\_{k}\left(m, I\_{k}, C\right)$,
where *C* is the non-central parameter for the Wishart distribution and
*C* = *P*<sup>⊤</sup>*B*<sup>⊤</sup>*X*<sup>⊤</sup>(*I*<sub>*n*</sub>−*X**A*<sub>1</sub>(*A*<sub>1</sub><sup>⊤</sup>*X*<sup>⊤</sup>*X**A*<sub>1</sub>)<sup>−1</sup>*A*<sub>1</sub><sup>⊤</sup>*X*<sup>⊤</sup>)*X**B**P*.

Hence we prove that *H*<sub>*P*</sub> follows the same asymptotic
distribution in both the Gaussian and non-Gaussian cases under
Conditions 1, 2, and 3. ◻

# Estimation of Optimal Projection Direction

Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>
explains how the optimal projection direction is determined by *Σ* and
*B*. We shall devise an estimation technique for the optimal direction.
We recommend setting *k* = *m* for the mentioned hypothesis above with
small *m*, as we did for the one-sample and two-sample mean testing
issues, and the optimal projection direction is
*Σ*<sup>−1</sup>*B*<sup>⊤</sup>*A*<sub>0</sub><sup>⊤</sup> according to
Theorems 1(*A*2) and (*A*3). As seen below, the ideal projection matrix
can be computed column by column.
*P*<sub>0, *k*</sub> = *Σ*<sup>−1</sup>*B*<sup>⊤</sup>*A*<sub>0, *k*</sub><sup>⊤</sup>  ∀*k* = 1, 2, ⋯, *m*
where *A*<sub>0, *k*</sub> is the *k*<sup>*t**h*</sup> row of
*A*<sub>0</sub> and *P*<sub>0, *k*</sub> is the *k*<sup>*t**h*</sup>
column of *P*<sub>0</sub>. So without loss of generality we can assume,
*m* = 1. Now we have to estimate the optimal projection direction
*P*<sub>0</sub> = *Σ*<sup>−1</sup>*B*<sup>⊤</sup>*A*<sub>0</sub><sup>⊤</sup>
So basically the estimate can be obtained as,
*P̂* = *Σ̂*<sup>−1</sup>*B̂*<sup>⊤</sup>*A*<sub>0</sub><sup>⊤</sup>
However, when *p* is larger than *m*, then *Σ̂* may become singular. To
deal with this problem, we will use an shrinkage estimator like ridge
regression to find,
*P̂*<sub>ridge</sub> = (*λ*<sub>0</sub>*I*<sub>*p*</sub>+*Σ̂*)<sup>−1</sup>*B̂*<sup>⊤</sup>*A*<sub>0</sub><sup>⊤</sup>
where *λ*<sub>0</sub> \> 0. Now to find the estimate of and *B̂*, we use
the ordinary least squares i.e.,
*B̂*<sub>LS</sub> = (*X*<sup>⊤</sup>*X*)<sup>−1</sup>*X*<sup>⊤</sup>*Y*
and estimate *Σ̂* by the sample covariance of the residuals i.e.,
$$\hat{\Sigma}\_{\text{LS}} = \dfrac{1}{n-d}(Y-X\hat{B}\_{\text{LS}})^\top(Y-X\hat{B}\_{\text{LS}})$$
Although having prior knowledge of *B* and *Σ*, using specific
estimators for them can increase the power of the test, we advise using
general estimators of *B* and *Σ* in this article to reflect the general
case. Now from *B̂*<sub>LS</sub> and *Σ̂*<sub>LS</sub> we have,
*P̂*<sub>LS</sub> = (*λ*<sub>0</sub>*I*<sub>*p*</sub>+*Σ̂*<sub>LS</sub>)<sup>−1</sup>*B̂*<sub>LS</sub><sup>⊤</sup>*A*<sub>0</sub><sup>⊤</sup>

# U-Projection test

Consider the general multivariate linear model
(<a href="#model" data-reference-type="ref" data-reference="model"><span
class="math display"><em>m</em><em>o</em><em>d</em><em>e</em><em>l</em></span></a>)
and the hypothesis *A*<sub>0</sub>*B* = 0. From Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>,
*P*<sub>0</sub> = *Σ*<sup>−1</sup>*B*<sup>⊤</sup>*A*<sub>0</sub><sup>⊤</sup>
is the asymptotic optimal projection matrix of dimension *p* × *m*. The
null hypothesis of the projection test with direction *P*<sub>0</sub> is
*H*<sub>0*P*</sub> : *A*<sub>0</sub>*B**P*<sub>0</sub> = *A*<sub>0</sub>*B**Σ*<sup>−1</sup>*B*<sup>⊤</sup>*A*<sub>0</sub><sup>⊤</sup> = 0.
Since
*A*<sub>0</sub>*B**Σ*<sup>−1</sup>*B*<sup>⊤</sup>*A*<sub>0</sub><sup>⊤</sup>
is positive semi-definite, it is equivalent to the test
*H*<sub>0*P*</sub> : *t**r*(*A*<sub>0</sub>*B**P*<sub>0</sub>) = *t**r*(*A*<sub>0</sub>*B**Σ*<sup>−1</sup>*B*<sup>⊤</sup>*A*<sub>0</sub><sup>⊤</sup>) = 0.
Following the sample-splitting test procedure, we split the data
(*X*,*Y*) into two parts (*X*<sub>1</sub>,*Y*<sub>1</sub>) and
(*X*<sub>2</sub>,*Y*<sub>2</sub>). The sample sizes of
*X*<sub>1</sub>, *Y*<sub>1</sub> are *k*, and the sample sizes of
*X*<sub>2</sub>, *Y*<sub>2</sub> are *n* − *k*. The first sample is used
an estimator *P̂*<sub>(*X*<sub>1</sub>,*Y*<sub>1</sub>)</sub> of the
*m*-dimensional projection and sample-splitting test statistic on the
projected second sample is calculated as follows:
*t**r*(*A*<sub>0</sub>*B̂*<sub>(*X*<sub>2</sub>,*Y*<sub>2</sub>)</sub>*P̂*<sub>(*X*<sub>1</sub>,*Y*<sub>1</sub>)</sub>) = *t**r*(*A*<sub>0</sub>(*X*<sub>2</sub><sup>⊤</sup>*X*<sub>2</sub>)<sup>−1</sup>*X*<sub>2</sub><sup>⊤</sup>*Y*<sub>2</sub>*P̂*<sub>(*X*<sub>1</sub>,*Y*<sub>1</sub>)</sub>)
So the subsequent U-projection statistic can be written as,
$$U_p = \dfrac{1}{\|\Gamma\|} \sum\_{\gamma\in\Gamma} tr\left( A_0(X\_{-\gamma}^\top X\_{-\gamma})^{-1}X\_{-\gamma}^\top Y\_{-\gamma}\hat{P}\_{(X\_{\gamma},Y\_{\gamma})} \right) \tag{3}$$
where *Γ* = {*γ*\|*γ*⊂{1,2,⋯,*n*},\|*γ*\| = *k*,
rank(*X*<sub>*γ*</sub>) \> *d*, rank(*X*<sub>−*γ*</sub>) ≥ *d*}, and
*X*<sub>*γ*</sub>, *Y*<sub>*γ*</sub>, *X*<sub>−*γ*</sub>, *Y*<sub>−*γ*</sub>
are subsamples of *X* and *Y* with and without index *γ*, respectively.
Now we can have the asymptotic distribution of *U*<sub>*p*</sub> using
the following theorem:

**Theorem 2**. *Consider *H*<sub>0</sub> : *A*<sub>0</sub>*B* = 0 under
model (1). Let
$$h\_{P}\left(Z\_{1,1}, \ldots, Z\_{k+d, 1} ; Z\_{1,2}, \ldots, Z\_{k+d, 2}\right)
= \frac{1}{\|\Gamma\|} \sum\_{\gamma \in \Gamma} \operatorname{diag}\left(A\_{0}\left(Z\_{-\gamma, 1}^\top Z\_{-\gamma, 1}\right)^{-1} Z\_{-\gamma, 1}^\top Z\_{-\gamma, 2} \hat{P}\_{\left(Z\_{\gamma, 1}, Z\_{\gamma, 2}\right)}\right)$$
where  − *γ* = {1, …, *k* + *d*} ∖ *γ*, *Γ* = {*γ*, *γ* ⊂ {1, 2, …, *k*+
*d*},\|*γ*\|=*k*,rank(*Z*<sub>*γ*, 1</sub>)\>*d*,rank(*Z*<sub>−*γ*, 1</sub>)=*d*}
; *Z*<sub>*γ*, *i*</sub>, *Z*<sub>−*γ*, *i*</sub> are samples of
*Z*<sub>…, *i*</sub> with and without index *γ* for *i* = 1, 2; and
*P̂*<sub>(*Z*<sub>*γ*, 1</sub>,*Z*<sub>*γ*, 2</sub>)</sub> is the
projection direction estimated with
(*Z*<sub>*γ*, 1</sub>,*Z*<sub>*γ*, 2</sub>). Then with the conditions
that *k* is fixed and that
𝔼\[tr(cov(*h*<sub>*P*</sub>(*Z*<sub>1, 1</sub>,…,*Z*<sub>*k* + *d*, 1</sub>;*Z*<sub>1, 2</sub>,…,*Z*<sub>*k* + *d*, 2</sub>)))\] \< *C*<sub>0</sub>
for some fixed *C*<sub>0</sub> \> 0 and
(*Z*<sub>*i*, 1</sub>,*Z*<sub>*i*, 2</sub>), *i* = 1, …, *k* + *d*
i.i.d. from the distribution of (*X*,*Y*), we have asymptotic normality
of the general U-projection statistic (3):
$$\sqrt{n}\left(U\_{P}-\mathbb{E} \[U\_{P}\]\right) \stackrel{d}{\rightarrow} N\left(0,(k+d)^{2} 1\_{m}^\top \Xi\_{1} 1\_{m}\right)$$
where 1<sub>*m*</sub> is a vector full of one with length *m*, and
*Ξ*<sub>1</sub> = cov (*h*<sub>*P*</sub>(*Z*<sub>1, 1</sub>,*Z*<sub>2, 1</sub>,…,*Z*<sub>*k* + *d*, 1</sub>;*Z*<sub>1, 2</sub>,…,*Z*<sub>*k* + *d*, 2</sub>),*h*<sub>*P*</sub>(*Z*<sub>1, 1</sub>,*Z*<sub>2, 1</sub><sup>′</sup>,…,*Z*<sub>*k* + *d*, 1</sub><sup>′</sup>;*Z*<sub>1, 2</sub>,*Z*<sub>2, 2</sub><sup>′</sup>,…,*Z*<sub>*k* + *d*, 2</sub><sup>′</sup>))
for (*Z*<sub>*i*, 1</sub>,*Z*<sub>*i*, 2</sub>), *i* = 1, …, *k* + *d*
and
(*Z*<sub>*i*, 1</sub><sup>′</sup>,*Z*<sub>*i*, 2</sub><sup>′</sup>), *i* = 2, …, *k* + *d*
i.i.d. from the distribution of (*X*,*Y*). Furthermore, under the null
hypothesis, we have 𝔼\[*U*<sub>*P*</sub>\] = 0. *

*Proof.* Let
$$M\_{P}=\frac{1}{\|\Gamma\|} \sum\_{\gamma \in \Gamma} \operatorname{diag}\left(A\_{0}\left(X\_{-\gamma}^{T} X\_{-\gamma}\right)^{-1} X\_{-\gamma}^{T} Y\_{-\gamma} \hat{P}\_{\left(X\_{\gamma}, Y\_{\gamma}\right)}\right)$$
where
*Γ* = {*γ*\|*γ*⊂{1,2,⋯,*n*},\|*γ*∣=*k*,rank(*X*<sub>*γ*</sub>)\>*d*,rank(*X*<sub>−*γ*</sub>)≥*d*},
and *X*<sub>*γ*</sub>, *Y*<sub>*γ*</sub>,
*X*<sub>−*γ*</sub>*Y*<sub>−*γ*</sub> are sub-samples of *X* and *Y* with
and without index *γ* respectively. Then we have
*U*<sub>*P*</sub> = 1<sub>*m*</sub><sup>*T*</sup>*M*<sub>*P*</sub>.
*M*<sub>*P*</sub> can be rewritten as a U-statistic from the kernel
*h*<sub>*P*</sub>. From the asymptotic normality result of U-statistics,
we have
$$\sqrt{n}\left(M\_{P}-\mathbb{E} \[M\_{P}\]\right) \stackrel{d}{\rightarrow} N\left(0,(k+d)^{2} \Xi\_{1}\right)$$
where *Ξ*<sub>1</sub> is as defined in the theorem. Hence the asymptotic
result of *U*<sub>*P*</sub> follows from the fact that
*U*<sub>*P*</sub> = 1<sub>*m*</sub><sup>*T*</sup>*M*<sub>*P*</sub>.
Moreover, it is easy to get
𝔼<sub>*H*<sub>0</sub></sub>\[*U*<sub>*P*</sub>\] = 0 under the null
hypothesis. This completes the proof of the theorem.◻

Although the asymptotic distribution of the U-projection statistic can
be obtained under certain assumptions, as shown in Theorem
<a href="#thm2" data-reference-type="ref" data-reference="thm2">2</a>,
it may not be very useful here for a number of reasons:

-   If the sample size *n* is small, the asymptotic distribution will be
    far from the real distribution.

-   In general, the relationship between the asymptotic variance under
    the null and the covariance is rather complex.

-   The (asymptotic) distribution of the statistic is likewise affected
    by the projection direction estimation method. The distribution will
    change if we use alternative estimation methods for the projection
    direction.

So using randomization, we can calculate the p-value for the test using
the following algorithm:

**P-value Calculation using Randomization**

**Input** *n* × *d* dimensional matrix *X*, *n* × *p* dimensional matrix
*Y*, *m* × *d* dimensional matrix *A*<sub>0</sub>,
*d* \< *k* ≤ *n* − *d*, and randomization times *N*.

Calculate the *U*-projection statistic *U*<sub>0</sub> on the original
dataset *X* and *Y*.<br> Let *i* = 0<br> **For** (j=1:N){<br>   Do
randomization, get *X*′ and *Y*′<br>   Calculate the *U*-projection
statistic *U*′ on the randomized dataset *X*′ and *Y*′.<br>   **If**
(*U*′≥*U*){ <br>     *i* = *i* + 1<br>   }<br> }<br> Calculate p-value
by $p=\dfrac{i}{N}$

## Examples

### One Sample Mean Testing:

For one-sample mean testing, suppose we have *n* samples *Y* from a
*p*-dimensional multivariate distribution with mean *μ* and covariance
*Σ*, and we want to test whether *μ* = 0. Now for this setting we can
rewrite our model as,
*Y* = 1<sub>*n*</sub>*B* + *E*
where 1<sub>*n*</sub> is a column vector full of 1 of length *n*,
*B* = *μ*<sup>⊤</sup>, *E* is the *n* × *p* random error matrix, and the
null hypothesis becomes
*H*<sub>0</sub> : *B* = 0
i.e., *A*<sub>0</sub> = 1 and is of rank 1 in this problem. By Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>,
the optimal projection dimension is 1 and the optimal projection
direction is *Σ*<sup>−1</sup>*μ*. Now we can find the estimates of *Σ*
and *μ* through least squares and subsequently estimate the optimal
projection direction. So the estimated projection matrix becomes,
*P̂* = (*λ*<sub>0</sub>*I*<sub>*p*</sub>+*Σ̂*)<sup>−1</sup>*B̂*<sub>LS</sub><sup>⊤</sup>1<sub>*n*</sub><sup>⊤</sup> = (*λ*<sub>0</sub>*I*<sub>*p*</sub>+*Σ̂*)<sup>−1</sup>*Ȳ*
Thus, the U-projection statistic becomes
$$U\_{\mathrm{P}}=  \dfrac{1}{\binom{n}{k}}\sum\_{\gamma \in \Gamma}\text{tr}\left(\bar{Y}\_{-\gamma}^\top \times\left\\\lambda\_{0} I\_{p}+S\_{Y\_{\gamma}}\right\\^{-1} \times \bar{Y}\_{\gamma}\right)$$
where for *k* is the number of independent samples from *Y* to estimate
the projection direction, *Γ* is collections of all subsets of {1,…,*n*}
with size *k*, and *Y*<sub>*γ*</sub> and *Y*<sub>−*γ*</sub> are subsets
of *Y* with and without index *γ* correspondingly.

### Testing of Significance of Predictors

Consider the general linear model (1). Testing the significance of a
particular predictor, *X*<sup>(*i*)</sup>, is of interest. In other
words, take
*H*<sub>0</sub> : *A*<sub>0</sub>*B* = 0
where *A*<sub>0</sub> = *e*<sub>*i*</sub><sup>⊤</sup> and
*e*<sub>*i*</sub> is a column vector in *d* dimensions, with the
*i*<sup>*t**h*</sup> element being one and the remaining elements being
zero. The above hypothesis testing problem is also testing for the
conditional independence of *Y* and *X*<sup>(*i*)</sup> given other
(*d*−1) predictors under the normality assumption. This is because, in
the normal case, testing whether the regression coefficients equal zero
or not is equivalent to testing of conditional independence. According
to Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>,
the ideal projection dimension is 1 and the optimal projection direction
is *Σ*<sup>−1</sup>*B*<sup>⊤</sup>*e*<sub>*i*</sub>. Now using the
previously demonstrated estimation method, we can estimate the
projection matrix by,
*P̂* = (*λ*<sub>0</sub>*I*<sub>*p*</sub>+*Σ̂*)<sup>−1</sup>*B̂*<sub>LS</sub><sup>⊤</sup>*e*<sub>*i*</sub> = (*λ*<sub>0</sub>*I*<sub>*p*</sub>+*Σ̂*)<sup>−1</sup>*B̂*<sub>*i*</sub>
Thus, the U-projection statistic becomes
$$U\_{\mathrm{P}}=  \dfrac{1}{\binom{n}{k}}\sum\_{\gamma \in \Gamma}\text{tr}\left(\hat{B}\_i^\top(Y\_{-\gamma}) \times\left\\\lambda\_{0} I\_{p}+S\_{Y\_{\gamma}}\right\\^{-1} \times \hat{B}\_i(Y\_{\gamma}) \right)$$
where for *k* is the number of independent samples from *Y* to estimate
the projection direction, *Γ* is collections of all subsets of {1,…,*n*}
with size *k*, and *Y*<sub>*γ*</sub> and *Y*<sub>−*γ*</sub> are subsets
of *Y* with and without index *γ* correspondingly.

# Future Direction: Multiple Sample Mean Testing

For general unbalanced multiple-sample mean testing problem, suppose for
*k* = 1, 2, ⋯, *K*, we have *n*<sub>*k*</sub> samples *Y*<sub>*k*</sub>
from a *p*-dimensional multivariate distribution *F*<sub>*k*</sub> with
mean *μ*<sub>*k*</sub> and covariance *Σ*, and we want to test whether
*μ*<sub>1</sub> = *μ*<sub>2</sub> = ⋯ = *μ*<sub>*k*</sub>. The problem
can be reformulated as
$$\begin{pmatrix}
    Y_1\\
    \cdots\\
    Y_k
\end{pmatrix} = \begin{pmatrix}
    1\_{n_1} & 0 & \cdots\\
    \cdots & \cdots & \cdots\\
    0 & \cdots & 1\_{n_k}
\end{pmatrix} B + E$$
where
*B* = (*μ*<sub>1</sub>,*μ*<sub>2</sub>,⋯,*μ*<sub>*k*</sub>)<sup>⊤</sup>,
and *E* is the $\left( \sum\_{i=1}^k n_i\right)\times p$ random error
matrix. We want to test
*H*<sub>0</sub> : *A*<sub>0</sub>*B* = 0
where
*A*<sub>0</sub> = {*a*<sub>*i**j*</sub>}<sub>1 ≤ *i* ≤ *K* − 1, 1 ≤ *j* ≤ *K*</sub>
with *a*<sub>*i**i*</sub> = *K* − 1 and
*a*<sub>*i**j*</sub> =  − 1, ∀*i* ≠ *j*. In this problem,
*A*<sub>0</sub> is of rank *K* − 1. By Theorem
<a href="#thm1" data-reference-type="ref" data-reference="thm1">1</a>,
the optimal projection dimension is less than or equal to *K* − 1 and
the optimal *p* × (*K*−1) projection direction matrix is
*K**Σ*<sup>−1</sup>(*μ*<sub>1</sub>−*μ̄*,⋯,*μ*<sub>*k* − 1</sub>−*μ̄*),
where $\bar{\mu} = \frac{1}{k}\sum\_{k=1}^K \mu_k$. In particular,when
*K* = 2, the optimal projection dimension is 1 and the optimal
projection direction is
*Σ*<sup>−1</sup>(*μ*<sub>1</sub>−*μ*<sub>2</sub>). Now using the
previously demonstrated estimation method, we can estimate the
projection matrix by,
*P̂*<sub>LS</sub> = (*λ*<sub>0</sub>*I*<sub>*p*</sub>+*Σ̂*)<sup>−1</sup>*B̂*<sub>LS</sub><sup>⊤</sup>*A*<sub>0</sub><sup>⊤</sup> = (*λ*<sub>0</sub>*I*<sub>*p*</sub>+*Σ̂*)<sup>−1</sup>(*Ȳ*<sub>1</sub>−*Ȳ*<sub>2</sub>)
Thus, the U-projection statistic becomes
$$U\_{\mathrm{LS}}=  \dfrac{1}{\binom{n_1}{k_1}\binom{n_2}{k_2}} \sum\_{\gamma\_{1} \in \Gamma\_{1}} \sum\_{\gamma\_{2} \in \Gamma\_{2}}\left(\bar{Y}\_{-\gamma\_{1}, 1}^\top-\bar{Y}\_{-\gamma\_{2}, 2}^\top\right) \times\left\\\lambda\_{0} I\_{p}+\frac{\left(k\_{1}-1\right) S\_{Y\_{\gamma\_{1}, 1}}+\left(k\_{2}-1\right) S\_{Y\_{\gamma\_{2}, 2}}}{k\_{1}+k\_{2}-2}\right\\^{-1} \times\left(\bar{Y}\_{\gamma\_{1}, 1}-\bar{Y}\_{\gamma\_{2}, 2}\right)$$
where for *i* = 1, 2, *k*<sub>*i*</sub> is the number of independent
samples from *Y*<sub>*i*</sub> to estimate the projection direction,
*Γ*<sub>*i*</sub> is collections of all subsets of
{1,…,*n*<sub>*i*</sub>} with size *k*<sub>*i*</sub>, and
*Y*<sub>*γ*, *i*</sub> and *Y*<sub>−*γ*, *i*</sub> are subsets of
*Y*<sub>*i*</sub> with and without index *γ* correspondingly. Now
consider the following condition and theorems:

**Condition 4:**
*E*<sub>*i*</sub> = *Γ**Z*<sub>*i*</sub>   for *i* = 1, …, *n*, where
*Γ* is a *p* × *t* matrix with some *t* ≥ *p* such that
*Γ**Γ*<sup>*T*</sup>= *Σ*, and
*Z*<sub>*i*</sub> = (*Z*<sub>*i*, 1</sub>,…,*Z*<sub>*i*, *t*</sub>) are
*t*-variate independent and identically distributed random vectors
satisfying 𝔼(*Z*<sub>*i*</sub>) = 0,
var (*Z*<sub>*i*</sub>) = *I*<sub>*t*</sub>, 𝔼(*Z*<sub>*i*, *k*</sub><sup>3</sup>) = 0, 𝔼(*Z*<sub>*i*, *k*</sub><sup>6</sup>)
is uniformly bounded, and
𝔼(*Z*<sub>*i*, *l*<sub>1</sub></sub><sup>*α*<sub>1</sub></sup>*Z*<sub>*i*, *l*<sub>2</sub></sub><sup>*α*<sub>2</sub></sup>⋯*Z*<sub>*i*, *l*<sub>*s*</sub></sub><sup>*α*<sub>*s*</sub></sup>) = 𝔼(*Z*<sub>*i*, *l*<sub>1</sub></sub><sup>*α*<sub>1</sub></sup>)𝔼(*Z*<sub>*i*, *l*<sub>2</sub></sub><sup>*α*<sub>2</sub></sup>)⋯𝔼(*Z*<sub>*i*, *l*<sub>*s*</sub></sub><sup>*α*<sub>*s*</sub></sup>)
for a positive integer *s* such that
$\sum\_{l=1}^{s} \alpha\_{l} \leq 8$ and
*l*<sub>1</sub> ≠ *l*<sub>2</sub>≠ ⋯ ≠ *l*<sub>*s*</sub>.

**Condition 5:**
∥*A*<sub>0</sub>*B*∥<sub>F</sub><sup>2</sup> = *o*(*p*/*n*).

**Condition C1:** There exists a uniformly bounded positive integer
*q* \< *p* such that
$\frac{\sqrt{n} \lambda\_{q}}{\operatorname{tr}(\Sigma)} \rightarrow \infty$
and *λ*<sub>*q* + 1</sub> is uniformly bounded from above. The smallest
eigenvalue *λ*<sub>*p*</sub> is uniformly bounded from below.

**Condition C2:**
tr (*Σ*<sup>4</sup>) = *o*(tr<sup>2</sup>(*Σ*<sup>2</sup>)).

**Theorem 3**. *Suppose the covariance *Σ* satisfies Condition C1. Under
high-dimensional setting
*n*<sub>1</sub> + *n*<sub>2</sub> = *o*(tr(*Σ*)), Conditions 3, 4, 5,
and conditions that
*k*<sub>*i*</sub>/*n*<sub>*i*</sub> → *γ*<sub>*i*</sub> ∈ (0,1) for
*i* = 1, 2,
*n*<sub>1</sub>/(*n*<sub>1</sub>+*n*<sub>2</sub>) → *κ* ∈ (0,1), then
*U*<sub>LS</sub> has an asymptotically normal distribution with a
uniformly bounded positive *λ*<sub>0</sub>. More specifically, we have*

*
$$\frac{\lambda\_{0}\left(U\_{\mathrm{LS}}-\mathbb{E}\left(U\_{\mathrm{LS}}\right)\right)}{\sigma\_{n}} \stackrel{d}{\rightarrow} N(0,1), \quad \text { and } \quad \lambda\_{0} \mathbb{E}\left(U\_{\mathrm{LS}}\right)-\left\\W\_{q+1} \mu\_{n}\right\\\_{2}^{2}=o\left(\left\\W\_{q+1} \mu\_{n}\right\\\_{2}^{2}\right),$$
where *q* is the number of divergent eigenvalues of *Σ* as specified in
Condition C1, *W*<sub>*q* + 1</sub> is the projection matrix onto the
linear span of eigenspaces of *Σ* corresponding to the smallest
*p* − *q* eigenvalues,
*λ*<sub>*q* + 1</sub>, *λ*<sub>*q* + 2</sub>, …, *λ*<sub>*p*</sub>, *μ*<sub>*n*</sub>
is the mean difference of the two populations, and
$\sigma\_{n}^{2}=\left(\frac{2}{n\_{1}^{2}}+\frac{2}{n\_{2}^{2}}+\frac{4}{n\_{1} n\_{2}}\right) \sum\_{i=q+1}^{p} \lambda\_{i}^{2}$.
*

**Theorem 4**. *Suppose the covariance *Σ* satisfies Condition C2. Under
high-dimensional setting
*n*<sub>1</sub> + *n*<sub>2</sub> = *O*(tr(*Σ*)), Conditions 3-5, and
conditions that
*k*<sub>*i*</sub>/*n*<sub>*i*</sub> → *γ*<sub>*i*</sub> ∈ \[0, 1) for
*i* = 1, 2,
*n*<sub>1</sub>/(*n*<sub>1</sub>+*n*<sub>2</sub>) → *κ* ∈ (0,1), and*

*
(*k*<sub>1</sub>+*k*<sub>2</sub>) = *o*(tr(*Σ*)/*λ*<sub>max</sub>(*Σ*)<sup>2</sup>)
where *λ*<sub>max</sub>(*Σ*) is the largest eigenvalue of
*Σ*, *λ*<sub>0</sub>*U*<sub>LS</sub> is asymptotically normally
distributed and has the same asymptotic variance and similar expectation
with *T*<sub>CQ</sub>,*

*
$$\begin{aligned}
& \frac{\lambda\_{0}\left(U\_{\mathrm{LS}}-\mathbb{E} U\_{\mathrm{LS}}\right)}{\sqrt{\operatorname{var} T\_{\mathrm{CQ}}}} \stackrel{d}{\rightarrow} N(0,1), \quad \text { and } \\
& \left\|\mathbb{E}\left(\lambda\_{0} U\_{\mathrm{LS}}-T\_{\mathrm{CQ}}\right)\right\|=o\left(\mathbb{E} T\_{\mathrm{CQ}}\right),
\end{aligned}$$
and *U*<sub>LS</sub> and *T*<sub>CQ</sub> have the same asymptotic
power, *β*<sub>*U*<sub>LS</sub></sub>(*μ*<sub>1</sub>−
*μ*<sub>2</sub>) − *β*<sub>*T*<sub>CQ</sub></sub>(*μ*<sub>1</sub>−*μ*<sub>2</sub>) → 0.
*

Above *T*<sub>*C**Q*</sub> denote the two-sample statistic in Chen and
Qin (2010) (Chen and Qin 2010). Now it can be observed that for high
correlation covariance *Σ*, Theorem
<a href="#thm3" data-reference-type="ref" data-reference="thm3">3</a>
gives the asymptotic distribution of *U*<sub>*L**S*</sub> and for low
correlation covariance *Σ*, Theorem
<a href="#thm4" data-reference-type="ref" data-reference="thm4">4</a>
does the same. Note that, Theorem
<a href="#thm3" data-reference-type="ref" data-reference="thm3">3</a>
and
<a href="#thm4" data-reference-type="ref" data-reference="thm4">4</a>
provide different formulas to calculate the asymptotic variance under
different structures of covariances. So we need to calculate the
p-values using randomization only. Now when the number of samples
becomes greater i.e., *K* \> 2, although the p-values can be obtained
through randomization, there are no proper documentation of the
U-projection statistic. We think this is an area which is still under
research and can be developed further.

# References

Anderson, TW. 2003. “An Introduction to Statistical Multivariate
Analysis.” John-Wiley.

Bai, Zhidong, and Hewa Saranadasa. 1996. “Effect of High Dimension: By
an Example of a Two Sample Problem.” *Statistica Sinica*, 311–29.

Chen, Song Xi, and Ying-Li Qin. 2010. “A Two-Sample Test for
High-Dimensional Data with Applications to Gene-Set Testing.”

Muirhead, Robb J. 2009. *Aspects of Multivariate Statistical Theory*.
John Wiley & Sons.

Runze Li, Changcheng Li. 2022. “Linear Hypothesis Testing in Linear
Models with High-Dimensional Responses.” *Journal of the American
Statistical Association* 117 (540): 1738–50.
