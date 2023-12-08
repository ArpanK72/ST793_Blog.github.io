---
author: <center> Shuvrarghya Ghosh &emsp; &emsp; Arpan Kumar </center>
date: <center> December 07 2023 </center> 
bibliography:
- biblio.bib
title: Linear Hypothesis Testing in Linear Models With High-Dimensional
  Responses - A Review
header-includes: \usepackage[ruled,vlined]{algorithm2e}
output:
  html_document:
    toc: true
    toc_float:
      collapsed: true
      smooth_scroll: true
---
<hr>
<style>
body{
text-align: justify
}
</style>
# Introduction

This project studies the paper titled ``Linear Hypothesis Testing in Linear Models With High-Dimensional Responses" by Li and Li (2022) [@runze2022linear] where the authors talk about the high-dimensional $(p\gg n)$ scenarios where Hotelling $T^2$ cannot be used to do one-sample and two-sample testing because of the singularity of sample covariance matrix. Although, there have several papers on solving this problem and its subsequent problems, the paper of our interest develops general linear hypothesis testing frameworks using projection matrices, which can be applied to one-sample, two-sample and multiple sample means testing in high-dimensions. The main goal of this project is to understand, simplify the concepts and lucidly explain this new projection test and its theoretical properties through a blog. As we proceed further in this project, we both shall elucidate the theorems in the main paper as far as possible.

One sample and two sample tests have gained a lot of focus in the last
two decades. The Hotelling $T^2$ test for the one-sample and two-sample
mean problem is the uniformly most powerful test invariant under affine
transformations for fixed or low-dimensional data. The singularity of
the sample covariance matrix renders the Hotelling $T^2$ test invalid
when the dimension $p$ of observations exceeds the sample size $n$.
Numerous tests for high-dimensional one and two-sample mean testing have
been proposed to address the singularity. When examining the effects of
dimensionality on two-sample mean testing, Bai and Saranadasa (1996)
[@bai1996effect] promoted sum-of-squares-type statistics for the mean
difference and disregarded correlation matrices. Their test circumvents
the issue brought on by the sample covariance matrices' rank deficiency.
There are also various other tests that have been created to address any
potential issues that may arise in this area. To create the methods,
some of them used supremum-type statistics, U-statistics, and so on. The
authors of this paper offer a projection test for a specific setting.
The correlation information is used by the projection test to increase
the test's power. They build the projection test by first projecting the
data into low dimension and then performing the test on the projected
data. They demonstrated that there are low-dimensional projection
matrices that leverage covariance information to maximize the power of
projection tests. They further demonstrate that the dimension of the
optimal projections is not more than $m$, and they derive the $m$
optimal projection directions to maximize power. They adopt a
sample-splitting approach to perform projection tests, in which they use
the first portion of the sample for estimate of projection matrices and
the second part for testing. They also suggest U-projection tests to
increase the test's power based on the sample-splitting technique by
creating a U-type test statistic from all samples. They have also
established the theoretical properties of the U-projection test.

# Notations

Let us consider the classical linear model setup:
$$Y=XB+E \label{model}\tag{1}$$ 
where $Y$ is an $n\times p$ response matrix,
$X$ is an $n\times d$ design matrix whose dimension $d$ is fixed, $B$ is
a $d\times p$ regression coefficient matrix, and $E$ is an $n\times p$
error matrix with mean zero and covariance matrix $I_n \otimes \Sigma$
for a positive definite matrix $\Sigma$. We denote
$\text{vec}(A)=(a_1^\top,\cdots,a_n^\top)^\top$ for
$A = (a_1,a_2,\cdots,a_n)^\top$ from here onward. We are interested in
testing the following hypotheses: $$\begin{aligned}
    H_0&: A_0B = 0 \\
    H_1 & : A_0B\ne 0 \notag
\end{aligned}$$ where $A_0$ is an $m\times d$ known constant matrix with
both $m$ and $d$ fixed and $A_0$ is of full row rank which implies
$m\leq d$. Under usual assumptions of dimensions, i.e., $n < p$, the
likelihood ratio test (LRT) and its asymptotic convergence to
chi-squared distribution holds good. However, in high dimensions when
$p >> n$, the usual properties of LRT fails. In this blog, we elucidate
a methods of testing hypotheses in such scenarios using projection
matrices that project the high-dimension problem to low dimensions such
that a tractable solution is achieved using the conventional LRT.

# Optimal Projection Direction

Firstly, let us redefine the problem. Let $P$ be a $p \times r$ full
column rank projection matrix with $r \ll p$ and $r<n$, consider the
$P$-projected model:

$$Y^{\star}=X B^{\star}+E^{\star}$$

where $Y^{\star}=Y P$, $B^{\star}=B P$, and $E^{\star}=E P$, which is
the $n \times r$ projected error matrix with mean zero and
$\operatorname{cov}\left(\operatorname{vec}\left(E^{\star}\right)\right)=I_{n} \otimes$
$\left(P^\top \Sigma P\right)$

The corresponding projected hypothesis becomes

$$H_{0 P}: A_{0} B^{\star}=0 \quad \text { versus } \quad H_{1 P}: A_{0} B^{\star} \neq 0 \tag{2}$$

and the projected test statistic is

$$\Lambda=\frac{\left|G_{P}\right|}{\left|G_{P}+H_{P}\right|}$$

where $G_{P}$ is the residual sum of squares under $H_{1 P}$,
$G_{P}+H_{P}$ is the residual sum of squares under $H_{0 P}$.
$G_{P}=P^\top Y^\top\left(I_{n}-P_{H_{1}}\right) Y P$,
$P_{H_{1}}=X\left(X^\top X\right)^{-1} X^\top$ and
$G_{P}+H_{P}= P^\top Y^\top\left(I_{n}-P_{H_{0}}\right) Y P$,
$P_{H_{0}}=X A_{1}\left(A_{1}^\top X^\top X A_{1}\right)^{-1} A_{1}^\top X^\top$,
and $A_{1}$ is a $d \times(d-m)$ matrix defined by
$A_{1}=\left(A_{0}^\top\right)^{\perp}$, which means that
$\left(A_{0}^\top, A_{1}\right)$ forms an orthogonal matrix. The test is
the LRT under normality assumption and is equivalent to many useful test
in various cases, like Hotelling $T^{2}$ in the one-sample mean testing.

The one-sample and two-sample mean tests can be written as special
instances of the mentioned hypotheses test if $X$ and $A_0$ are
appropriately specified. As a result, the approach suggested in this
report applies directly to the one-sample, two-sample mean and even
multiple-sample mean problems. The projected null hypothesis $H_{0 P}$
is rejected for small $\Lambda$, and $H_{0}$ is rejected if $H_{0 P}$ is
rejected. In general, $H_{0}$ and $H_{0 P}$ are not equivalent. In this
article, we will show that there exists an optimal projection direction
$P$ which makes $H_{0 P}$ equivalent to $H_{0}$ and also maximizes the
power of the test.

The problem how to construct an optimal direction can be divided into
two sub-problems:

-   to find the dimension $r$ of the optimal projection direction $P$,
    and

-   to find the optimal $P$ of a particular dimension.

These issues are addressed in Theorem [1](#thm1){reference-type="ref"
reference="thm1"}. We assume the following conditions on the
multivariate linear model (1) to derive the asymptotically optimal
projection direction matrix.

**Condition 1**:
$\exists C_{1} > 0 \ni \left\|X X^\top\right\|_{\infty}<C_{1}$, where
$\|A\|_{\infty}=\max _{i, j}\left|a_{i, j}\right|$ for
$A=\left(a_{i, j}\right)_{i, j}$.

**Condition 2**: The limit
$\lim _{n \rightarrow \infty} n^{-1}\left(X^\top X\right)=M_{X}$ exists,
where $M_{X}$ is a nonsingular $d \times d$ matrix.

**Condition 3**:
$\exists C_{2} > 0 \ni \mathbb{E}\left(E_{i j}^{4}\right)<C_{2}$,
$i=1, \ldots, n, j=1, \ldots, p$ for elements $E_{i, j}$ in the error
matrix.

Conditions 1-3 are quite mild and are used to guarantee asymptotic
normality of least squares estimate in linear models. For example,
Conditions 1 and 2 hold in the one-sample mean testing problem
automatically, and they also hold in the two-sample mean testing problem
if $\frac{n_{1}}{n_{1}+n_{2}} \rightarrow \kappa \in(0,1)$, where
$n_{1}$ and $n_{2}$ are sample sizes of the first and the second
samples, respectively.

::: {#thm1 .theorem}
**Theorem 1** (Optimal Projection Direction). *Suppose that the error
matrix $E$ follows a matrix normal distribution. Under Conditions 1-3,
the following statements are valid.*

*(A1) Let
$W=\Sigma^{-1 / 2} B^\top X^\top\left(I_{n}-P_{H_{0}}\right) X B \Sigma^{-1 / 2}$,
where
$P_{H_{0}}=X A_{1}\left(A_{1}^\top X^\top X A_{1}\right)^{-1} A_{1}^\top X^\top$,
and $A_{1}$ is a $d \times$ $(d-m)$ matrix which is defined by
$A_{1}=\left(A_{0}^\top\right)^{\perp}$. $W$ is a $p \times p$ matrix of
rank $m$. Suppose $W$ admits the eigenvalue decomposition
$W=\sum_{i=1}^{m} \lambda_{0 i} R_{i} R_{i}^\top$, where
$0<\lambda_{0 m} \leq \cdots \leq \lambda_{01}$ are $m$ nonzero
eigenvalues, and $R_{i}, i=1, \ldots, m$, are the corresponding
orthogonal eigenvectors. Then for any given $r \leq m$,
$P_{0}=\Sigma^{-1 / 2}\left(R_{1}, \ldots, R_{r}\right)$ is the
projection direction matrix which is optimal among all $p \times r$
matrices for the hypothesis testing problem (2).*

*(A2) The optimal dimension of projection matrix is less than or equal
to $m$ for the hypothesis testing problem (2).*

*(A3) Set $r=m, \Sigma^{-1} B^\top A_{0}^\top$ is the projection
direction matrix which is optimal among all $p \times m$ matrices for
the hypothesis testing problem (2).*

*(B) Without normality assumption, the statements in (A1)-(A3) are still
valid in asymptotic sense by changing optimal projection to
asymptotically optimal projection which maximizes the asymptotic local
power. *
:::

::: proof
*Proof.* To prove the theorem, we use the unexplained variability or
variability due to random error components. For the ordinary linear
model, orthogonal split of sum of squares is given by
$Y^2 = Y^TP_HY + Y^T(I_n-P_H)Y$ where the later term in the summation is
accountable for the sum of squares due to errors.

We first prove Theorem [1](#thm1){reference-type="ref"
reference="thm1"}(A). For linear model and under $H_{1}$, the estimated
$Y$ is the projection of $Y$ onto the column space of $X$, denoted by
$P_{H_{1}} Y$; and under $H_{0}$, the statement that $A_{0} B=0$ is
equivalent to that $B$ is of the form $A_{1}^\top D$, where
$\left(A_{0}, A_{1}\right)$ forms an orthogonal basis and $D$ is a
$(d-m) \times p$ matrix. Hence the estimated $Y$ under $H_{0}$ is the
projection of $Y$ onto the column space of $X A_{1}^\top$, denoted by
$P_{H_{0}} Y$. So
$G=\left(Y-P_{H_{1}} Y\right)^\top\left(Y-P_{H_{1}} Y\right)=$
$Y^\top\left(I_{n}-P_{H_{1}}\right) Y, G+H=\left(Y-P_{H_{0}} Y\right)^\top\left(Y-P_{H_{0}} Y\right)=Y^\top\left(I_{n}-P_{H_{0}}\right) Y, H=Y^\top\left(P_{H_{1}}-P_{H_{0}}\right) Y$.
Under normality assumption, LRT statistic for $H_{0}$ defined in (1) is
$\frac{|G|}{|G+H|}$.

Moreover, in the projection test, $Y$ is replaced by $Y P$. The LRT
statistic in the projection test can be written as
$\frac{\left|G_{P}\right|}{\left|G_{P}+H_{P}\right|}$, where $G_{P}$ is
the residual sum of squares under $H_{1 P}$, and $G_{P}+H_{P}$ is the
residual sum of squares under $H_{0 P}$. Without loss of generality, we
assume that $\Sigma^{1 / 2} P$ is a projection matrix $Q$, i.e., $Q$ is
symmetric and idempotent, so $P=\Sigma^{-1 / 2} Q$ and

$$G_{P}=P^\top Y^\top\left(I_{n}-P_{H_{1}}\right) Y P=Q^\top \Sigma^{-1 / 2} Y^\top\left(I_{n}-P_{H_{1}}\right) Y \Sigma^{-1 / 2} Q$$

Since
$\left(I_{n}-P_{H_{1}}\right) Y \sim N\left(0,\left(I_{n}-P_{H_{1}}\right) \otimes \Sigma\right)$,
we have
$\Sigma^{-1 / 2} Y^\top\left(I_{n}-P_{H_{1}}\right) Y \sim W_{p}\left(n-d, I_{p}\right)$.
Furthermore, since $Q$ is a $p \times k$ projection matrix,
$G_{P}=Q^\top \Sigma^{-1 / 2} Y^\top\left(I_{n}-P_{H_{1}}\right) Y \Sigma^{-1 / 2} Q \sim$
$W_{k}\left(n-d, I_{k}\right)$. We have
$H_{P}=P^\top Y^\top\left(P_{H_{1}}-P_{H_{0}}\right) Y P=Q^\top \Sigma^{-1 / 2} Y^\top\left(P_{H_{1}}-P_{H_{0}}\right) Y \Sigma^{-1 / 2} Q$.
Since
$\left(P_{H_{1}}-P_{H_{0}}\right) Y \sim N\left(\left(I_{n}-P_{H_{0}}\right) X B,\left(P_{H_{1}}-P_{H_{0}}\right) \otimes \Sigma\right)$,
we have

$$\Sigma^{1 / 2} Y^\top\left(P_{I_{1}}-P_{I_{0}}\right) Y \Sigma^{-1 / 2} \sim W_{p}\left(m, I_{p}, \Sigma^{-1 / 2} B^\top X^\top\left(I_{n}-P_{I_{0}}\right) X B \Sigma^{-1 / 2}\right)$$

Hence we have

$$\begin{aligned}
H_{P} & =Q^\top \Sigma^{-1 / 2} Y^\top\left(P_{H_{1}}-P_{H_{0}}\right) Y \Sigma^{-1 / 2} Q \\
& \sim W_{k}\left(m, I_{k}, Q^\top \Sigma^{-1 / 2} B^\top X^\top\left(I_{n}-P_{H_{0}}\right) X B \Sigma^{-1 / 2} Q\right) .
\end{aligned}$$

By Cochran's Theorem, we know that $G_{P}$ is independent to $H_{P}$
given $X$ and $P$. Proof of Theorem [1](#thm1){reference-type="ref"
reference="thm1"}(A1). Here we prove a more general Theorem
[1](#thm1){reference-type="ref" reference="thm1"}(A1) for any $r \leq p$
by defining $R_{i}, i=m+1, \cdots, p$ as orthogonal eigenvectors
corresponding to eigenvalue 0 . We first show that
$P_{0}=\Sigma^{-1 / 2} Q_{0}$ is the optimal $p \times r$ projection
direction matrix under normality assumption, where
$Q_{0}=\left(R_{1}, \cdots, R_{r}\right)$ and $R_{i}$ is defined as in
Theorem [1](#thm1){reference-type="ref" reference="thm1"} .

Given $r \leq m$ fixed, it is easy to see that
$\lambda_{P}=\frac{\left|G_{P}\right|}{\left|G_{P}+H_{P}\right|}$ will
be of the same distribution under $H_{0}$ no matter what $Q$ is, so the
critical value will be the same for any $p \times r$ projection matrix
$Q$. So to compare the power of $p \times r$ projection tests, we only
need to compare
$P\left(\frac{\left|G_{P}\right|}{\left|G_{P}+H_{P}\right|}<c\right)$
for any positive $c$.

For any $p \times r$ projection matrix $Q$, it is easy to know that
$\lambda_{i} \leq \lambda_{0 i}$, where
$\lambda_{r} \leq \lambda_{r-1} \leq$ $\cdots \leq \lambda_{1}$ are the
nonzero eigenvalues of
$Q^\top \Sigma^{-1 / 2} B^\top X^\top\left(I_{n}-P_{H_{0}}\right) X B \Sigma^{-1 / 2} Q$,
and $\lambda_{0 r} \leq$
$\lambda_{0(r-1)} \leq \cdots \leq \lambda_{01}$ are the $r$ largest
eigenvalues of
$\Sigma^{-1 / 2} B^\top X^\top\left(I_{n}-P_{H_{0}}\right) X B \Sigma^{-1 / 2}$.
Let
$Q^\top \Sigma^{-1 / 2} B^\top X^\top\left(I_{n}-P_{H_{0}}\right) X B \Sigma^{-1 / 2} Q=U^\top \operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{r}\right) U$,
where the right-hand side is the eigenvalue decomposition of the left
and $U$ is an $r \times r$ orthogonal matrix. We have

$$\begin{aligned}
        Q^\top \Sigma^{-1 / 2} B^\top X^\top\left(I_{n}-P_{H_{0}}\right) X B \Sigma^{-1 / 2} Q & =U^\top \operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{r}\right) U \\
        & \leq U^\top \operatorname{diag}\left(\lambda_{01}, \ldots, \lambda_{0 r}\right) U \\
        & =U^\top Q_{0}^\top \Sigma^{-1 / 2} B^\top X^\top\left(I_{n}-P_{H_{0}}\right) X B \Sigma^{-1 / 2} Q_{0} U .
    \end{aligned}$$

So

$$\begin{aligned}
        M= & U^\top Q_{0}^\top \Sigma^{-1 / 2} B^\top X^\top\left(I_{n}-P_{H_{0}}\right) X B \Sigma^{-1 / 2} Q_{0} U \\
        & -Q^\top \Sigma^{-1 / 2} B^\top X^\top\left(I_{n}-P_{H_{0}}\right) X B \Sigma^{-1 / 2} Q
    \end{aligned}$$

is a semi-definite matrix. Generate
$K \sim W_{k}\left(0, I_{r}, M\right)$, which is a $r \times r$ random
matrix of non-central Wishart distribution with zero degree of freedom
and non-central parameter $M$ independent to $(X, Y)$. Note that $K$ is
semi-definite almost surely. From the property of non-central Wishart
distribution from Muirhead (2009) [@muirhead2009aspects], we know that

$$H_{\Sigma^{-1 / 2} Q}+K \stackrel{d}{=} H_{\Sigma^{-1 / 2} Q_{0} U} \stackrel{d}{=} H_{\Sigma^{-1 / 2} Q_{0}}$$

Since
$G_{\Sigma^{-1 / 2} Q} \stackrel{d}{=} G_{\Sigma^{-1 / 2} Q_{0}}, G_{\Sigma^{-1 / 2} Q}$
is independent to $H_{\Sigma^{-1 / 2} Q}+K$ and
$G_{\Sigma^{-1 / 2} Q_{0}}$ is independent to
$H_{\Sigma^{-1 / 2} Q_{0}}$, we have

$$\frac{\left|G_{\Sigma^{-1 / 2} Q}\right|}{\left|G_{\Sigma^{-1 / 2} Q}+H_{\Sigma^{-1 / 2} Q}+K\right|} \stackrel{d}{=} \frac{\left|G_{\Sigma^{-1 / 2} Q_{0}}\right|}{\left|G_{\Sigma^{-1 / 2} Q_{0}}+H_{\Sigma^{-1 / 2} Q_{0}}\right|}$$

Thus, we have

$$\begin{aligned}
    P\left(\lambda_{\Sigma^{-1 / 2} Q}<c\right) & =P\left(\frac{\mid G_{\Sigma^{-1 / 2} Q}}{\left|G_{\Sigma^{-1 / 2} Q}+H_{\Sigma^{-1 / 2} Q}\right|}<c\right) \\
        & \leq P\left(\frac{\mid G_{\Sigma^{-1 / 2} Q}}{\left|G_{\Sigma^{-1 / 2} Q}+H_{\Sigma^{-1 / 2} Q}+K\right|}<c\right) \\
        & =P\left(\frac{\mid G_{\Sigma^{-1 / 2} Q_{0}}}{\left|G_{\Sigma^{-1 / 2} Q}+H_{\Sigma^{-1 / 2} Q_{0}}\right|}<c\right)=P\left(\lambda_{\Sigma^{-1 / 2} Q_{0}}<c\right) .
\end{aligned}$$

Thus, $\Sigma^{-1 / 2} Q_{0}$ is the optimal $p \times r$ projection
matrix.

Proof of Theorem [1](#thm1){reference-type="ref" reference="thm1"}(A2)
For any $p \times r$ projection matrix $Q$ and $r>m$ and from Theorem
[1](#thm1){reference-type="ref" reference="thm1"}(A2),
$P_{0}=\Sigma^{-1 / 2} Q_{0}$ is the optimal $p \times r$ projection
direction matrix, where $Q_{0}=$ $\left(R_{1}, \cdots, R_{r}\right)$ and
$R_{i}$ is defined as in Theorem [1](#thm1){reference-type="ref"
reference="thm1"} . Note that
$P_{0}=\left(\begin{array}{ll}P_{01} & P_{02}\end{array}\right)$, where
$P_{01}=\Sigma^{-1 / 2}\left(R_{1}, \cdots, R_{m}\right)$ is the optimal
$p \times m$ projection direction matrix and $P_{02}=$
$\Sigma^{-1 / 2}\left(R_{m+1}, \cdots, R_{r}\right)$.

The LRT statistic for the $P_{0}$ projection test can be decomposed
correspondingly as follows:
$$\dfrac{|G_{P_0}|}{|G_{P_0}+H_{P_0}|}=\dfrac{|G_{22}|}{|G_{22}+H_{22}|}\dfrac{\left|\begin{pmatrix}
    G_{11} & G_{12}\\
   G_{21} & G_{22}
\end{pmatrix}\right|/|G_{22}|}{\left|\left(\begin{array}{cc}
    G_{11}+H_{11} & G_{12}+H_{12} \\
    G_{21}+H_{21} & G_{22}+H_{22}
    \end{array}\right)\right| /\left|G_{22}+H_{22}\right|}$$ where
$G_{11}, G_{11}+H_{11}$ are of dimension $m \times m$, and
$G_{22}, G_{22}+H_{22}$ are of dimension $(r-m) \times$ $(r-m) . G_{22}$
and $G_{22}+H_{22}$ can be seen as the residual sums of squares under
$H_{1}$ and $H_{0}$ for $P_{02}$ projection test. It is easy to show
that $\left(G_{22}, H_{22}\right)$ has the same distribution under
$H_{1}$ and $H_{0}$. In fact, under $H_{1}$ and
$H_{0}, G_{22} \sim W_{r-m}\left(n-d, I_{r-m}\right), H_{22} \sim W_{r-m}\left(m, I_{m}\right)$
and they are independent to each other. So
$\left|G_{22}\right| /\left|G_{22}+H_{22}\right|$ has the same
distribution under $H_{1}$ and $H_{0}$. Moreover,
$\left|\left(\begin{array}{ll}G_{11} & G_{12} \\ G_{21} & G_{22}\end{array}\right)\right| /\left|G_{22}\right|=\left|G_{11}-G_{12} G_{22}^{-1} G_{21}\right|$
can be seen as the generalized residual sum of squares of $Y P_{01}$ on
$Y P_{02}$ and $X$.
$W_{1}=G_{11}-G_{12} G_{22}^{-1} G_{21} \sim W_{m}\left(n-(r-m)-d, I_{m}\right)$.

$$\left|\left(\begin{array}{cc}
    G_{11}+H_{11} & G_{12}+H_{12} \\
    G_{21}+H_{21} & G_{22}+H_{22}
    \end{array}\right)\right| /\left|G_{22}+H_{22}\right|=\mid G_{11}+H_{11}-\left(G_{12}+H_{12}\right)\left(G_{22}+H_{22}\right)^{-1}\left(G_{21}+H_{21}\right)$$

can be seen as the generalized residual sum of squares of $Y P_{01}$ on
$Y P_{02}$ and $X A_{1}^\top$. $W_{0}=$
$G_{11}+H_{11}-\left(G_{12}+H_{12}\right)\left(G_{22}+H_{22}\right)^{-1}\left(G_{21}+H_{21}\right)=W_{1}+W_{2}$,
where
$W_{2} \sim W_{m}\left(m, I_{m}, P_{01}^\top B^\top X^\top\left(I_{n}-\right.\right.$
$\left.\left.P_{H_{0}}\right) X B \Sigma^{-1 / 2}\right)$ and $W_{2}$ is
independent to $W_{1}$.

In sum,

$$\frac{\left|\left(\begin{array}{ll}
    G_{11} & G_{12} \\
    G_{21} & G_{22}
    \end{array}\right)\right| /\left|G_{22}\right|}{\left|\left(\begin{array}{ll}
    G_{11}+H_{11} & G_{12}+H_{12} \\
    G_{21}+H_{21} & G_{22}+H_{22}
    \end{array}\right)\right| /\left|G_{22}+H_{22}\right|} \stackrel{d}{=} \frac{\left|W_{1}\right|}{\left|W_{1}+W_{2}\right|},$$

where
$W_{1} \sim W_{m}\left(n-(r-m)-d, I_{m}\right), W_{2} \sim W_{m}\left(m, I_{m}, P_{01}^\top B^\top X^\top\left(I_{n}-P_{H_{0}}\right) X B P_{01}\right)$,
and $W_{1}$ is dependent to $W_{2}$.

Note that the distribution of
$\frac{\left|W_{1}\right|}{\left|W_{1}+W_{2}\right|}$ does not depend on
$Y P_{02}$. So it is independent to
$\frac{\left|G_{22}\right|}{\left|G_{22}+H_{22}\right|}$. The proof here
is similar to the proof in Lemma 8.4.3, Lemma 8.4.4 in Anderson (2003)
[@anderson2003introduction], which decomposes the test statistic into
the product of a series of independent Beta distributed variables.

Since $\frac{\left|G_{22}\right|}{\left|G_{22}+H_{22}\right|}$ has the
same distribution under $H_{0}$ and $H_{1}$, the test statistic

$$\frac{\left|\left(\begin{array}{cc}
    G_{11} & G_{12} \\
    G_{21} & G_{22}
    \end{array}\right)\right| /\left|G_{22}\right|}{\left|\left(\begin{array}{cc}
    G_{11}+H_{11} & G_{12}+H_{12} \\
    G_{21}+H_{21} & G_{22}+H_{22}
    \end{array}\right)\right| /\left|G_{22}+H_{22}\right|}$$

is better than

$$\frac{\left|G_{22}\right|}{\left|G_{22}+H_{22}\right|} \frac{\left|\left(\begin{array}{ll}
    G_{11} & G_{12} \\
    G_{21} & G_{22}
    \end{array}\right)\right| /\left|G_{22}\right|}{\left|\left(\begin{array}{ll}
    G_{11}+H_{11} & G_{12}+H_{12} \\
    G_{21}+H_{21} & G_{22}+H_{22}
    \end{array}\right)\right| /\left|G_{22}+H_{22}\right|} .$$

Moreover,
$\frac{\left|\left(\begin{array}{ll}G_{11} & G_{12} \\ G_{21} & G_{22}\end{array}\right)\right| /\left|G_{22}\right|}{\left|\left(\begin{array}{ll}G_{11}+H_{11} & G_{12}+H_{12} \\ G_{21}+H_{21} & G_{22}+H_{22}\end{array}\right)\right| /\left|G_{22}+H_{22}\right|}$
has the same distribution with $\Lambda_{P_{01}}$ only with sample size
reduced from $n$ to $n-(r-m)$.

Thus, $P_{01}$ projection test is powerful than
$\left(\begin{array}{ll}P_{01} & P_{02}\end{array}\right)$ projection
test. So the optimal $m$ dimension projection test is better than the
optimal $r$ dimension projection test with $r>m$, which means that the
dimension for the optimal projection direction is less than or equal to
$m$ under normality assumption.

Proof of Theorem [1](#thm1){reference-type="ref" reference="thm1"}(A3).
Theorem [1](#thm1){reference-type="ref" reference="thm1"}(A3) is just a
special case of Theorem [1](#thm1){reference-type="ref"
reference="thm1"}(A1) and can be proved using the same techniques.

Proof of Theorem [1](#thm1){reference-type="ref" reference="thm1"} (B).
It is sufficient to prove that the test statistic
$\Lambda=\frac{\left|G_{P}\right|}{\left|G_{P}+H_{P}\right|}$ follows
the same asymptotic distribution in both Gaussian and non-Gaussian
multivariate linear models, thus the asymptotic results in the Gaussian
case also hold in the non-Gaussian cases.

From Conditions 1, 2, and 3 of Theorem [1](#thm1){reference-type="ref"
reference="thm1"}, we have
$\frac{G_{P}}{n} \stackrel{p}{\rightarrow} P^\top \Sigma P$ with
convergence rate $O_{p}\left(n^{-1 / 2}\right)$. So with the delta
method, we have

$$\begin{aligned}
        -\log (\Lambda) & =-\log \frac{\left|G_{P}\right|}{\left|G_{P}+H_{P}\right|}=\log \left(\left|I_{k}+H_{P} G_{P}^{-1}\right|\right) \\
        & =\log \operatorname{det}\left(I_{k}+H_{P}\left(P^\top \Sigma P\right)^{-1} / n\right)\left(1+O_{p}\left(n^{-1 / 2}\right)\right)
    \end{aligned}$$

Hence the theorem holds if $H_{P}$ follows the same asymptotic
distribution under both the Gaussian and non-Gaussian cases.

For any projection direction matrix $P$ of dimension $p \times k$, we
have $\hat{B}_{I_{1}}^*=\hat{B}^*=$
$\left(X^\top X\right)^{-1} X^\top Y P$, and
$\hat{B}_{H_{0}}^*=A_{1}\left(A_{1}^\top X^\top X A_{1}\right)^{-1} A_{1}^\top X^\top Y P$,
where $A_{1}=\left(A_{0}^\top\right)^{\perp}$. Hence,

$$d_{B}=\hat{B_{H_{1}}^*}-\hat{B_{H_{0}}^*}=\left(\left(X^\top X\right)^{-1}-A_{1}\left(A_{1}^\top X^\top X A_{1}\right)^{-1} A_{1}^\top\right) X^\top Y P$$

Then from Conditions 1, 2 and 3 of Theorem
[1](#thm1){reference-type="ref" reference="thm1"}, we have
$d_{B}=\hat{B_{H_{1}}^*}-\hat{B_{H_{0}}^*}$ has a limiting normal
distribution with

$$\begin{aligned}
        \mathbb{E} d_{B} & =\left(\left(X^\top X\right)^{-1}-A_{1}\left(A_{1} X^\top X A_{1}\right)^{-1} A_{1}^\top\right) X^\top X B P \\
        & =\left(I_{d}-A_{1}\left(A_{1}^\top X^\top X A_{1}\right)^{-1} A_{1}^\top X^\top X\right) B P
    \end{aligned}$$

and
$\sqrt{n} \operatorname{vec}\left(d_{B}-\mathbb{E} d_{B}\right) \stackrel{d}{\rightarrow} N_{d k}(0, V)$,
where

$$\begin{aligned}
        V & =\lim _{n \rightarrow \infty}\left(\left(X^\top X / n\right)^{-1}-A_{1}\left(A_{1}^\top X^\top X A_{1} / n\right)^{-1} A_{1}^\top\right) \otimes\left(P^\top \Sigma P\right) \\
        & =\left(M_{X}^{-1}-A_{1}\left(A_{1}^\top M_{X} A_{1}\right)^{-1} A_{1}^\top\right) \otimes\left(P^\top \Sigma P\right) .
    \end{aligned}$$

Hence we have
$H_{P}=d_{B}^\top X^\top X d_{B} \stackrel{d}{\rightarrow} W_{k}\left(m, I_{k}, C\right)$,
where $C$ is the non-central parameter for the Wishart distribution and
$C=P^\top B^\top X^\top\left(I_{n}-X A_{1}\left(A_{1}^\top X^\top X A_{1}\right)^{-1} A_{1}^\top X^\top\right) X B P$.

Hence we prove that $H_{P}$ follows the same asymptotic distribution in
both the Gaussian and non-Gaussian cases under Conditions 1, 2, and 3. ◻
:::

# Estimation of Optimal Projection Direction

Theorem [1](#thm1){reference-type="ref" reference="thm1"} explains how
the optimal projection direction is determined by $\Sigma$ and $B$. We
shall devise an estimation technique for the optimal direction. We
recommend setting $k = m$ for the mentioned hypothesis above with small
$m$, as we did for the one-sample and two-sample mean testing issues,
and the optimal projection direction is $\Sigma^{-1}B^\top A_0^\top$
according to Theorems $1(A2)$ and $(A3)$. As seen below, the ideal
projection matrix can be computed column by column.
$$P_{0,k} = \Sigma^{-1}B^\top A_{0,k}^\top\quad \forall k=1,2,\cdots,m$$
where $A_{0,k}$ is the $k^{th}$ row of $A_0$ and $P_{0,k}$ is the
$k^{th}$ column of $P_0$. So without loss of generality we can assume,
$m=1$. Now we have to estimate the optimal projection direction
$$P_0 = \Sigma^{-1}B^\top A_{0}^\top$$ So basically the estimate can be
obtained as, $$\hat{P} = \hat{\Sigma}^{-1} \hat{B}^\top A_0^\top$$
However, when $p$ is larger than $m$, then $\hat{\Sigma}$ may become
singular. To deal with this problem, we will use an shrinkage estimator
like ridge regression to find,
$$\hat{P}_{\text{ridge}} = (\lambda_0I_p + \hat{\Sigma})^{-1}  \hat{B}^\top A_0^\top$$
where $\lambda_0>0$. Now to find the estimate of and $\hat{B}$, we use
the ordinary least squares i.e.,
$$\hat{B}_{\text{LS}} = (X^\top X)^{-1}X^\top Y$$ and estimate
$\hat{\Sigma}$ by the sample covariance of the residuals i.e.,
$$\hat{\Sigma}_{\text{LS}} = \dfrac{1}{n-d}(Y-X\hat{B}_{\text{LS}})^\top(Y-X\hat{B}_{\text{LS}})$$
Although having prior knowledge of $B$ and $\Sigma$, using specific
estimators for them can increase the power of the test, we advise using
general estimators of $B$ and $\Sigma$ in this article to reflect the
general case. Now from $\hat{B}_{\text{LS}}$ and
$\hat{\Sigma}_{\text{LS}}$ we have,
$$\hat{P}_{\text{LS}} = (\lambda_0I_p + \hat{\Sigma}_{\text{LS}})^{-1} \hat{B}_{\text{LS}}^\top A_0^\top$$

# U-Projection test

Consider the general multivariate linear model
([\[model\]](#model){reference-type="ref" reference="model"}) and the
hypothesis $A_0B=0$. From Theorem [1](#thm1){reference-type="ref"
reference="thm1"}, $P_0 = \Sigma^{-1}B^\top A_{0}^\top$ is the
asymptotic optimal projection matrix of dimension $p\times m$. The null
hypothesis of the projection test with direction $P_0$ is
$H_{0P} : A_0BP_0 = A_0B\Sigma^{-1}B^\top A_{0}^\top = 0$. Since
$A_0B\Sigma^{-1}B^\top A_{0}^\top$ is positive semi-definite, it is
equivalent to the test
$H_{0P} : tr(A_0BP_0) = tr(A_0B\Sigma^{-1}B^\top A_{0}^\top) = 0$.
Following the sample-splitting test procedure, we split the data
$(X, Y)$ into two parts $(X_1, Y_1)$ and $(X_2, Y_2)$. The sample sizes
of $X_1, Y_1$ are $k$, and the sample sizes of $X_2, Y_2$ are $n-k$. The
first sample is used an estimator $\hat{P}_{(X_1,Y_1)}$ of the
$m$-dimensional projection and sample-splitting test statistic on the
projected second sample is calculated as follows:
$$tr\left(A_0\hat{B}_{(X_2,Y_2)}\hat{P}_{(X_1,Y_1)}\right) = tr\left( A_0(X_2^\top X_2)^{-1}X_2^\top Y_2\hat{P}_{(X_1,Y_1)} \right)$$
So the subsequent U-projection statistic can be written as,
$$U_p = \dfrac{1}{|\Gamma|} \sum_{\gamma\in\Gamma} tr\left( A_0(X_{-\gamma}^\top X_{-\gamma})^{-1}X_{-\gamma}^\top Y_{-\gamma}\hat{P}_{(X_{\gamma},Y_{\gamma})} \right) \tag{3}$$
where
$\Gamma = \{\gamma | \gamma \subset \{1,2,\cdots,n\},|\gamma| = k$,
rank$(X_{\gamma})> d$, rank($X_{-\gamma} ) \geq d\}$, and
$X_\gamma , Y_\gamma , X_{-\gamma}, Y_{-\gamma}$ are subsamples of $X$
and $Y$ with and without index $\gamma$, respectively. Now we can have
the asymptotic distribution of $U_p$ using the following theorem:

::: {#thm2 .theorem}
**Theorem 2**. *Consider $H_{0}: A_{0} B=0$ under model (1). Let
$$h_{P}\left(Z_{1,1}, \ldots, Z_{k+d, 1} ; Z_{1,2}, \ldots, Z_{k+d, 2}\right)
= \frac{1}{|\Gamma|} \sum_{\gamma \in \Gamma} \operatorname{diag}\left(A_{0}\left(Z_{-\gamma, 1}^\top Z_{-\gamma, 1}\right)^{-1} Z_{-\gamma, 1}^\top Z_{-\gamma, 2} \hat{P}_{\left(Z_{\gamma, 1}, Z_{\gamma, 2}\right)}\right)$$
where
$-\gamma=\{1, \ldots, k+d\} \backslash \gamma, \Gamma=\{\gamma, \gamma \subset\{1,2, \ldots, k+$
$\left.d\},|\gamma|=k, \operatorname{rank}\left(Z_{\gamma, 1}\right)>d, \operatorname{rank}\left(Z_{-\gamma, 1}\right)=d\right\}$
; $Z_{\gamma, i}, Z_{-\gamma, i}$ are samples of $Z_{\ldots, i}$ with
and without index $\gamma$ for $i=1,2$; and
$\hat{P}_{\left(Z_{\gamma, 1}, Z_{\gamma, 2}\right)}$ is the projection
direction estimated with $\left(Z_{\gamma, 1}, Z_{\gamma, 2}\right)$.
Then with the conditions that $k$ is fixed and that
$\mathbb{E}[ \operatorname{tr}\left(\operatorname{cov}\left(h_{P}\left(Z_{1,1}, \ldots, Z_{k+d, 1} ; Z_{1,2}, \ldots, Z_{k+d, 2}\right)\right)\right)]<C_{0}$
for some fixed $C_{0}>0$ and
$\left(Z_{i, 1}, Z_{i, 2}\right), i=1, \ldots, k+d$ i.i.d. from the
distribution of $(X, Y)$, we have asymptotic normality of the general
U-projection statistic (3):
$$\sqrt{n}\left(U_{P}-\mathbb{E} [U_{P}]\right) \stackrel{d}{\rightarrow} N\left(0,(k+d)^{2} 1_{m}^\top \Xi_{1} 1_{m}\right)$$
where $1_{m}$ is a vector full of one with length $m$, and
$$\Xi_{1}=\operatorname{cov}\left(h_{P}\left(Z_{1,1}, Z_{2,1}, \ldots, Z_{k+d, 1} ; Z_{1,2}, \ldots, Z_{k+d, 2}\right),
h_{P}\left(Z_{1,1}, Z_{2,1}^{\prime}, \ldots, Z_{k+d, 1}^{\prime} ; Z_{1,2}, Z_{2,2}^{\prime}, \ldots, Z_{k+d, 2}^{\prime}\right) \right)$$
for $\left(Z_{i, 1}, Z_{i, 2}\right), i=1, \ldots, k+d$ and
$\left(Z_{i, 1}^{\prime}, Z_{i, 2}^{\prime}\right), i=2, \ldots, k+d$
i.i.d. from the distribution of $(X, Y)$. Furthermore, under the null
hypothesis, we have $\mathbb{E} [U_{P}]=0$. *
:::

::: proof
*Proof.* Let
$$M_{P}=\frac{1}{|\Gamma|} \sum_{\gamma \in \Gamma} \operatorname{diag}\left(A_{0}\left(X_{-\gamma}^{T} X_{-\gamma}\right)^{-1} X_{-\gamma}^{T} Y_{-\gamma} \hat{P}_{\left(X_{\gamma}, Y_{\gamma}\right)}\right)$$
where
$\Gamma=\left\{\gamma|\gamma \subset\{1,2, \cdots, n\},| \gamma \mid=k, \operatorname{rank}\left(X_{\gamma}\right)>d, \operatorname{rank}\left(X_{-\gamma}\right) \geq d\right\}$,
and $X_{\gamma}, Y_{\gamma}$, $X_{-\gamma} Y_{-\gamma}$ are sub-samples
of $X$ and $Y$ with and without index $\gamma$ respectively. Then we
have $U_{P}=1_{m}^{T} M_{P}$. $M_{P}$ can be rewritten as a U-statistic
from the kernel $h_{P}$. From the asymptotic normality result of
$\mathrm{U}$-statistics, we have
$$\sqrt{n}\left(M_{P}-\mathbb{E} [M_{P}]\right) \stackrel{d}{\rightarrow} N\left(0,(k+d)^{2} \Xi_{1}\right)$$
where $\Xi_{1}$ is as defined in the theorem. Hence the asymptotic
result of $U_{P}$ follows from the fact that $U_{P}=1_{m}^{T} M_{P}$.
Moreover, it is easy to get $\mathbb{E}_{H_{0}} [U_{P}]=0$ under the
null hypothesis. This completes the proof of the theorem.◻
:::

Although the asymptotic distribution of the U-projection statistic can
be obtained under certain assumptions, as shown in Theorem
[2](#thm2){reference-type="ref" reference="thm2"}, it may not be very
useful here for a number of reasons:

-   If the sample size $n$ is small, the asymptotic distribution will be
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
    
**Input** $n\times d$ dimensional matrix $X$, $n\times p$ dimensional matrix $Y$, $m\times d$ dimensional matrix $A_0$, $d<k\leq n-d$, and randomization times $N$.
    
Calculate the $U$-projection statistic $U_0$ on the original dataset $X$ and $Y$.<br>
Let $i=0$<br>
**For** (j=1:N)$\{$<br>
&emsp; Do randomization, get $X'$ and $Y'$<br>
&emsp; Calculate the $U$-projection statistic $U'$ on the randomized dataset $X'$ and $Y'$.<br>
&emsp; **If** $(U'\geq U)\{$ <br>
&emsp; &emsp; $i=i+1$<br>
&emsp; $\}$<br>
$\}$<br>
Calculate p-value by $p=\dfrac{i}{N}$

## Examples

### One Sample Mean Testing:

For one-sample mean testing, suppose we have $n$ samples $Y$ from a
$p$-dimensional multivariate distribution with mean $\mu$ and covariance
$\Sigma$, and we want to test whether $\mu = 0$. Now for this setting we
can rewrite our model as, $$Y=1_nB + E$$ where $1_n$ is a column vector
full of $1$ of length $n$, $B = \mu^\top$, $E$ is the $n \times p$
random error matrix, and the null hypothesis becomes $$H_0 : B = 0$$
i.e., $A_0 =1$ and is of rank $1$ in this problem. By Theorem
[1](#thm1){reference-type="ref" reference="thm1"}, the optimal
projection dimension is $1$ and the optimal projection direction is
$\Sigma^{-1}\mu$. Now we can find the estimates of $\Sigma$ and $\mu$
through least squares and subsequently estimate the optimal projection
direction. So the estimated projection matrix becomes,
$$\hat{P}=\left(\lambda_{0} I_{p}+\hat{\Sigma}\right)^{-1} \hat{B}_{\mathrm{LS}}^\top 1_n^\top=\left(\lambda_{0} I_{p}+\hat{\Sigma}\right)^{-1}\bar{Y}$$
Thus, the U-projection statistic becomes
$$U_{\mathrm{P}}=  \dfrac{1}{\binom{n}{k}}\sum_{\gamma \in \Gamma}\text{tr}\left(\bar{Y}_{-\gamma}^\top \times\left\{\lambda_{0} I_{p}+S_{Y_{\gamma}}\right\}^{-1} \times \bar{Y}_{\gamma}\right)$$
where for $k$ is the number of independent samples from $Y$ to estimate
the projection direction, $\Gamma$ is collections of all subsets of
$\left\{1, \ldots, n\right\}$ with size $k$, and $Y_{\gamma}$ and
$Y_{-\gamma}$ are subsets of $Y$ with and without index $\gamma$
correspondingly.

### Testing of Significance of Predictors

Consider the general linear model (1). Testing
the significance of a particular predictor, $X^{(i)}$, is of interest.
In other words, take $$H_0:A_0B=0$$ where $A_0=e^\top_i$ and $e_i$ is a
column vector in $d$ dimensions, with the $i^{th}$ element being one and
the remaining elements being zero. The above hypothesis testing problem
is also testing for the conditional independence of $Y$ and $X^{(i)}$
given other $(d-1)$ predictors under the normality assumption. This is
because, in the normal case, testing whether the regression coefficients
equal zero or not is equivalent to testing of conditional independence.
According to Theorem [1](#thm1){reference-type="ref" reference="thm1"},
the ideal projection dimension is 1 and the optimal projection direction
is $\Sigma^{-1}B^\top e_i$. Now using the previously demonstrated
estimation method, we can estimate the projection matrix by,
$$\hat{P}=\left(\lambda_{0} I_{p}+\hat{\Sigma}\right)^{-1} \hat{B}_{\mathrm{LS}}^\top e_i=\left(\lambda_{0} I_{p}+\hat{\Sigma}\right)^{-1}\hat{B}_i$$
Thus, the U-projection statistic becomes
$$U_{\mathrm{P}}=  \dfrac{1}{\binom{n}{k}}\sum_{\gamma \in \Gamma}\text{tr}\left(\hat{B}_i^\top(Y_{-\gamma}) \times\left\{\lambda_{0} I_{p}+S_{Y_{\gamma}}\right\}^{-1} \times \hat{B}_i(Y_{\gamma}) \right)$$
where for $k$ is the number of independent samples from $Y$ to estimate
the projection direction, $\Gamma$ is collections of all subsets of
$\left\{1, \ldots, n\right\}$ with size $k$, and $Y_{\gamma}$ and
$Y_{-\gamma}$ are subsets of $Y$ with and without index $\gamma$
correspondingly.

# Future Direction: Multiple Sample Mean Testing

For general unbalanced multiple-sample mean testing problem, suppose for
$k = 1, 2,\cdots, K$, we have $n_k$ samples $Y_k$ from a $p$-dimensional
multivariate distribution $F_k$ with mean $\mu_k$ and covariance
$\Sigma$, and we want to test whether $\mu_1=\mu_2=\cdots=\mu_k$. The
problem can be reformulated as $$\begin{pmatrix}
    Y_1\\
    \cdots\\
    Y_k
\end{pmatrix} = \begin{pmatrix}
    1_{n_1} & 0 & \cdots\\
    \cdots & \cdots & \cdots\\
    0 & \cdots & 1_{n_k}
\end{pmatrix} B + E$$ where $B=(\mu_1,\mu_2,\cdots, \mu_k)^\top$, and
$E$ is the $\left( \sum_{i=1}^k n_i\right)\times p$ random error matrix.
We want to test $$H_0: A_0B = 0$$ where
$A_0 = \{a_{ij}\}_{1\leq i\leq K-1,1\leq j\leq K}$ with $a_{ii}=K-1$ and
$a_{ij}=-1,\forall i\ne j$. In this problem, $A_0$ is of rank $K-1$. By
Theorem [1](#thm1){reference-type="ref" reference="thm1"}, the optimal
projection dimension is less than or equal to $K-1$ and the optimal
$p\times(K-1)$ projection direction matrix is
$K\Sigma^{-1}(\mu_1-\bar{\mu},\cdots,\mu_{k-1}-\bar{\mu})$, where
$\bar{\mu} = \frac{1}{k}\sum_{k=1}^K \mu_k$. In particular,when $K=2$,
the optimal projection dimension is $1$ and the optimal projection
direction is $\Sigma^{-1}(\mu_1 - \mu_2)$. Now using the previously
demonstrated estimation method, we can estimate the projection matrix
by,
$$\hat{P}_{\mathrm{LS}}=\left(\lambda_{0} I_{p}+\hat{\Sigma}\right)^{-1} \hat{B}_{\mathrm{LS}}^\top A_{0}^\top=\left(\lambda_{0} I_{p}+\hat{\Sigma}\right)^{-1}\left(\bar{Y}_{1}-\bar{Y}_{2}\right)$$
Thus, the U-projection statistic becomes
$$U_{\mathrm{LS}}=  \dfrac{1}{\binom{n_1}{k_1}\binom{n_2}{k_2}} \sum_{\gamma_{1} \in \Gamma_{1}} \sum_{\gamma_{2} \in \Gamma_{2}}\left(\bar{Y}_{-\gamma_{1}, 1}^\top-\bar{Y}_{-\gamma_{2}, 2}^\top\right) \times\left\{\lambda_{0} I_{p}+\frac{\left(k_{1}-1\right) S_{Y_{\gamma_{1}, 1}}+\left(k_{2}-1\right) S_{Y_{\gamma_{2}, 2}}}{k_{1}+k_{2}-2}\right\}^{-1} \times\left(\bar{Y}_{\gamma_{1}, 1}-\bar{Y}_{\gamma_{2}, 2}\right)$$
where for $i=1,2, k_{i}$ is the number of independent samples from
$Y_{i}$ to estimate the projection direction, $\Gamma_{i}$ is
collections of all subsets of $\left\{1, \ldots, n_{i}\right\}$ with
size $k_{i}$, and $Y_{\gamma, i}$ and $Y_{-\gamma, i}$ are subsets of
$Y_{i}$ with and without index $\gamma$ correspondingly. Now consider
the following condition and theorems:

**Condition 4:**
$E_{i}=\Gamma Z_{i} \quad \text { for } i=1, \ldots, n$, where $\Gamma$
is a $p \times t$ matrix with some $t \geq p$ such that
$\Gamma \Gamma^{T}=$ $\Sigma$, and
$Z_{i}=\left(Z_{i, 1}, \ldots, Z_{i, t}\right)$ are $t$-variate
independent and identically distributed random vectors satisfying
$\mathbb{E}\left(Z_{i}\right)=0$,
$\operatorname{var}\left(Z_{i}\right)=I_{t}, \mathbb{E}\left(Z_{i, k}^{3}\right)=0, \mathbb{E}\left(Z_{i, k}^{6}\right)$
is uniformly bounded, and
$$\mathbb{E}\left(Z_{i, l_{1}}^{\alpha_{1}} Z_{i, l_{2}}^{\alpha_{2}} \cdots Z_{i, l_{s}}^{\alpha_{s}}\right)=\mathbb{E}\left(Z_{i, l_{1}}^{\alpha_{1}}\right) \mathbb{E}\left(Z_{i, l_{2}}^{\alpha_{2}}\right) \cdots \mathbb{E}\left(Z_{i, l_{s}}^{\alpha_{s}}\right)$$
for a positive integer $s$ such that $\sum_{l=1}^{s} \alpha_{l} \leq 8$
and $l_{1} \neq l_{2} \neq$ $\cdots \neq l_{s}$.

**Condition 5:** $\left\|A_{0} B\right\|_{\mathrm{F}}^{2}=o(p / n)$.

**Condition C1:** There exists a uniformly bounded positive integer
$q<p$ such that
$\frac{\sqrt{n} \lambda_{q}}{\operatorname{tr}(\Sigma)} \rightarrow \infty$
and $\lambda_{q+1}$ is uniformly bounded from above. The smallest
eigenvalue $\lambda_{p}$ is uniformly bounded from below.

**Condition C2:**
$\operatorname{tr}\left(\Sigma^{4}\right)=o\left(\operatorname{tr}^{2}\left(\Sigma^{2}\right)\right)$.

::: {#thm3 .theorem}
**Theorem 3**. *Suppose the covariance $\Sigma$ satisfies Condition C1.
Under high-dimensional setting
$n_{1}+n_{2}=o(\operatorname{tr}(\Sigma))$, Conditions $3,4,5$, and
conditions that $k_{i} / n_{i} \rightarrow \gamma_{i} \in(0,1)$ for
$i=1,2$, $n_{1} /\left(n_{1}+n_{2}\right) \rightarrow \kappa \in(0,1)$,
then $U_{\mathrm{LS}}$ has an asymptotically normal distribution with a
uniformly bounded positive $\lambda_{0}$. More specifically, we have*

*$$\frac{\lambda_{0}\left(U_{\mathrm{LS}}-\mathbb{E}\left(U_{\mathrm{LS}}\right)\right)}{\sigma_{n}} \stackrel{d}{\rightarrow} N(0,1), \quad \text { and } \quad \lambda_{0} \mathbb{E}\left(U_{\mathrm{LS}}\right)-\left\|W_{q+1} \mu_{n}\right\|_{2}^{2}=o\left(\left\|W_{q+1} \mu_{n}\right\|_{2}^{2}\right),$$
where $q$ is the number of divergent eigenvalues of $\Sigma$ as
specified in Condition C1, $W_{q+1}$ is the projection matrix onto the
linear span of eigenspaces of $\Sigma$ corresponding to the smallest
$p-q$ eigenvalues,
$\lambda_{q+1}, \lambda_{q+2}, \ldots, \lambda_{p}, \mu_{n}$ is the mean
difference of the two populations, and
$\sigma_{n}^{2}=\left(\frac{2}{n_{1}^{2}}+\frac{2}{n_{2}^{2}}+\frac{4}{n_{1} n_{2}}\right) \sum_{i=q+1}^{p} \lambda_{i}^{2}$.
*
:::

::: {#thm4 .theorem}
**Theorem 4**. *Suppose the covariance $\Sigma$ satisfies Condition C2.
Under high-dimensional setting
$n_{1}+n_{2}=O(\operatorname{tr}(\Sigma))$, Conditions 3-5, and
conditions that $k_{i} / n_{i} \rightarrow \gamma_{i} \in[0,1)$ for
$i=1,2$, $n_{1} /\left(n_{1}+n_{2}\right) \rightarrow \kappa \in(0,1)$,
and*

*$$\left(k_{1}+k_{2}\right)=o\left(\operatorname{tr}(\Sigma) / \lambda_{\max }(\Sigma)^{2}\right)$$
where $\lambda_{\max }(\Sigma)$ is the largest eigenvalue of
$\Sigma, \lambda_{0} U_{\mathrm{LS}}$ is asymptotically normally
distributed and has the same asymptotic variance and similar expectation
with $T_{\mathrm{CQ}}$,*

*$$\begin{aligned}
& \frac{\lambda_{0}\left(U_{\mathrm{LS}}-\mathbb{E} U_{\mathrm{LS}}\right)}{\sqrt{\operatorname{var} T_{\mathrm{CQ}}}} \stackrel{d}{\rightarrow} N(0,1), \quad \text { and } \\
& \left|\mathbb{E}\left(\lambda_{0} U_{\mathrm{LS}}-T_{\mathrm{CQ}}\right)\right|=o\left(\mathbb{E} T_{\mathrm{CQ}}\right),
\end{aligned}$$ and $U_{\mathrm{LS}}$ and $T_{\mathrm{CQ}}$ have the
same asymptotic power, $\beta_{U_{\mathrm{LS}}}\left(\mu_{1}-\right.$
$\left.\mu_{2}\right)-\beta_{T_{\mathrm{CQ}}}\left(\mu_{1}-\mu_{2}\right) \rightarrow 0$.
*
:::

Above $T_{CQ}$ denote the two-sample statistic in Chen and Qin (2010)
[@chen2010two]. Now it can be observed that for high correlation
covariance $\Sigma$, Theorem [3](#thm3){reference-type="ref"
reference="thm3"} gives the asymptotic distribution of $U_{LS}$ and for
low correlation covariance $\Sigma$, Theorem
[4](#thm4){reference-type="ref" reference="thm4"} does the same. Note
that, Theorem [3](#thm3){reference-type="ref" reference="thm3"} and
[4](#thm4){reference-type="ref" reference="thm4"} provide different
formulas to calculate the asymptotic variance under different structures
of covariances. So we need to calculate the p-values using randomization
only. Now when the number of samples becomes greater i.e., $K>2$,
although the p-values can be obtained through randomization, there are
no proper documentation of the U-projection statistic. We think this is
an area which is still under research and can be developed further.

# References

<div id = "refs"></div>
