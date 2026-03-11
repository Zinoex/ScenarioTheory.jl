# ScenarioTheory.jl

[![Build Status](https://github.com/Zinoex/ScenarioTheory.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Zinoex/ScenarioTheory.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Zinoex/ScenarioTheory.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Zinoex/ScenarioTheory.jl)
[![Fuzzy Tests](https://raw.githubusercontent.com/Seelengrab/Supposition.jl/main/badge.svg)](https://github.com/Seelengrab/Supposition.jl)

This package provides computation of various bounds of the Scenario Theory pioneered by Marco Campi and Simone Garatti [1]. Scenario theory is commonly framed as optimization [1] or compression [2]; _this package supports neither_. Instead, the purpose is, given the number of samples, the desired confidence, and the number of decision variables/support constraints/compressed size, to compute the level of violation defined as the probability of a change of compression. 

The computation of the violation is _numerically unstable_ with both large binomial coefficients and large exponents of values in the range `[0, 1]`, which are hard to represent with floats. I strongly recommend that you do _not_ try to do it manually - this is exactly why I created this package.

A goal is to also support inverse problems: given a desired confidence and (a bound on) compression cardinality, compute the required number of samples. 


## Installation

To install the package, run the following command from a Julia REPL.
```julia
(v1.11) pkg> add ScenarioTheory
```

## Usage

The package is used as follows:
```julia
using ScenarioTheory

samples = 1500
decision_variables = 30
β = 1e-6  # the confidence is 1 - β

theory = ScenarioOptimization(samples, decision_variables)

ϵ_lower, ϵ_upper = violation(theory, β)

# Output: ϵ_lower = 0.0, ϵ_upper = 0.041878994612488896
```

## Theory types

### Scenario optimization

This theory is the original based around convex optimization with chance-constraints [1]. That is, define the violation probability $V(x) = \mathbb{P}[\delta \in \Delta : x \notin X_\delta]$ where $\delta$ is a random variable with support $\Delta$ and $X_\delta$ is a convex set conditional on $\delta$. Then, the corresponding chance-constraint can be formulated as $V(x) < \epsilon$.
Now, assume a given dataset $D = \{\delta^1, \ldots, \delta^N\}$ of size $N$ sampled according to $\mathbb{P}^N$ and let
```math
    \begin{aligned}
        x^\star(D) := \mathop{arg\,min}_{x \in R^d} & f(x)
                & x \in X_{\delta^i}, \quad i = 1, \ldots, N.
    \end{aligned}
```
Then, $\mathbb{P}^N[V(x^\star(D)) < \epsilon] \geq 1 - \beta$ where $\beta$ is a function of $N$ and $d$. For a more in-depth explanation, please refer to [1].

> [!WARNING]
> The theory requires that almost surely the problem is feasible and the solution is unique. This is not checked by this code (it does not have access to the distribution nor the optimization problem); it is your responsibility to check that. 

To compute $\epsilon$ given $N$, $d$, and $\beta$, run the following:
```julia
using ScenarioTheory

samples = 1500
decision_variables = 30
β = 1e-6  # the confidence is 1 - β

theory = ScenarioOptimization(samples, decision_variables)

ϵ_lower, ϵ_upper = violation(theory, β)
```
`ϵ_lower` is always zero for `ScenarioOptimization`.


### Wait-and-judge scenario optimization

This theory is similar to vanilla scenario optimization except that, rather than relying on the number of decision variables, it relies on the number of support constraints (from $D$) [3]. A constraint is a support constraint if removing it changes the optimal value $ f(x^\star(D)) $. Let the number of support constraints be $k$.
Then, $\mathbb{P}^N[V(x^\star(D)) < \epsilon] \geq 1 - \beta$ where $\beta$ is a function of $N$ and $k$. For a more in-depth explanation, please refer to [3].

> [!WARNING]
> The theory requires that almost surely the problem is feasible, the solution is unique, and the solution with only the support constraints coincides with the solution with all constraints (non-degeneracy). This is not checked by this code (it does not have access to the distribution nor the optimization problem); it is your responsibility to check that. 

To compute $\epsilon$ given $N$, $d$, and $\beta$, run the following:
```julia
using ScenarioTheory

samples = 1500
support_constraints = 10
β = 1e-6  # the confidence is 1 - β

theory = WaitAndJudge(samples, support_constraints)

ϵ_lower, ϵ_upper = violation(theory, β)
```
`ϵ_lower` is always zero for `WaitAndJudge`.


### One-tail change of compression

This theory is a modern reformulation and generalization in terms of compression. To this end, let $D$ be a multi-set sampled according to $\mathbb{P}^N$ and assume a given compression function $c$ such that $c(D) \subset D$ for all $D$. The change of compression is defined as $\phi(D) = \mathbb{P}[c(c(D), \delta) != c(D) | D]$. Let $k$ be cardinality of $D$. Then, $\mathbb{P}^N[\phi(D) < \epsilon] \geq 1 - \beta$ where $\epsilon$ is a function of $N$, $k$, and $\beta$. For a more in-depth explanation, please refer to [3].

> [!WARNING]
> The theory requires that compression function satisfies a preference property (see [3]). This is not checked by this code (it does not have access to the compression function); it is your responsibility to check that. 

To compute $\epsilon$ given $N$, $d$, and $\beta$, run the following:
```julia
using ScenarioTheory

samples = 1500
compression = 10
β = 1e-6  # the confidence is 1 - β

theory = CompressionOneTail(samples, compression)

ϵ_lower, ϵ_upper = violation(theory, β)
```
`ϵ_lower` is always zero for `CompressionOneTail`.


### Two-tail change of compression

This theory is an extension to no only upper but also lower bounds on the violation probability.

Then, $\mathbb{P}^N[\underline{\epsilon} < \phi(D) < \overline{\epsilon}] \geq 1 - \beta$ where $\underline{\epsilon}$ and $\overline{\epsilon}$ are functions of $N$, $k$, and $\beta$. For a more in-depth explanation, please refer to [3].

> [!WARNING]
> The theory requires that compression function satisfies a preference property, the distribution has no concentrated mass, and the compression function is almost surely non-associative (see [3]). This is not checked by this code (it does not have access to the distribution nor the compression function); it is your responsibility to check that. 

To compute $\epsilon$ given $N$, $d$, and $\beta$, run the following:
```julia
using ScenarioTheory

samples = 1500
compression = 10
β = 1e-6  # the confidence is 1 - β

theory = CompressionTwoTail(samples, compression)

ϵ_lower, ϵ_upper = violation(theory, β)
```


## References

[1] Campi, M. C., & Garatti, S. (2008). The exact feasibility of randomized solutions of uncertain convex programs. SIAM Journal on Optimization, 19(3), 1211-1230.

[2] Campi, M. C., & Garatti, S. (2023). Compression, generalization and learning. Journal of Machine Learning Research, 24(339), 1-74.

[3] Campi, M. C., & Garatti, S. (2018). Wait-and-judge scenario optimization. Mathematical Programming, 167(1), 155-189.
