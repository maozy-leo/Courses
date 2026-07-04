# Inference in Linear Models (Small Sample Case)

For finite sample settings(i.e. with small $n$), classical asymptotic inference may not be reliable. Two widely used alternatives are:

1. **Bootstrap**
2. **Permutation Test**

## Bootstrap

Consider the linear model:

$$
Y_i = x_i' \beta + \varepsilon_i, \quad i = 1, \dots, n
$$
We are interested in inference on a linear functional:

$$
c' \beta
$$
The estimator is given by

$$
\widehat{c' \beta} = c' (X' X)^{-1} X' Y
$$

We now introduce a bootstrap method called **residual bootstrap procedure**.

### Step 1: Fit the model

Estimate:

$$
\hat{\beta}, \quad \hat{Y}_i = x_i' \hat{\beta}
$$

Residuals:

$$
\hat{\varepsilon}_i = Y_i - x_i' \hat{\beta}
$$

### Step 2: Resample residuals

Sample with replacement from:

$$
{\hat{\varepsilon}_1, \dots, \hat{\varepsilon}_n}
$$

to obtain:

$$
{\hat{\varepsilon}_1^*, \dots, \hat{\varepsilon}_n^*}
$$

### Step 3: Generate bootstrap sample

$$
Y_i^* = x_i' \hat{\beta} + \hat{\varepsilon}_i^*
$$

Construct new dataset:

$$
{(x_i, Y_i^*)}_{i=1}^n
$$

### Step 4: Refit model

Compute:

$$
\hat{\beta}^*
\quad \Rightarrow \quad c' \hat{\beta}^*
$$

### Step 5: Repeat

Repeat steps 2–4 for ( B ) times:

$$
c' \hat{\beta}^{*(1)}, \dots, c' \hat{\beta}^{*(B)}
$$

### Step 6: Construct confidence interval

Use empirical quantiles:

- Percentile CI:
  $$
  \left[ q_{\alpha/2}, ; q_{1-\alpha/2} \right]
  $$
  where $q_\cdot$ are quantiles of the bootstrap samples.

## Remarks

- Many bootstrap variants exist:
  - Residual bootstrap (used here)
  - Pairs bootstrap
  - Wild bootstrap (for heteroscedasticity)
- Works well for **confidence intervals and standard errors**
- Reliability depends on:
  - model correctness
  - residual exchangeability

# Permutation Test

We test:

$$
H_0: \beta = 0
\quad \text{vs} \quad
H_1: \beta \ne 0
$$

Under $H_0$, the relationship between $X$ and $Y$ disappears.
Thus, the pairing between $X$ and $Y$ is **exchangeable**.

### Step 1: Original data

$$
(x_1, Y_1), (x_2, Y_2), \dots, (x_n, Y_n)
$$

Compute estimator:

$$
\hat{\beta}
$$

### Step 2: Permute response

Shuffle $Y$:

$$
(Y_1, \dots, Y_n) \rightarrow (Y_1^*, \dots, Y_n^*)
$$

Construct permuted dataset:

$$
(x_i, Y_i^*)
$$

### Step 3: Recompute estimator

$$
\hat{\beta}^*
$$

### Step 4: Repeat

Repeat permutation $N$ times:

$$
\hat{\beta}^{*(1)}, \dots, \hat{\beta}^{*(N)}
$$

- In theory: $N = n!$
- In practice: use random permutations

## p-value

Using a norm-based statistic:

$$
p = \frac{1}{N} \sum_{k=1}^N
\mathbf{1}\left( ||\hat{\beta}^{*(k)}||_2 \ge ||\hat{\beta}||_2 \right)
$$

