The problem of estimability comes from the scenario: given data and model, in what case can we get the estimators of the model?

In the following part we propose the condition of estimability of general linear model.



<u>Thm</u>

Equivalent definition of estimability for a vector $c\in\mathbb R^p$:

- $c'\hat\beta$ in unbiased for $c'\beta$, where $\hat\beta=(X'X)^{-}X'Y$ is a solution of the normal equation $X'X\hat\beta=X'Y$;
- $c'\beta$ is unique for all solutions of normal equation;
- $c\in\mathscr C(X')=\mathscr C(X'X)$;
- there exists a linear function of $Y$ which is an unbiase estimation for $c'\beta$.

<u>Proof</u>

"$1\Rightarrow2$"：
$$
\mathbb E{c'\hat\beta}=c'(X'X)^{-}X'X\beta=c'\beta,\forall\beta
$$
Then
$$
c'=c'(X'X)^{-}X'X
$$
Let $\hat\beta_1,\hat\beta_2$ be two solutions of the Normal Equation.
$$
c'\hat\beta_1-c'\hat\beta_2=c'(\hat\beta_1-\hat\beta_2)=c'(X'X)^{-}X'X(\hat\beta_1-\hat\beta_2)=c'(X'X)^{-}(X'Y-X'Y)=0
$$
"$2\Rightarrow3$":

It suffices to show $c\perp \lambda,\forall \lambda\in\mathscr C^\perp(X')$.
$$
\lambda'X'=0\Rightarrow (X'X)\lambda=0
$$
Let $\hat\beta$ be a solution to the normal equation.
$$
(X'X)\hat\beta=X'Y=X'X(\hat\beta+\lambda)
$$
Thus, $\hat\beta+\lambda$ is also a solution to the normal equation.

From the uniqueness, $c'(\hat\beta+\lambda)=c'\hat\beta$.

Thus, we have $c'\lambda=0$.

"$3\Rightarrow4$":

$\exists a\,{\rm s.t.}c=X'a$
$$
\mathbb E(a'Y)=a'X\beta=c'\beta
$$
"$4\Rightarrow1$":

Let $a'Y$ be unbiased for $c'\beta$.

We then have
$$
a'X\beta=c'\beta,\forall\beta
$$
Thus, $a'X=c'$.

And we have $c\in\mathscr C(X')=\mathscr C(X'X)$, i.e. $\exists \lambda,{\rm s.t.}c=(X'X)\lambda$.

Then, 
$$
\mathbb Ec'\hat\beta=\mathbb E\lambda'(X'X)(X'X)^-X'Y=\lambda'(X'X)(X'X)^-(X'X)\beta=c'\beta
$$


Remark for $3$:

Clearly $\mathscr C(X'X)\subset\mathscr C(X')$.

Also we have $rank(X'X)=rank(X')$.

Thus $\mathscr C(X'X)=\mathscr C (X')$