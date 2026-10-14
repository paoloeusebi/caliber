# Category probabilities and expected score of one GRM item

\`grm_probs()\` gives \\P(C = 0..3 \mid \eta)\\ for one item with
\\\mathrm{logit}\\ P(C \ge k \mid \eta) = a \eta + d_k\\;
\`grm_expected()\` gives \\E\[C \mid \eta\]\\.

## Usage

``` r
grm_probs(eta, a, d)

grm_expected(eta, a, d)
```

## Arguments

- eta:

  latent values

- a:

  slope

- d:

  intercepts \\d_1 \> d_2 \> d_3\\

## Value

\`grm_probs()\`: length(eta) x 4 matrix; \`grm_expected()\`: vector
