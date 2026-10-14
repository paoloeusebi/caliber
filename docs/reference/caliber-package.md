# caliber: monotone calibration of a digital measure to a multi-item PRO

One latent trait \\\eta \sim N(0, 1)\\ drives an ordinal multi-item PRO
(logistic graded response model, items scored 0 to 3) and a continuous
digital or performance measure \\y\\: \$\$\mathrm{logit}\\ P(C_j \ge k
\mid \eta) = a_j \eta + d\_{jk}, \qquad y = \mu(\eta) + e, \quad e \sim
N(0, \sigma(\eta)^2), \quad \log \sigma(\eta) = g_0 + g_1 \eta,\$\$ with
\\\mu(\eta) = \sum_m \beta_m B_m(\eta)\\ a non-increasing cubic B-spline
(24 coefficients, 20 interior knots on \\\[-4, 4\]\\). Higher \\\eta\\
means more of the PRO construct (worse sleep, more mobility difficulty);
the digital measure falls as \\\eta\\ rises.

## Details

The curve is fitted by monotonic likelihood boosting of B-splines
(Leitenstorfer and Tutz 2007) with the AIC stop, in one of three ways
that differ only in how \\\eta\\ is handled (see \[fit_ms()\]):

- MS-Joint:

  concurrent calibration: GRM items and curve estimated together,
  \\\eta\\ integrated out over \\p(\eta \mid C, y)\\;

- MS-EAP:

  EAP plug-in: GRM on the PRO alone, EAP scores treated as known;

- MS-Fixed:

  fixed-item calibration: GRM on the PRO alone, items then fixed,
  \\\eta\\ integrated out (the recommended primary analysis).

Throughout the code the digital measure is called \`sqi\` (the
sleep-quality index of the worked example); any continuous measure that
falls with the latent trait can be used (gait speed in the NHANES
example).

## Author

**Maintainer**: Paolo Eusebi <paoloeusebi@gmail.com>

Authors:

- Paolo Eusebi <paoloeusebi@gmail.com>
