# Substituent Effects on Chloride Binding in Thiourea Derivatives (DFT)

Linear regression analysis of Hammett σ substituent constants vs. computed
binding energies (kcal/mol), with a publication-ready scatter plot generated
in Python (NumPy, SciPy, Matplotlib).

## Data

| Substituent | Hammett σ | Binding Energy (kcal/mol) |
| ----------- | --------: | -------------------------: |
| H           |      0.00 |                      -29.1 |
| CH₃         |     -0.17 |                      -43.2 |
| OCH₃        |     -0.27 |                      -41.1 |
| Cl          |     +0.23 |                      -35.5 |

## What the script does

- Fits an ordinary least-squares (OLS) linear regression:
  `Binding Energy = slope × σ + intercept`
- Reports slope, intercept, R², Pearson r, p-value, and 95% CI for the slope
- Prints observed vs. fitted values and residuals for each substituent
- Produces a scientific-report-style scatter plot with:
  - labelled data points (H, CH₃, OCH₃, Cl)
  - fitted trendline
  - equation and R² annotation box
  - 300 dpi PNG output suitable for reports/publications

## Requirements

```
python >= 3.9
numpy
scipy
matplotlib
```

Install with:

```bash
pip install numpy scipy matplotlib
```

## Usage

```bash
python hammett_binding_energy_analysis.py
```

This prints the regression statistics to the console and saves
`hammett_binding_energy.png` in the working directory.

## Results (current dataset, n = 4)

| Parameter | Value |
|---|---|
| Slope | 17.1133 kcal/mol per σ unit (SE = 16.44) |
| Intercept | -36.3266 kcal/mol (SE = 3.23) |
| R² | 0.3515 |
| Pearson r | 0.5929 |
| p-value | 0.407 |
| Equation | y = 17.1133x - 36.3266 |

## Notes

With only 4 data points, R² = 0.35 and p = 0.41 indicate a weak, statistically
inconclusive correlation between σ and binding energy. The H substituent
carries the largest residual and appears to drive most of the scatter;
excluding it raises R² to ~0.81 (3 points only — for diagnostic purposes,
not a valid conclusion). A broader substituent set spanning a wider σ range,
and/or dual-parameter (F, R) or σ⁺ scales, would give a more robust test of
the correlation.

