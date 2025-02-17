## GJR-GARCH Simulation and Estimation

The program `xgjr_garch.f90` simulates a 10000 observations of a GJR-GARCH(1,1) process and then fits a GJR-GARCH model to the simulated returns. The [GJR-GARCH](https://frds.io/algorithms/gjr-garch/) model "extends the basic GARCH(1,1) by accounting for leverage effects, where bad news (negative returns) has a greater impact on volatility than good news." Other main programs fit ARCH and GARCH models to financial returns, computed from daily ETF closing prices in `spy_efa_eem_tlt.csv`.

### Code Overview

- **Simulation:**  
  The program generates 10,000 observations from a GJR-GARCH process using the true parameters.

- **Estimation:**  
  The model is then fitted to the simulated returns using the Nelder-Mead algorithm. The estimated parameters, negative log-likelihood, and other diagnostic statistics are output.

- **Diagnostics:**  
  The output includes basic statistics (mean, standard deviation, skewness, kurtosis, min, and max) for both the true conditional standard deviations (`sigma`) and the estimated standard deviations (`sigma_est`).  
  Additional diagnostics include:
  - Kurtosis of returns and standardized returns.
  - Autocorrelation function (ACF) of squared returns and of the squared standardized returns.
  - Correlation between the true and estimated volatilities.

## Explanation of Sample Output in `results.txt`

### 1. Model Specification

- **GARCH Model:**  
  The selected model is printed as **gjr_garch**.
  
- **Conditional Distribution:**  
  The returns are modeled using a **normal** distribution.

### 2. True and Estimated Parameters

The following table compares the true parameters with the estimated parameters obtained from fitting the model:

| Parameter | True Value | Estimated Value |
|-----------|------------|-----------------|
| **mu**    | 0.000000   | -0.017736       |
| **omega** | 0.100000   | 0.117329        |
| **alpha** | 0.100000   | 0.087711        |
| **gamma** | 0.100000   | 0.101340        |
| **beta**  | 0.800000   | 0.800561        |

This close agreement indicates that the estimation procedure is performing well.

### 3. Log-Likelihood

- The log-likelihood is reported as **-16677.550339**.  
- The value computed by the negative log-likelihood function matches this value, confirming the consistency of the likelihood calculation.

### 4. Volatility Statistics

- **True Volatility (`sigma`):**  
  Basic statistics (mean, standard deviation, skew, kurtosis, min, and max) are computed for the simulated conditional standard deviation. For example, the mean is approximately 1.321965 with a standard deviation of 0.363676.

- **Estimated Volatility (`sigma_est`):**  
  The corresponding statistics for the estimated volatilities are nearly identical to the true values (e.g., mean ≈ 1.319146), indicating a good model fit.

### 5. Additional Diagnostics

- **Kurtosis:**  
  The kurtosis of the raw returns is near unity, while the kurtosis of the standardized returns (both ret/sigma and ret/sigma_est) is close to zero. This suggests that the model has successfully normalized the returns.

- **Correlation:**  
  The correlation between `sigma` and `sigma_est` is extremely high (≈ 0.999808), demonstrating that the estimated volatility closely tracks the true volatility.

- **Autocorrelation Function (ACF):**  
  - The ACF of squared returns shows moderate autocorrelation at low lags, a common sign of volatility clustering.  
  - In contrast, the ACF of the squared standardized returns (ret/sigma_est)^2 is near zero, confirming that the model has effectively removed the time-dependence in the volatility.
