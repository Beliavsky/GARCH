### GARCH estimation results for some ETFs

Looking at the table of likelihoods at the end, the GJR-GARCH model, which incorporates volatility asymmetry,
has a higher likelihood for the the stock ETFs (SPY, EFA, EEM) but not the government bond ETF (TLT).

```
                  program: xgeneral_garch_data_gfort.exe
                   infile: general_garch_data.dat
              prices file: spy_efa_eem_tlt.csv
           return scaling:           100.000000
              max_iter_nm:                 1000
           log10(lik_tol):           -12.000000
                   alpha0:             0.050000
                    beta0:             0.800000
                  #assets:                    4
#rows, columns: 5490 4
first, last indices: 20030414 20250205
first, last columns: SPY TLT

             returns        SPY        EFA        EEM        TLT
mean                     0.0492     0.0363     0.0469     0.0180
sd                       1.1693     1.3104     1.7245     0.9159
skew                    -0.0872    -0.0650     0.5054     0.0797
kurt                    15.2201    14.3123    18.2294     3.4888
min                    -10.9424   -11.1633   -16.1661    -6.6682
max                     14.5198    15.8880    22.7700     7.5195

return correlations
                            SPY        EFA        EEM        TLT
SPY                      1.0000     0.8832     0.8183    -0.3127
EFA                      0.8832     1.0000     0.8679    -0.2848
EEM                      0.8183     0.8679     1.0000    -0.2662
TLT                     -0.3127    -0.2848    -0.2662     1.0000

garch log likelihoods
    symbol best_dist    normal   laplace        t4        t5        t6       t10       t20       t30
       SPY        t6   -7179.5   -7082.2   -7060.5   -7045.1   -7042.4   -7056.5   -7093.6   -7114.2
       EFA        t6   -7947.0   -7950.0   -7879.5   -7857.1   -7849.5   -7852.4   -7878.8   -7895.0
       EEM       t10   -9433.2   -9515.3   -9428.9   -9400.2   -9386.5   -9374.4   -9385.3   -9395.0
       TLT       t20   -6772.9   -6958.8   -6839.1   -6801.6   -6777.5   -6751.2   -6748.3   -6752.3

gjr_garch log likelihoods
    symbol best_dist    normal   laplace        t4        t5        t6       t10       t20       t30
       SPY        t6   -7084.4   -7012.7   -6970.9   -6954.9   -6950.1   -6960.5   -6994.3   -7014.7
       EFA       t10   -7898.3   -7919.4   -7839.2   -7816.8   -7807.2   -7806.1   -7829.4   -7844.9
       EEM       t10   -9393.5   -9498.7   -9399.2   -9370.1   -9355.3   -9340.5   -9348.8   -9357.4
       TLT       t10   -6774.2   -6956.0   -6839.0   -6798.1   -6779.1   -6751.3   -6752.6   -6753.2

                                    task    cpu_time   wall_time
   read prices and computed return stats    0.078125    0.156000
                     fit 64 GARCH models    9.734375   11.922000
                                   TOTAL    9.812500   12.078000

(10) finished xgeneral_garch_data
```
