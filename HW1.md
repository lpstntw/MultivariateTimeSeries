HW1
================
Po-Sheng Lee
2019/4/14

``` r
library(quantmod)
library(MTS)
library(tidyverse)
library(magrittr)
library(lubridate)
source("Corner.R")
```

Question 2
==========

First, specify an ARMA model for Xt using function `arima`, specifying order at (3, 0, 0), which means three lags are included and no deadtime. Then, use `tsdiag` to diagnose the model fit. Ljung-Box statistics seems quite good in this model.

In another model not reported here, I build an ARMA model at the order of (2, 0, 0) and the Ljung-Box seems fairly good enough. However, given the function `tsdiag` do not use the right degree of freedom, the p-value of Ljung-Box is inflated. After taking this into account, it is more appropriate to build an ARMA model of xt at the order of (3, 0, 0).

The results are shown as follows.

``` r
hw1 <- read.table("hw1-2.txt", skip = 1, col.names = c("yt", "xt"))

m1 <- arima(x = hw1$xt, order = c(3, 0, 0), include.mean = FALSE)

tsdiag(m1, gof.lag = 20)
```

![](HW1_files/figure-markdown_github/question%202-1-1.png)

``` r
m1
```

    ## 
    ## Call:
    ## arima(x = hw1$xt, order = c(3, 0, 0), include.mean = FALSE)
    ## 
    ## Coefficients:
    ##           ar1     ar2      ar3
    ##       -0.2454  0.1512  -0.0798
    ## s.e.   0.0446  0.0454   0.0446
    ## 
    ## sigma^2 estimated as 1.023:  log likelihood = -715.19,  aic = 1438.39

For the transfer function, we use the pre-whitening method to obtain a preliminary estimation of impulse response function. Then, apply Corner method to obtain the best guess of (r,s,b) of transfer function. It is shown as follows and the vector of (r,s,b) derived from Corner method is (1,1,1)

``` r
xf <- m1$residuals
f1 <- c(1, -m1$coef)
yf <- stats::filter(hw1$yt, f1, method = "convolution", sides = 1)
zt <- cbind(xf[4:500], yf[4:500])
colnames(zt) <- c("xf", "yf")
MTSplot(zt)
```

![](HW1_files/figure-markdown_github/question%202-2-1.png)

``` r
ccm(zt)
```

    ## [1] "Covariance matrix:"
    ##       xf     yf
    ## xf 1.028  0.274
    ## yf 0.274 73.230
    ## CCM at lag:  0 
    ##        [,1]   [,2]
    ## [1,] 1.0000 0.0315
    ## [2,] 0.0315 1.0000
    ## Simplified matrix: 
    ## CCM at lag:  1 
    ## . . 
    ## + + 
    ## CCM at lag:  2 
    ## . . 
    ## + + 
    ## CCM at lag:  3 
    ## . . 
    ## + + 
    ## CCM at lag:  4 
    ## . . 
    ## + + 
    ## CCM at lag:  5 
    ## . . 
    ## + + 
    ## CCM at lag:  6 
    ## . . 
    ## + + 
    ## CCM at lag:  7 
    ## . . 
    ## + + 
    ## CCM at lag:  8 
    ## . . 
    ## + + 
    ## CCM at lag:  9 
    ## . . 
    ## + + 
    ## CCM at lag:  10 
    ## . . 
    ## + + 
    ## CCM at lag:  11 
    ## . . 
    ## + + 
    ## CCM at lag:  12 
    ## . . 
    ## . .

![](HW1_files/figure-markdown_github/question%202-2-2.png)

    ## Hit Enter for p-value plot of individual ccm:

![](HW1_files/figure-markdown_github/question%202-2-3.png)

``` r
Corner(zt[,2], zt[,1])
```

    ## Corner Table:  
    ##       r->     1      2     3      4     5     6     7
    ##  [1,]   0 0.033  0.001 0.000  0.000 0.000 0.000 0.000
    ##  [2,]   1 0.692  0.445 0.308  0.213 0.147 0.102 0.070
    ##  [3,]   2 1.000  0.441 0.429  0.429 0.429 0.429 0.429
    ##  [4,]   3 0.808 -0.001 0.027  0.021 0.017 0.014 0.011
    ##  [5,]   4 0.653 -0.024 0.030  0.003 0.002 0.001 0.001
    ##  [6,]   5 0.558  0.034 0.014  0.001 0.001 0.001 0.000
    ##  [7,]   6 0.424 -0.016 0.009 -0.003 0.002 0.000 0.000
    ##  [8,]   7 0.351 -0.007 0.003 -0.001 0.000 0.000 0.000
    ##  [9,]   8 0.306 -0.007 0.001  0.000 0.000 0.000 0.000
    ## [10,]   9 0.288 -0.010 0.000  0.000 0.000 0.000 0.000
    ## [11,]  10 0.304  0.028 0.005  0.001 0.000 0.000 0.000
    ## 
    ## Simplified Table: 2/sqrt(T):  
    ##       r->  1   2   3   4   5   6   7  
    ##  [1,] "0"  "O" "O" "O" "O" "O" "O" "O"
    ##  [2,] "1"  "X" "X" "X" "X" "X" "X" "O"
    ##  [3,] "2"  "X" "X" "X" "X" "X" "X" "X"
    ##  [4,] "3"  "X" "O" "O" "O" "O" "O" "O"
    ##  [5,] "4"  "X" "O" "O" "O" "O" "O" "O"
    ##  [6,] "5"  "X" "O" "O" "O" "O" "O" "O"
    ##  [7,] "6"  "X" "O" "O" "O" "O" "O" "O"
    ##  [8,] "7"  "X" "O" "O" "O" "O" "O" "O"
    ##  [9,] "8"  "X" "O" "O" "O" "O" "O" "O"
    ## [10,] "9"  "X" "O" "O" "O" "O" "O" "O"
    ## [11,] "10" "X" "O" "O" "O" "O" "O" "O"

Then, use the `tfm` function to determine the proper ARMA model for Nt. The partial autocorrelation of Nt shows that the proper lag should be 4. This is corrobrated by autocorrelation of residual. The final model will be $Y\_{t} = 5.03 + \\frac{3+2.02B}{1-0.8B}BX\_{t} + \\frac{a\_{t}}{1-0.19B-0.52B^2+0.03B^3-0.002B^4}$

``` r
m2 <- tfm(hw1$yt, hw1$xt, b = 1, s = 1, p = 1)
```

    ## ARMA coefficients & s.e.: 
    ##              ar1
    ## coef.arma 0.9112
    ## se.arma   0.0181
    ## Transfer function coefficients & s.e.: 
    ##      intercept      X       
    ## v         4.47 2.0262 1.8328
    ## se.v      1.23 0.0949 0.0951

``` r
pacf(m2$nt)
```

![](HW1_files/figure-markdown_github/question%202-3-1.png)

``` r
m2_1 <- tfm(hw1$yt, hw1$xt, b = 1, s = 1, p = 4)
```

    ## ARMA coefficients & s.e.: 
    ##             ar1   ar2    ar3    ar4
    ## coef.arma 1.109 0.174 -0.729 0.3546
    ## se.arma   0.044 0.065  0.061 0.0421
    ## Transfer function coefficients & s.e.: 
    ##      intercept      X      
    ## v         4.47 1.7038 1.875
    ## se.v      1.05 0.0689 0.069

``` r
acf(m2_1$residuals)
```

![](HW1_files/figure-markdown_github/question%202-3-2.png)

``` r
m2_2 <- tfm1(hw1$yt, hw1$xt, orderN = c(4,0,0), orderX = c(1,1,1))
```

    ## Delay:  1 
    ## Transfer function coefficients & s.e.: 
    ## in the order: constant, omega, and delta: 1 2 1 
    ##       [,1]   [,2]   [,3]    [,4]
    ## v    5.033 3.0087 2.0170 0.80157
    ## se.v 0.131 0.0415 0.0463 0.00411
    ## ARMA order: 
    ## [1] 4 0 0
    ## ARMA coefficients & s.e.: 
    ##             [,1]   [,2]    [,3]    [,4]
    ## coef.arma 0.1857 0.5225 -0.0323 0.00162
    ## se.arma   0.0449 0.0459  0.0458 0.04525

Question 3
==========

Import the data from FRED using `quantmod`. The unemployment is under the tag *UNRATE* and initial claim is under *ICSA*. Data preparation is as follows

``` r
getSymbols("UNRATE", src = "FRED")
```

    ## [1] "UNRATE"

``` r
chartSeries(UNRATE)
```

![](HW1_files/figure-markdown_github/question%203%20data%20preperation-1.png)

``` r
getSymbols("ICSA", src = "FRED")
```

    ## [1] "ICSA"

``` r
chartSeries(ICSA)
```

![](HW1_files/figure-markdown_github/question%203%20data%20preperation-2.png)

``` r
## The accessed data is not structured as a data frame, tidy it for further analysis 
  
ICSA <- cbind(index(ICSA), as_data_frame(ICSA))
colnames(ICSA) <- c("date", "claims")
ICSA %<>%
  separate(date, into = c("year", "mon", "day"), sep = "-") %>%
  group_by(year, mon) %>%
  summarize(claims = mean(claims)) %>%
  ungroup()

UNRATE <- cbind(index(UNRATE), as_data_frame(UNRATE))
colnames(UNRATE) <- c("date", "unrate")
UNRATE %<>%
  separate(date, into = c("year", "mon", "day"), sep = "-")

unrate_claim_ts <- ICSA %>%
  left_join(UNRATE, by = c("year", "mon")) %>%
  select(claims, unrate) %>%
  drop_na()

MTSplot(unrate_claim_ts)
```

![](HW1_files/figure-markdown_github/question%203%20data%20preperation-3.png)

Then, employ the `VARorder` to determine the AR order. Based on the AIC, the proper order for claims should be 7 and 13 for the unemployment rate. Then use the arima to estimate and refine the model. The scalar model can be written as follows

For claims: $x\_{t}= 1.1x\_{t-1}-0.085x\_{t-2}-0.076x\_{t-5}+0.15x+{t-6}-0.12x+{t-7}+343287 $ For unemployment rate: *y*<sub>*t*</sub> = 1.04*y*<sub>*t* − 1</sub> + 0.1*y*<sub>*t* − 2</sub> − 0.072*y*<sub>*t* − 5</sub> + 0.069*y*<sub>*t* − 11</sub> − 0.19*y*<sub>*t* − 12</sub> + 0.11*y*<sub>*t* − 13</sub> + 5.92

``` r
## Use VARorder to determine the AR model
VARorder(unrate_claim_ts$claims, maxp = 20)
```

    ## selected order: aic =  7 
    ## selected order: bic =  2 
    ## selected order: hq =  2 
    ## Summary table:  
    ##        p     AIC     BIC      HQ      M(p) p-value
    ##  [1,]  0 22.6926 22.6926 22.6926    0.0000  0.0000
    ##  [2,]  1 19.6411 19.6482 19.6439 1846.5659  0.0000
    ##  [3,]  2 19.6279 19.6421 19.6334    9.9008  0.0017
    ##  [4,]  3 19.6261 19.6474 19.6344    2.9938  0.0836
    ##  [5,]  4 19.6287 19.6570 19.6397    0.3808  0.5372
    ##  [6,]  5 19.6282 19.6636 19.6420    2.1819  0.1396
    ##  [7,]  6 19.6311 19.6736 19.6476    0.1948  0.6590
    ##  [8,]  7 19.6205 19.6701 19.6398    8.2521  0.0041
    ##  [9,]  8 19.6215 19.6781 19.6435    1.3239  0.2499
    ## [10,]  9 19.6247 19.6884 19.6494    0.0009  0.9763
    ## [11,] 10 19.6277 19.6986 19.6552    0.0674  0.7952
    ## [12,] 11 19.6306 19.7085 19.6609    0.1763  0.6746
    ## [13,] 12 19.6334 19.7184 19.6664    0.2434  0.6217
    ## [14,] 13 19.6363 19.7284 19.6721    0.1445  0.7038
    ## [15,] 14 19.6395 19.7386 19.6780    0.0338  0.8542
    ## [16,] 15 19.6426 19.7489 19.6839    0.0193  0.8895
    ## [17,] 16 19.6388 19.7521 19.6828    4.1645  0.0413
    ## [18,] 17 19.6349 19.7553 19.6817    4.1664  0.0412
    ## [19,] 18 19.6374 19.7649 19.6870    0.3596  0.5487
    ## [20,] 19 19.6403 19.7749 19.6926    0.1781  0.6730
    ## [21,] 20 19.6423 19.7839 19.6973    0.7365  0.3908

``` r
VARorder(unrate_claim_ts$unrate, maxp = 20)
```

    ## selected order: aic =  13 
    ## selected order: bic =  5 
    ## selected order: hq =  7 
    ## Summary table:  
    ##        p     AIC     BIC      HQ      M(p) p-value
    ##  [1,]  0  0.9568  0.9568  0.9568    0.0000  0.0000
    ##  [2,]  1 -3.4797 -3.4727 -3.4770 2683.8064  0.0000
    ##  [3,]  2 -3.5098 -3.4957 -3.5043   20.0884  0.0000
    ##  [4,]  3 -3.5691 -3.5478 -3.5608   37.6132  0.0000
    ##  [5,]  4 -3.6057 -3.5773 -3.5946   23.9141  0.0000
    ##  [6,]  5 -3.6307 -3.5953 -3.6170   16.9730  0.0000
    ##  [7,]  6 -3.6362 -3.5937 -3.6197    5.1963  0.0226
    ##  [8,]  7 -3.6399 -3.5903 -3.6206    4.1235  0.0423
    ##  [9,]  8 -3.6367 -3.5801 -3.6147    0.0001  0.9913
    ## [10,]  9 -3.6335 -3.5698 -3.6088    0.0101  0.9200
    ## [11,] 10 -3.6304 -3.5596 -3.6029    0.0349  0.8518
    ## [12,] 11 -3.6272 -3.5493 -3.5970    0.0156  0.9008
    ## [13,] 12 -3.6298 -3.5448 -3.5968    3.4309  0.0640
    ## [14,] 13 -3.6404 -3.5483 -3.6046    8.1189  0.0044
    ## [15,] 14 -3.6373 -3.5382 -3.5988    0.0897  0.7646
    ## [16,] 15 -3.6363 -3.5301 -3.5951    1.3100  0.2524
    ## [17,] 16 -3.6343 -3.5209 -3.5902    0.6609  0.4162
    ## [18,] 17 -3.6359 -3.5155 -3.5891    2.8181  0.0932
    ## [19,] 18 -3.6332 -3.5057 -3.5836    0.2849  0.5935
    ## [20,] 19 -3.6306 -3.4960 -3.5783    0.3544  0.5516
    ## [21,] 20 -3.6310 -3.4894 -3.5760    2.1172  0.1457

``` r
m3 <- arima(unrate_claim_ts$claims, order = c(7,0,0))
m3
```

    ## 
    ## Call:
    ## arima(x = unrate_claim_ts$claims, order = c(7, 0, 0))
    ## 
    ## Coefficients:
    ##          ar1      ar2      ar3     ar4      ar5     ar6      ar7
    ##       1.0983  -0.0728  -0.0350  0.0352  -0.0888  0.1486  -0.1151
    ## s.e.  0.0397   0.0589   0.0588  0.0588   0.0589  0.0589   0.0398
    ##       intercept
    ##       343268.28
    ## s.e.   23456.32
    ## 
    ## sigma^2 estimated as 320675904:  log likelihood = -7031.54,  aic = 14081.08

``` r
m3_r <- arima(unrate_claim_ts$claims, order = c(7,0,0), fixed = c(NA, NA, 0, 0, NA, NA, NA, NA))
m3_r
```

    ## 
    ## Call:
    ## arima(x = unrate_claim_ts$claims, order = c(7, 0, 0), fixed = c(NA, NA, 0, 0, 
    ##     NA, NA, NA, NA))
    ## 
    ## Coefficients:
    ##          ar1      ar2  ar3  ar4      ar5     ar6      ar7  intercept
    ##       1.0982  -0.0848    0    0  -0.0762  0.1487  -0.1157  343287.93
    ## s.e.  0.0396   0.0464    0    0   0.0465  0.0588   0.0398   23506.47
    ## 
    ## sigma^2 estimated as 320916496:  log likelihood = -7031.77,  aic = 14077.54

``` r
tsdiag(m3_r, gof.lag = 24)
```

![](HW1_files/figure-markdown_github/question%203-1-1.png)

``` r
m3_1 <- arima(unrate_claim_ts$unrate, order = c(13,0,0))
m3_1
```

    ## 
    ## Call:
    ## arima(x = unrate_claim_ts$unrate, order = c(13, 0, 0))
    ## 
    ## Coefficients:
    ##          ar1     ar2     ar3     ar4      ar5      ar6      ar7    ar8
    ##       1.0360  0.0967  0.0048  0.0083  -0.0765  -0.0068  -0.0795  0.019
    ## s.e.  0.0396  0.0568  0.0569  0.0569   0.0569   0.0571   0.0570  0.057
    ##          ar9     ar10    ar11     ar12    ar13  intercept
    ##       0.0009  -0.0118  0.0726  -0.1915  0.1154     5.9169
    ## s.e.  0.0569   0.0570  0.0571   0.0569  0.0398     0.4907
    ## 
    ## sigma^2 estimated as 0.0249:  log likelihood = 265.35,  aic = -500.7

``` r
m3_1_r <- arima(unrate_claim_ts$unrate, order = c(13,0,0), fixed = c(NA, NA, 0, 0, NA, 0, NA, 0, 0, 0, NA, NA, NA, NA))
m3_1_r
```

    ## 
    ## Call:
    ## arima(x = unrate_claim_ts$unrate, order = c(13, 0, 0), fixed = c(NA, NA, 0, 
    ##     0, NA, 0, NA, 0, 0, 0, NA, NA, NA, NA))
    ## 
    ## Coefficients:
    ##          ar1     ar2  ar3  ar4      ar5  ar6      ar7  ar8  ar9  ar10
    ##       1.0359  0.1023    0    0  -0.0719    0  -0.0715    0    0     0
    ## s.e.  0.0393  0.0483    0    0   0.0401    0   0.0383    0    0     0
    ##         ar11     ar12    ar13  intercept
    ##       0.0691  -0.1911  0.1147     5.9178
    ## s.e.  0.0473   0.0567  0.0397     0.4915
    ## 
    ## sigma^2 estimated as 0.02491:  log likelihood = 265.24,  aic = -512.48

``` r
tsdiag(m3_1_r, gof.lag = 24)
```

![](HW1_files/figure-markdown_github/question%203-1-2.png)

The Granger Causality can be tested by `GrangerTest` in `MTS`. Using AR(13) model, we find that both direction has a high p-value, thus resulting in no unidirectional causality.

``` r
VARorder(unrate_claim_ts)
```

    ## selected order: aic =  12 
    ## selected order: bic =  2 
    ## selected order: hq =  7 
    ## Summary table:  
    ##        p     AIC     BIC      HQ      M(p) p-value
    ##  [1,]  0 22.7177 22.7177 22.7177    0.0000  0.0000
    ##  [2,]  1 15.8184 15.8467 15.8294 4219.8218  0.0000
    ##  [3,]  2 15.7622 15.8188 15.7842   41.9452  0.0000
    ##  [4,]  3 15.7371 15.8220 15.7701   22.9820  0.0001
    ##  [5,]  4 15.7244 15.8377 15.7684   15.3690  0.0040
    ##  [6,]  5 15.7113 15.8530 15.7664   15.5447  0.0037
    ##  [7,]  6 15.6881 15.8581 15.7542   21.6093  0.0002
    ##  [8,]  7 15.6611 15.8594 15.7382   23.7942  0.0001
    ##  [9,]  8 15.6641 15.8908 15.7522    5.8081  0.2139
    ## [10,]  9 15.6564 15.9114 15.7555   12.1670  0.0162
    ## [11,] 10 15.6562 15.9395 15.7663    7.6837  0.1039
    ## [12,] 11 15.6654 15.9771 15.7865    2.1015  0.7171
    ## [13,] 12 15.6431 15.9830 15.7751   20.6727  0.0004
    ## [14,] 13 15.6477 16.0160 15.7908    4.7382  0.3152

``` r
GrangerTest(unrate_claim_ts, p = 12, locInput = c(1))
```

    ## Number of targeted zero parameters:  12 
    ## Chi-square test for Granger Causality and p-value:  30.64864 0.002228052

``` r
GrangerTest(unrate_claim_ts, p = 12, locInput = c(2))
```

    ## Number of targeted zero parameters:  12 
    ## Chi-square test for Granger Causality and p-value:  173.8459 0

Question 4
==========

Get the time series data with `quantmod`. The price report time of crude oil is three days ahead (or four days lag) of the regular gas. I take the crude oil price as being ahead of the regular oil; therefore I add three days to the pricing day of crude oil.

``` r
getSymbols("WCOILWTICO", src = "FRED")
```

    ## [1] "WCOILWTICO"

``` r
chartSeries(WCOILWTICO)
```

![](HW1_files/figure-markdown_github/question%204%20data%20preparation-1.png)

``` r
getSymbols("GASREGCOVW", src = "FRED")
```

    ## [1] "GASREGCOVW"

``` r
chartSeries(GASREGCOVW)
```

![](HW1_files/figure-markdown_github/question%204%20data%20preparation-2.png)

``` r
crude_oil <- cbind(index(WCOILWTICO), as_data_frame(WCOILWTICO))
colnames(crude_oil) <- c("date", "crude")
crude_oil %<>%
  mutate(date = date + 3)

reg_gas <- cbind(index(GASREGCOVW), as_data_frame(GASREGCOVW))
colnames(reg_gas) <- c("date", "reg")

i <- ymd("1996-01-01") %--% ymd("2018-12-31")

oil_ts <- crude_oil %>%
  left_join(reg_gas) %>%
  drop_na() %>%
  mutate(date = as_date(date)) %>%
  filter(date %within% i) %>%
  mutate(crude = log10(crude), reg = log10(reg)) %>%
  select(crude, reg)

MTSplot(oil_ts)
```

![](HW1_files/figure-markdown_github/question%204%20data%20preparation-3.png)

Use `VARorder` and `GrangerTest` to see whether there is Granger Caulsality among these tow variables. Use AR(9) to test the Granger Causality. No unidirectional causality is detected.

``` r
VARorder(oil_ts)
```

    ## selected order: aic =  9 
    ## selected order: bic =  3 
    ## selected order: hq =  3 
    ## Summary table:  
    ##        p      AIC      BIC       HQ      M(p) p-value
    ##  [1,]  0  -9.5203  -9.5203  -9.5203    0.0000  0.0000
    ##  [2,]  1 -17.6585 -17.6416 -17.6521 9647.6059  0.0000
    ##  [3,]  2 -18.0279 -17.9940 -18.0151  444.7043  0.0000
    ##  [4,]  3 -18.0520 -18.0012 -18.0329   36.3242  0.0000
    ##  [5,]  4 -18.0566 -17.9888 -18.0310   13.2245  0.0102
    ##  [6,]  5 -18.0554 -17.9706 -18.0235    6.4191  0.1700
    ##  [7,]  6 -18.0567 -17.9550 -18.0184    9.3773  0.0523
    ##  [8,]  7 -18.0547 -17.9360 -18.0100    5.4922  0.2404
    ##  [9,]  8 -18.0502 -17.9146 -17.9991    2.5126  0.6424
    ## [10,]  9 -18.0599 -17.9073 -18.0024   19.0610  0.0008
    ## [11,] 10 -18.0598 -17.8903 -17.9959    7.6971  0.1033
    ## [12,] 11 -18.0538 -17.8673 -17.9836    0.7786  0.9413
    ## [13,] 12 -18.0552 -17.8517 -17.9786    9.3568  0.0528
    ## [14,] 13 -18.0491 -17.8287 -17.9661    0.6507  0.9573

``` r
GrangerTest(oil_ts, p = 9, locInput = c(1))
```

    ## Number of targeted zero parameters:  9 
    ## Chi-square test for Granger Causality and p-value:  37.90163 1.815407e-05

``` r
GrangerTest(oil_ts, p = 9, locInput = c(2))
```

    ## Number of targeted zero parameters:  9 
    ## Chi-square test for Granger Causality and p-value:  60.44978 1.098027e-09

Question 5
==========

The plot and ccm is shown below. Multivariate Ljung-Box statistic rejects the null hypothesis that *ρ*<sub>*i*</sub> = 0

``` r
c <- matrix(c(.9, -.3, .4, .6), 2, 2)
s <- matrix(c(1.5, .3, .3, 1), 2, 2)

m5 <- VARMAsim(200, arlags=c(1) , phi= c , sigma= s )
zt5 <- m5$series

MTSplot(zt5)
```

![](HW1_files/figure-markdown_github/VARMAsim-1.png)

``` r
zt5_ccm <- ccm(zt5, lag = 3)
```

    ## [1] "Covariance matrix:"
    ##       [,1]  [,2]
    ## [1,]  6.43 -1.09
    ## [2,] -1.09  3.39
    ## CCM at lag:  0 
    ##        [,1]   [,2]
    ## [1,]  1.000 -0.233
    ## [2,] -0.233  1.000
    ## Simplified matrix: 
    ## CCM at lag:  1 
    ## + . 
    ## - + 
    ## CCM at lag:  2 
    ## + + 
    ## - + 
    ## CCM at lag:  3 
    ## + + 
    ## - .

![](HW1_files/figure-markdown_github/VARMAsim-2.png)

    ## Hit Enter for p-value plot of individual ccm:

![](HW1_files/figure-markdown_github/VARMAsim-3.png)

``` r
zt5_ccm$ccm
```

    ##            [,1]       [,2]       [,3]       [,4]
    ## [1,]  1.0000000  0.8323839  0.5756921  0.2862727
    ## [2,] -0.2329955 -0.5558873 -0.6981885 -0.7090978
    ## [3,] -0.2329955  0.1009845  0.3434335  0.4296676
    ## [4,]  1.0000000  0.7333385  0.4505815  0.1405265

``` r
mq(zt5, lag = 5)
```

    ## Ljung-Box Statistics:  
    ##        m       Q(m)     df    p-value
    ## [1,]     1       290       4        0
    ## [2,]     2       512       8        0
    ## [3,]     3       674      12        0
    ## [4,]     4       789      16        0
    ## [5,]     5       868      20        0

![](HW1_files/figure-markdown_github/VARMAsim-4.png)
