<center> <h4>SISPA: Method for Sample Integrated Set Profile Analysis</h4> </center>
========================================================
<center> <h6>Bhakti Dwivedi$^1$ and Jeanne Kowalski$^1$}</h6></center>
\begin{document}

*Sample Integrated Set Profile Analysis (SISPA) generates sample profile identifiers based on sample z-scores*


```r
#load required packages
require(GSVA)
require(changepoint)
require(data.table)
require(plyr)
require(ggplot2)
```

### Gene-set GSVA enrichment analysis

```r
#load the input data
setwd("/Users/bhaktidwivedi/Documents/R-code_GitHub/R_packages/SISPA/data")
load("SISPA_data.Rda")
```




```r
#perform GSVA with "zscore" method
gsva_results <- callGSVA(expr,genes)
```


```r
#read the gene-set enrichment scores
head(gsva_results)
```

```
##          samples      zscore
## 1 MMRF_1024_1_BM -1.29341981
## 2 MMRF_1184_1_BM -0.57416846
## 3 MMRF_1293_1_BM  0.02372068
## 4 MMRF_1423_1_BM -0.18449918
## 5 MMRF_1512_1_BM -0.47431570
## 6 MMRF_1520_1_BM -1.55127672
```



### Sample Profile Identifier Analysis




#### Samples with Increased profile activity

```r
#identification of sample profiles with 'increased' change using changepoint model
cpt_on_samples <- cptSamples(gsva_results,dir="up",cpt_data="var",cpt_method="BinSeg",cpt_max=6)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 


```r
head(cpt_on_samples[[2]])
```

```
##           samples   zscore changepoints sample_groups
## 13 MMRF_1613_1_BM 3.378475            1             1
## 10 MMRF_1584_1_BM 2.927757            1             1
## 98 MMRF_1823_1_BM 2.900674            1             1
## 46 MMRF_1683_1_BM 2.867514            1             1
## 21 MMRF_1630_1_BM 2.830496            1             1
## 97 MMRF_1822_1_BM 2.611923            1             1
```




```r
wfplot <- waterfallplot(cpt_on_samples[[2]])
plot(wfplot)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png) 


```r
count_data <- count(cpt_on_samples[[2]],"sample_groups")
plot_sample_freq <- freqplot(count_data$freq)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png) 

#### Samples with Decreased profile activity

```r
#identification of sample profiles with 'decreased' change using changepoint model
cpt_on_samples <- cptSamples(gsva_results,dir="down",cpt_data="var",cpt_method="BinSeg",cpt_max=6)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png) 


```r
head(cpt_on_samples[[2]])
```

```
##           samples    zscore changepoints sample_groups
## 71 MMRF_1758_1_BM -3.838842            1             1
## 96 MMRF_1819_1_BM -2.740161            1             1
## 8  MMRF_1543_1_BM -2.671651            1             1
## 87 MMRF_1795_1_BM -2.313509            1             1
## 57 MMRF_1716_1_BM -2.117718            1             1
## 31 MMRF_1648_1_BM -1.993919            1             1
```


```r
wfplot <- waterfallplot(cpt_on_samples[[2]])
plot(wfplot)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png) 


```r
count_data <- count(cpt_on_samples[[2]],"sample_groups")
plot_sample_freq <- freqplot(count_data$freq)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png) 
