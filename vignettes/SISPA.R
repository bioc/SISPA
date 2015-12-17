## ---- message=F, warning=F-----------------------------------------------
library("SISPA")
load("SISPA_data.Rda")

## ---- message=F, warning=F, fig.width = 6, fig.height = 4.5, fig.cap="Identified changepoints on the underlying data zscores"----
sampleScores <- SISPA(feature=1,f1.df=ExpressionSet,f1.genes=GeneSet,f1.profile="up")
head(sampleScores)

## ---- message=F, warning=F, fig.width = 4.5, fig.height = 3, fig.cap="Waterfall plot of zscores for samples with (orange) and without profile activity (grey)"----
waterfallplot(sampleScores)

## ---- message=F, warning=F, fig.width = 4.5, fig.height = 3, fig.cap="Bar plot distribution of samples with (orange) and without profile activity"----
freqplot(sampleScores)

## ---- message=F, warning=F, fig.width = 6, fig.height = 4.5, fig.cap="Identified changepoints on the underlying data zscores"----
sampleScores <- SISPA(feature=2,f1.df=ExpressionSet,f1.genes=GeneSet,f1.profile="up",
                      f2.df=VariantSet,f2.genes=VariantGeneSet,f2.profile="up")
head(sampleScores)

## ---- message=F, warning=F, fig.width = 4.5, fig.height = 3, fig.cap="Waterfall plot of zscores for samples with (orange) and without profile activity (grey)"----
waterfallplot(sampleScores)

## ---- message=F, warning=F, fig.width = 4.5, fig.height = 3, fig.cap="Bar plot distribution of samples with (orange) and without profile activity"----
freqplot(sampleScores)

