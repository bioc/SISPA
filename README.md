<h4>SISPA: Method for Sample Integrated Set Profile Analysis</h4>
Bhakti Dwivedi, Ph.D and Jeanne Kowalski, Ph.D; 
Winship Cancer Institute, Emory University, Atlanta, 30322, USA

###Introduction
Sample Integrated Gene Set Analysis (SISPA) is a method designed to define sample groups with similar gene set enrichment profiles. The user specifies a gene list of interest and sample by gene molecular data (expression, methylation, variant, or copy change data) to obtain gene set enrichment scores by each sample. The score statistics is rank ordered by the desired profile (e.g., upregulated or downregulated) for samples. A change point model is then applied to the sample scores to identify groups of samples that show similar gene set profile patterns. Samples are ranked by desired profile activity score and grouped by samples with and without profile activity. Figure 1 shows the schematic representation of the SISPA method overview.

###Installation
To install the latest build directly from GitHub, run the following command from R console:

```{r}
if (!require("devtools"))
  install.packages("devtools")
  devtools::install_github("BhaktiDwivedi/SISPA")
```
###Getting Started
For an introduction and example, please see ./vignettes/SISPA.html

