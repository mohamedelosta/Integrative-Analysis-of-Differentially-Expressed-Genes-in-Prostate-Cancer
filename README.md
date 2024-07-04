# Integrative Analysis of Differentialy Expressed genes in Prostate Cancer

## Abstract:

Prostate cancer is a complex and heterogeneous disease with diverse molecular subtypes and clinical outcomes. In this thesis, we conducted a comprehensive analysis of prostate cancer using integrative omics approaches to unravel key biomarkers and pathways associated with disease progression.
The study utilized RNA-seq data from 553 samples and machine learning algorithms to identify 1224 differentially expressed genes (DEGs) between normal and tumor samples. Gene set enrichment analysis (GSEA) revealed significant associations with biological processes such as serine-type endopeptidase activity regulation and pathways linked to cancer growth and metastasis.
Machine learning models demonstrated high accuracy in predicting sample classification, with top-ranked genes including EPHA10, ABCG2, and GSTM3 identified as potential biomarkers for prostate cancer. Network analysis highlighted gene clusters and interactions, further validating our findings through publication enrichment and biological process enrichment.
Our thesis proposes strategies for experimental validation, multi-omics integration, and clinical translation to advance precision medicine in prostate cancer diagnosis and treatment. The findings contribute to a deeper understanding of prostate cancer biology and hold promise for personalized therapeutic interventions and improved patient outcomes.

# THE AIM

Unlike the current methods, the thesis aims to help detect new genetic biomarkers of Prostate Cancer. For this purpose, I used multiple approaches beginning with differential expression analysis to focus on the statistically significant different genes. those genes were used in advanced steps to gain biological insights by enrichment analysis and network analysis and to prioritize the most important genes some classification models were trained to know which are the most important genes that makes the cell cancerous.

But this READme aims to show exactly the steps I followed to conduct my analysis regarding Data collection, filtration exploration, and analysis including differential expression analysis and machine learning. in addition some comments on the methods I used.

# Data Collection

The gene expression data utilized in this study was sourced from the GDC Data Portal and imported directly into our development environment using the TCGAbiolinks R\
package. This method was way easier and faster to get your data with the desired features but it maybe computationaly consuming to do it in case of big samples.

```{r eval=FALSE}
##building query data

gdc.query = GDCquery(
  project = "TCGA-PRAD" ,
  data.category = 'Transcriptome Profiling' ,
  experimental.strategy = 'RNA-Seq',
  workflow.type = 'STAR - Counts',
  access = 'open')

results<- getResults(gdc.query)
```

Then we download the data and retrieve the count matrix and the metadata.

```{r eval=FALSE}
#downlaod the data 

GDCdownload(gdc.query)

#prepare the data

TCGA.data <- GDCprepare(gdc.query, summarizedExperiment = TRUE)

data.matrix <- assay(TCGA.data, 'unstranded')

rowData <- rowData(TCGA.data)

colData <-  colData(TCGA.data)
```

This data originated from the TCGA-Prostate Adenocarcinoma (PRAD) project and had already undergone essential preprocessing steps, converting raw sequencing data from Fastq files into raw count expression data. We opted for the RNA-Seq experimental strategy and chose the STAR-Counts workflow type, specifically focusing on "unstranded" read count analysis that considers both sense and antisense transcripts. Access to this data was granted through an open-access license, ensuring its availability for research purposes. The dataset comprised 553 cases, representing samples from individuals diagnosed with prostate adenocarcinoma (PRAD). Importantly, the data were obtained in their raw form without normalization or transformation, allowing for a thorough analysis of the gene expression values and metadata associated with the samples.

## Filtration and Normalization

\
Before normalization, RNAs with near-zero or zero read count values across all samples
were excluded from the study. Additionally, a filtration step based on protein-coding
information was applied to focus on actively transcribed genes relevant to the functional
aspects of the transcriptome. This filtration ensured the analysis centered on biologically
relevant genes with known protein-coding potential. Following the filtration process for
protein-coding genes, a refined Gene Expression Matrix comprising 18,611 protein-coding genes was constructed. This filtered dataset served as the basis for our exploratory analysis. Subsequently, the data was normalized using the variance-stabilizing transform (VST) through the DESeq2 R-project package. VST was employed to accommodate the different variances between RNAs, as it diminishes the dependency of variance on the mean value of the read counts. This normalization strategy aimed to reduce systematic errors and enhance the reliability of our subsequent exploratory analysis.

## Exploratory Data Analysis

Following the essential data preprocessing steps, the dataset obtained from TCGA was explored through visualization techniques, including Principal Component Analysis (PCA), density plots, and box plots. The subsequentt section elaborates on the outcomes and insights from each visualization method.

The study's initial phase, known as exploratory analysis, aims to run the data through several visualization and filtration processes for three reasons: to find RNAs that might not be relevant to our investigation; to find samples that might be incorrectly labeled, misplaced, or otherwise flawed and need to be excluded from the analysis; and, finally, to confirm that the available tumor and non-tumor samples show a distinct and obvious difference. These goals will mostly be achieved through data visualization, making the patterns displayed by the samples and RNAs easy to see.

### Density Plot

In this section, the log-transformed counts were used to produce density plots. Upon examining the density plots for both tumor and normal samples, a consistent distribution of expression values was observed across all samples in both categories. This uniformity suggests the absence of outliers or abnormalities within the sample set. This finding underscores the reliability and consistency of the data set, laying a strong foundation for further analysis and interpretation.

![Tumor samples density plots](images/tumor%20density.png){fig-align="left" width="1000"}

![Normal samples density plot](images/normal%20density.png)

## Box Plots

The utilization of box plots to visualize raw expression values showcased high comparability among the samples before normalization. Notably, these plots did not exhibit any peculiar or abnormal distributions within any of the samples. This uniformity across the dataset, even before normalization, signifies the robustness of the data and bolsters confidence in the subsequent analytical procedures.


## PCA

Principal Component Analysis (PCA) was employed on the DESeqTransform object, which contains normalized data processed using the DESeq2 package's VTS method. PCA revealed distinct clustering between tumor and normal samples with some overlap. This overlap may be attributed to prostate cancer's inherent heterogeneity and the intricate transcriptional patterns associated with it.

![PCA of normalized counts](images/PCAnormalized.png)

# Differential Expression Analysis

Differential analysis was conducted using the DESeq2 package with stringent criteria\
(\|log2 FoldChange\| \> 1.5 and padj \< 0.01). to identify genes with significant expression\
changes in prostate adenocarcinoma (PRAD). A total of 1224 differentially expressed\
genes (DEGs) were identified, with 675 upregulated and 549 downregulated in tumor\
samples compared to normal samples. A volcano plot was constructed to visualize these\
changes and identify the top genes based on adjusted P and log2 fold change values\
for further investigation. A volcano plot was constructed to visualize these changes\
and identify the top genes based on adjusted P and log2 fold change values for further\
investigation.

![Volcano Plot reveals the upregulated and downregulated genes](images/f.jpg)

## Heat map

Subsequently, a heatmap was generated to depict gene expression levels across samples, focusing on the top 100 DEGs. The heatmap highlighted the heterogeneity of\
tumor samples compared to the more consistent transcriptional profiles of normal sam-\
ples within this gene subset. Notably, the heatmap illustrated the downregulation of\
the first 13 DEGs in normal samples, while the remaining 87 genes were predominantly\
upregulated in most normal samples.

![Heatmap of Top 100 DEGs](images/Picture1.png)\

## PCA for Differentially Expressed Genes

\
To delve deeper into the transcriptomic profiles and visualize patterns and variations\
within the dataset, we conducted principal component analysis (PCA) on the expression\
data. We generated both 2D and 3D PCA plots. The PCA analysis yielded clear segregation between tumor and normal samples in the 2D and 3D plots. This observation underscores the robustness of our findings in discriminating between normal and cancerous tissues in the context of prostate adenocarcinoma\
(PRAD). Notably, this enhanced separation in the PCA plots with selected genes suggests improved discriminative power compared to PCA using all genes, further validating\
the efficacy of our analysis.

![PCA of DEGs](images/PCA2D.png)

```{r echo=FALSE}
# Load necessary libraries
library(magrittr)
library(plotly)
# Create a 3D plotly plot
pca3d <- plot_ly(
  data = pca_scores,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~class,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 5),
  text = ~paste("PC1:", round(PC1, 2), "<br>PC2:", round(PC2, 2), "<br>PC3:", round(PC3, 2))
) %>% layout(
    title = "3D PCA Plot",
    scene = list(
      xaxis = list(title = 'PC 1'),
      yaxis = list(title = 'PC 2'),
      zaxis = list(title = 'PC 3')
    )
  )

# Display the plot
pca3d

```

## Machine Learning

The machine learning part of the thesis aims to prioritize the most important genes differentiating normal and tumor tissues by training accurate classification models to shed light on the most important variables in the data.

![machine learning workflow](images/ml.jpg)

### 1. Data Pre-processing

The process begins with acquiring gene expression data from cancer and normal samples. This raw data undergoes pre-processing steps such as normalization, imputation of missing values, and transformation to ensure consistency and reliability. Our data set is the DEGs expression values, which are already preprocessed.

### 2. Feature selection

In the realm of biomarker discovery, identifying the most relevant features (genes) is crucial to reduce dimensionality and retain only the informative features. This step is done by obtaining DEGs between samples.

### 3. Data Splitting

The dataset is divided into training 70%, testing 20%, and validation 10% sets. The training set is used to train the machine learning models, the testing set evaluates model performance during development, and the validation set assesses the model's generalization on unseen data.

### 4. Model Selection

we chose 4 clustering algorithms for supervised learning in order to differentiate between samples depending on their transcription profiles of DEGs.

#### SVM (Support Vector Machines):

SVM is typically used for classification tasks that aim to find a hyperplane that best separates different classes. In clustering, SVM can be employed in a way known as one-class SVM, which essentially tries to isolate one class of data from the rest. This can be useful in outlier detection or anomaly detection within clusters.

#### RF (Random Forest):

RF is an ensemble learning method that constructs multiple decision trees during training and outputs the mode of the classes (for classification) or the average prediction (for regression). In clustering, RF can be used indirectly by training it on the data and then using the leaf node assignments of the trees as cluster labels. This approach is known as "using RF for clustering via leaf node assignments."

#### XGBoost:

XGBoost is another ensemble learning technique that is particularly powerful for classification and regression tasks. Like RF, XGBoost can be indirectly used for clustering by training it on the data and then extracting cluster assignments based on the leaf nodes of the boosted trees.

#### Neural Networks:

For the sake of the experiment neural networks were used to classify the samples. the architecture consisted of 2 hidden layers 200 nodes each in addition to 1224 input (number of DEGs) and 2 nodes output layers. ReLU activation function was used.

#### There is no doubt that NNs are not the best tool in our case due to the small size of the data and the lack of interpretability in DL models which is the aim of using AI for our data. 

### 5. Model Training

Various machine learning algorithms, such as Support Vector Machines (SVM), Random Forests, extreme Gradient Boosting (XGB), and Neural Networks, are trained using the training data set. These models learn patterns and relationships within the data to accurately classify cancer and normal samples.

### 6. Model Evaluation

The trained models are evaluated using the testing dataset to assess their performance metrics such as accuracy, precision, recall, and F1 score. This evaluation ensures the models generalize well and can reliably classify unseen samples.

### 7. Feature Importance

Post-training models provide insights into feature importance, highlighting the genes (variables) that contribute significantly to the classification task. This step is pivotal in biomarker discovery, as it identifies potential biomarkers for further validation and clinical translation.

## Model Evaluation

The trained models underwent evaluation post-validation, assessing their performance through diverse metrics based on predictions made on the test set. we estimated 4 metrics

![Models Estimated Metrics](images/metrics.png)
![metrics](https://github.com/mohamedelosta/Integrative-Analysis-of-Differentially-Expressed-Genes-in-Prostate-Cancer/assets/143076941/97a8e244-daa6-4e91-a01e-37f499bbaafb)

As depicted in the figure, all models exhibited notably high accuracy, underscoring\
their robust performance. After confirming their high predictive accuracy, we extracted\
the most influential variables, specifically the genes pivotal in sample classification. This\
iterative process of model evaluation, validation, and feature selection validates our pre-\
dictive models' efficacy and unveils key genetic determinants critical for accurate sample\
classification. The most important variables in each model have been obtained to gain\
insights about what could be the possible biomarkers.

![]()

The fact that all of our models had high accuracy encouraged us to choose to combine
all importance scores and rank top genes to obtain the best biomarker based on their
score. These models gave a lot of insights about the disease. The highest score gene
EPHA10 (erythropoietin-producing hepatocellular carcinoma receptor A10) was found
to be over-expressed in prostate cancer. Also ABCG2. GSTM3 has different polymorphisms in various tumor cells and regulates tumorigenesis, cell invasion, metastasis, chemoresistance, and oxidative stress.


