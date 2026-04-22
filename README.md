# MIRACLE: Machine learning Integration for RNA Cancer-specific LncRNA Exploration 
This research was conducted as a part of the **2025-2026 undergraduate research cohort** at **Albany State University**, Albany GA, USA 31705

Conducted by **Ayodeji Williams** and mentored by: **Dr. Olabisi Ojo & Dr. Wanjun Hu**

A joint research between the **Department of Math, Computer Science and Physics & Department of Natural Sciences**

🏆 Awarded 3rd Place for Oral Presentation & 3rd Place for Poster Presentation at the 2026 ASU Research Symposium.

We used the research paper [LncDC](https://www.nature.com/articles/s41598-022-22082-7) by Li, M., & Liang, C. (2022) as the basis of our research.

The Github repo for LncDC:https://github.com/lim74/LncDC

### Dataset Sources & Analysis
**We used data sets from the Target-OS project**

**Files are above**

**1.** We used the python script "mRNA-Analysis1.py" to perform data analysis and visualization for
  
  --->[mRNA Data: os_target_gdc_clinical_data.tsv](https://master.cbioportal.org/study/summary?id=os_target_gdc) (in folders Data Analysis and Datasets respectively
  
**2.** We downloaded Augmented Star Gene data "gdc_download_20260303_162919.355206.tar.gz" and perfomed and lncRNA focused analysis of it with the python script "augmentedStarGeneAnalysis.py"

**3.** We also have another Augmented Star Gene data in "gdc_download_20260303_170024.760832.tar.gz" and performed an analysis(lncRNA focused) in "2ndaugmentedStarGeneAnalysis.py".

**4.** The third Augmented Star Gene data in "gdc_download_20260303_162919.355206.tar.gz" and the python script in "3rdaugmentedStarGeneAnalysis.py".

**5.** The fourth Augmented Star Gene data in "gdc_download_20260304_201550.898969.tar.gz" and the python script in "4thaugmentedStarGeneAnalysis.py".

**6.** The fifth Augmented Star Gene data in "gdc_download_20260304_202621.467167.tar.gz" and python script in "5thaugmentedStarGeneAnalysis.py"


----
----

## ML Model

The first ML model has been built in the folder "ML-Model" and it's named "os_lncrna_detector.pkl" the detection code is "1stMLModelDetCode.py"

Here is the url for the GUI for the model: https://miraclerna.streamlit.app/

----------------------------------------------------------------------------------------------------------------------------------------------------------------

**The contents below comes from the research proposal, they contains words like "will" implying our plans rather than action already taken, this README was written during the research process.**
## Abstract
Long non-coding RNAs (lncRNA), which constitute a large portion of the human transcriptome, are increasingly recognized as critical regulators in cancer biology, particularly in osteosarcoma (OS). Despite their significance, many lncRNAs remain poorly characterized due to limitations in traditional experimental approaches. The MIRACLE framework integrates machine learning–based coding potential tools to improve the identification and functional interpretation of lncRNAs in OS transcriptomic data. Analysis of five Augmented STAR-aligned RNA-seq datasets revealed a consistent lncRNA-dominant transcriptional profile, with lncRNAs representing approximately 28% of detected genes. Key transcripts, including MALAT1, NEAT1, NORAD, and H19, were highly expressed across samples, with MALAT1 reaching levels above 12,000 TPM in aggressive profiles. Significant proportions of reads mapped to ambiguous and unannotated regions, suggesting the presence of novel lncRNAs. Comparative patterns across datasets indicate variability in expression linked to tumor progression and cellular stress adaptation, with NEAT1 particularly elevated in advanced conditions. These findings support the hypothesis that lncRNAs serve as central regulators in OS pathogenesis, functioning through mechanisms such as microRNA sponging, chromatin remodeling, and stress-response modulation. Their high expression, tissue specificity, and correlation with disease stage position them as superior biomarkers compared to traditional protein-based markers. Integration with machine learning tools enhances the discovery of both annotated and novel lncRNAs, enabling scalable and reproducible analysis. Collectively, this work demonstrates the potential of computational genomics frameworks like MIRACLE to advance biomarker discovery, improve diagnostic precision, and provide accessible research workflows for undergraduate training in cancer genomics. 

## Introduction
Long non-coding RNAs (lncRNAs) are RNA molecules longer than 200 nucleotides that, while lacking protein-coding capacity, serve essential regulatory roles in gene expression and developmental pathways and have been implicated in a wide array of diseases including cancer and cardiovascular disorders (Hu et al, 2018). Remarkably, these non-coding sequences comprise nearly 98% of the human transcriptome, yet a significant portion remains functionally uncharacterized and their biological importance in specific disease contexts is only beginning to be uncovered. Osteosarcoma (OS), the most common malignant bone tumor in youth, has a particularly dire prognosis when metastasized, and mounting evidence suggests that aberrant lncRNA expression is both a hallmark of OS and a potential source for novel biomarkers and therapeutic targets (Ghert, et al, 2019; Li and Liang, 2022). However, the systematic identification and characterization of these lncRNAs have been restricted by the limitations of traditional laboratory approaches. Recent advances in high-throughput sequencing and computational biology have revolutionized this field, enabling the systematic discovery and classification of lncRNAs using efficient, alignment-free machine learning methods. Tools such as CPAT, which applies logistic regression to features like open reading frame (ORF) length and hexamer usage, and PLEK, which uses k-mer patterns with a support vector machine (SVM), exemplify this paradigm. CPC2 also leverages SVMs combined with ORF and physicochemical properties, while CPPred utilizes nucleotide composition and CTD (Composition, Transition, and Distribution) features for refined coding potential classification. LncFinder integrates ORF, sequence, secondary structure, and physicochemical information within an R package framework, and COME broadens the analytic landscape by incorporating both genomic attributes and experimental data (such as GC content and histone modifications) for comprehensive lncRNA detection (Li and Liang, 2022; Singh and Roy, 2022 ).This project therefore aims to review, compare, and apply computational tools to osteosarcoma RNA-seq datasets, establishing clear and reproducible protocols tailored for undergraduates at Albany State University. It will also provide a foundation for students and faculty to pursue similar research on other cancers and diseases. The systematic use of machine learning–based coding potential calculators is expected to uncover previously unrecognized, disease-relevant lncRNAs in osteosarcoma, thereby deepening our understanding of its transcriptomics and advancing the search for clinically important biomarkers and therapeutic targets.

## Hypothesis
Modern machine learning-based alignment-free coding potential assessment tools can accurately distinguish lncRNA transcripts from protein-coding mRNAs in large-scale transcriptomic datasets, enabling the discovery of novel functional non-coding RNAs in Osteosarcoma transcriptome.

## Objectives & Activities
Review and Comparison: We will systematically review and compare the underlying algorithms, features, strengths and weaknesses of CPAT, PLEK, CPC2, CPPred, LncFinder, and COME.
Dataset Preparation: We will curate and preprocess benchmark RNA-Seq transcript data, including both annotated lncRNAs and mRNAs, from publicly available sources.
Tool Implementation: We will apply each tool to predict lncRNAs in the dataset, document prediction workflows, and compare outputs.
Feature Analysis: We will investigate the biological relevance of core sequence features (e.g., ORF length, hexamer composition, secondary structure) used by these tools.
Results Synthesis: We will compile findings to assess the utility of these tools for large-scale genome annotation and potential biomarker discovery.
Dissemination: We will eventually present the research at Albany State Undergraduate Research Conference (and potentially other regional and national conferences), and prepare a manuscript for possible publication.

## Methods
The methods for this study will involve a comprehensive workflow integrating data collection, preprocessing, tool deployment, and comparative analysis. High-quality RNA-Seq datasets from osteosarcoma and control tissue will be curated from public repositories (National Cancer Institute, 2025). Transcript assembly will be performed using standard genome alignment tools, such as HISAT2 or STAR, followed by transcript reconstruction via StringTie and the filtering of candidate transcripts to retain those greater than 200 nucleotides and with at least two exons. Non-lncRNA sequences (tRNAs, rRNAs, and snoRNAs) will be removed by cross-referencing established sequence databases. The remaining transcript candidates will then be formatted for compatibility and analyzed using CPAT, PLEK, CPC2, CPPred, LncFinder, and COME, each of which applies+mimenidistinctive machine learning models and sequence- or structure-based feature extraction techniques for lncRNA prediction. Tool outputs for each transcript will be collated, and their predictions compared to annotated controls to evaluate accuracy, sensitivity, and F1-score (Li and Liang, 2022). Comparative feature importance and misclassification analysis will be conducted to understand each tool's strengths and limitations. All code, workflows, and findings will be documented within reproducible Jupyter Notebooks and released via GitHub to support undergraduate research training and transparency in research.

## Analysis plan (from research proposal)
The analysis plan will systematically compare the predictive performance of CPAT, PLEK, CPC2, CPPred, LncFinder, and COME using benchmark osteosarcoma and control RNA-Seq datasets. Each tool’s predictions for lncRNA versus protein-coding transcripts will be evaluated against ground truth annotations, calculating sensitivity, specificity, precision, F1-score, and overall accuracy for quantitative assessment. Further, receiver operating characteristic (ROC) and precision-recall (PR) curves will be generated to visualize performance metrics and discriminate thresholds across tools. Discrepancies in predictions will be examined through cross-tool Venn diagrams and feature importance ranking to highlight agreement and divergence, as well as to identify features most influential in classification outcomes. Instances of misclassification or discordant results will be tracked to characterize potential sources of error (e.g., transcript length, sequence ambiguity, GC content). The statistical analysis will be complemented by computational efficiency, runtime, and user experience observations, guiding practical recommendations for undergraduate research use. All findings, scripts, and data visualizations will be integrated into Jupyter Notebooks and openly shared for training, reproducibility, and transparent evaluation of tool performance in disease-relevant lncRNA discovery.


🤖 Technical Stack
Language: Python 3.12
Machine Learning: Scikit-Learn (Random Forest Classifier)
Bioinformatics: GDC Data Portal (TARGET-OS Dataset), LncDC
Data Processing: Pandas, NumPy, Joblib

## References
Ghert, M., Mwita, W., & Mandari, F. N. (2019). Primary Bone Tumors in Children and Adolescents Treated at a Referral Center in Northern Tanzania. Journal of the American Academy of Orthopaedic Surgeons. Global research & reviews, 3(3), e045. https://doi.org/10.5435/JAAOSGlobal-D-17-00045

Hu, G., Niu, F., Humburg, B. A., Liao, K., Bendi, S., Callen, S., Fox, H. S., & Buch, S. (2018). Molecular mechanisms of long noncoding RNAs and their role in disease pathogenesis. Oncotarget, 9(26), 18648–18663. https://doi.org/10.18632/oncotarget.24307

Li, M., & Liang, C. (2022). LncDC: a machine learning-based tool for long non-coding RNA detection from RNA-Seq data. Scientific reports, 12(1), 19083. https://doi.org/10.1038/s41598-022-22082-7

National Cancer Institute. (2025, Sept, 16). Therapeutically Applicable Research to Generate Effective Treatments (TARGET). https://www.cancer.gov/ccg/research/genome-sequencing/target

Singh, D., & Roy, J. (2022). A large-scale benchmark study of tools for the classification of protein-coding and non-coding RNAs. Nucleic acids research, 50(21), 12094–12111. https://doi.org/10.1093/nar/gkac1092
