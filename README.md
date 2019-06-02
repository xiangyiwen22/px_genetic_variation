# Genotyping-by-Sequencing for the Analysis of the Genetic Variation of *Podosphaera xanthii*, Incitant of Cucurbit Powdery Mildew 
Study genetic variation of 109  *Podospheara xanthii* isolates, which causing cucurbits powdery mildew worldwide <br>

- `map2snp109.csv` contains SNPs data of 109 Px and their collection information including hosts and geographical locations. 
- `GBS_git.r` contains R code analyze and create tables and figures for this project including 
  - Plot "elbow method" to compute within cluster variance for K-means analysis
  - Dendrogram of hierarchical clustering with color coded bars of isolates from different hosts, locations and clustering groups by K-means analysis
  - Pairwise Fst calculation for isolates from different populations, host wise and location wise, respectively. 
