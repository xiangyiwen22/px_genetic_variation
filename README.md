### Genotyping-by-Sequencing for the Analysis of the Genetic Variation of *Podosphaera xanthii*, Incitant of Cucurbit Powdery Mildew 
Study genetic variation of 109  *Podospheara xanthii* isolates, which causing cucurbits powdery mildew worldwide <br>

- `map2snp109.csv` contains SNPs data of 109 Px and their collection information including hosts and geographical locations. 
- `GBS_git.r` analyze and create tables and figures for this project including 
  - `elbow.png`: Plot "elbow method" to compute within cluster variance for K-means analysis
  - `dend.png`: Dendrogram of hierarchical clustering with color coded bars of isolates from different hosts, locations and clustering groups by K-means analysis
  - `host_fst` and `loc_fst`: Pairwise Fst calculation for isolates from different populations, host wise and location wise, respectively. 
