### Genotyping-by-Sequencing for the Analysis of the Genetic Variation of *Podosphaera xanthii*, Incitant of Cucurbit Powdery Mildew

#### Introduction
This research was conducted to identify species causing powdery mildew on cucurbits and to determine genetic variations among isolates of the pathogen. We collected 109 isolates from six cucurbit species hosts (*Cucumis melo*, *Cucumis sativus*, *Cucurbita maxima*, *Cucurbita moschata*, *Cucurbita pepo*, and *Lagenaria siceraria*) in California, Illinois, Indiana, Michigan, New York, Texas, Washington, and Wisconsin in the United States and in Italy. By sequencing the internal transcribed spacer region of the nuclear rDNA of these 109 isolates, Podosphaera xanthii was found as the only species causing powdery mildew on cucurbits in the United States. Genotyping-by-sequencing was applied to these 109 isolates to investigate their genetic diversity, which showed a trend of isolates clustering from New York and Italy. 

#### Dataset
- [`map2snp109.csv`](/map2snp109.csv) contains SNPs data of 109 Px and their collection information including hosts and geographical locations. 

#### Analysis
- [`GBS_git.r`](/GBS_git.r) analyzes and creates tables and figures for this project including calculate genetic similarity using Simpson's index and Fst index, and plot visualization using clustering plot. 


#### Conclusion
We found that isolates collected in Illinois had high within-sample genetic diversity, which was shown by both clustering analysis and Simpson’s index values (Simpson 1949), but there was higher genetic similarity between Illinois and other locations based on Fst calculations

##### Fig. 1. Dendrogram of hierarchical clustering with color coded bars of isolates from different hosts, locations and clustering groups by K-means analysis
![GitHub Logo](/results/dend.png) <br>
Seven distinct groups (from left to right: red, brown, light green, green, light blue, blue, and pink) of 109 *Podosphaera xanthii* isolates with 2,266 single nucleotide polymorphisms, identified by hierarchical clustering on the distance matrix. Different groupings are shown for different values of k (k = 3 to 7). *P. xanthii* isolates were collected from six known hosts (*Cucumis melo*, *Cucumis sativus*, *Cucurbita maxima*, *Cucurbita moschata*, *Cucurbita pepo*, and *Lagenaria siceraria*) and one unknown cucurbit host in California (CA), Illinois (IL), Indiana (IN), Michigan (MI), New York (NY), Texas (TX), Washington (WA), Wisconsin (WI), and Italy with 16, 58, 3, 9, 3, 8, 6, 2, and 4 isolates, respectively.

![Table 2.](/results/table2.png)<br>
Combined results from raw output [fst index](/results/loc_fst.csv) and [Simpson's index](/results/px_diversity.csv) <br>
![Table 3.](/results/table3.png)<br>
From raw output [fst index calculated based on host](/results/host_fst.csv)<br>
