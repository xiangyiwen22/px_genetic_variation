library(adegenet) 
library(hierfstat)
library(RColorBrewer)
library(dendextend)

#convert SNPs from hapMap into 0,1,and 2 format 
#( 0 is homozygote, 1 is heterozygote and 2 # is another type of homozygote.)
hapMap2genlight <- function(file){
  require(adegenet)
  require(parallel)
  hapmap <- read.table(file, header=TRUE, row.names=1, sep="\t",
                       stringsAsFactors=FALSE)[,-(2:10)]
  samples <- names(hapmap)[-1]
  loci <- row.names(hapmap)
  
  s <- as.integer(c(0,1,2,NA))
  ac <- s
  ag <- s
  at <- s
  cg <- s
  ct <- s
  gt <- s
  names(ac) <- c("A","M","C","N")
  names(ag) <- c("A","R","G","N")
  names(at) <- c("A","W","T","N")
  names(cg) <- c("C","S","G","N")
  names(ct) <- c("C","Y","T","N")
  names(gt) <- c("G","K","T","N")
  conv <- list(ac,ac,ag,ag,at,at,cg,cg,ct,ct,gt,gt)
  names(conv) <- c("A/C","C/A","A/G","G/A","A/T","T/A","C/G","G/C",
                   "C/T","T/C","G/T","T/G")
  
  S <- length(samples)
  SBlist <- vector(mode="list",S)   
  for(i in 1:S){
    mygen <- mapply(function(type,gen) unname(conv[[type]][gen]),
                    type=hapmap[[1]], gen=hapmap[[i+1]],
                    SIMPLIFY=TRUE, USE.NAMES=FALSE)
    
    SBlist[[i]] <- new("SNPbin", mygen)
  }
  
  
  x <- new("genlight", SBlist)
  locNames(x) <- loci
  indNames(x) <- samples
  
  return(x)
}

####### Read dataset and cleanup########
# 109 isolates with information of collections 
# (iso names, hosts, location of collection)
map2snp109<-read.csv("./raw_data/map2snp109.csv", header = TRUE)

# fill missing hosts as unknown
# "^$" represents empty string "" (in regular expression ^-start, $-end) 
map2snp109[0:nrow(map2snp109), 4] <- sub("^$", "unknown", 
                                         map2snp109[0:nrow(map2snp109), 4]) 

# fill missing snps with col means 
for(i in 9:ncol(map2snp109)){
  map2snp109[is.na(map2snp109[,i]), i] <- mean(map2snp109[,i], na.rm = TRUE)
}

# subset snps cols with only
px_snps = map2snp109[0:nrow(map2snp109), 9:ncol(map2snp109)] 

##### Table 1. Count freq of each location ##### 
library(plyr)
isolate109 <- reshape(count(map2snp109, c("Location", "Host")), 
                      direction = "wide", idvar = "Host", timevar = "Location")
# replace NA with 0
isolate109[is.na(isolate109)] <- 0
# rename column names
colnames(isolate109)<-c('Host \ Location', 'CA', 'IL', 'IN', 'Italy', 'MI', 'NY', 'TX', 'WA', 'WI')
write.csv(isolate109, file = "./results/table1.csv") 

####### Table 2 FST calculation location ########
# convert df to genind format
px_genind = df2genind(px_snps, ploidy=1, ncode=1) # no isolate names 

# assign location as populations 
pop(px_genind)<-map2snp109[0:nrow(map2snp109), 5] 

# check freq of each location
as.data.frame(table(map2snp109[0:nrow(map2snp109), 5]))

#location fst calculation
#pairwise calculation takes a while to run
loc_fst <- pairwise.fst(px_genind, res.type="matrix") 

write.csv(loc_fst, file = "./results/loc_fst.csv") 

# Simpson
library("poppr")
# ref https://grunwaldlab.github.io/Population_Genetics_in_R/Genotypic_EvenRichDiv.html
# convert from genind to genclone, after assign location as populations 
px_genclone <- as.genclone(px_genind)

# lambda column calculated Simpsons index
px_diversity <- poppr(px_genclone)
write.csv(px_diversity, file = "./results/px_diversity.csv") 


####### Table 3 FST calculation host ########
# convert df to genind format
px_genind = df2genind(px_snps, ploidy=1, ncode=1) # no isolate names 

# assign hosts as populations 
pop(px_genind)<-map2snp109[0:nrow(map2snp109), 4] 

# check freq of each host
as.data.frame(table(map2snp109[0:nrow(map2snp109), 4]))

# fst calculation, host as population
host_fst <- pairwise.fst(px_genind, res.type="matrix") 
write.csv(host_fst, file = "host_fst.csv") 

####### Figure 2. Kmeans elbow test #######
#Elbow Method for finding the optimal number of clusters
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
wss <- sapply(1:k.max, function(k){kmeans(px_snps, k, nstart=50,iter.max = 15)$tot.withinss})
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


##### Figure 1. Dendrogram ##### 
# assign row names: isolate name - host
rownames(px_snps)<-do.call(paste,c(map2snp109[c("Row.names","Host")],sep="-"))
# check dimension
dim(px_snps)
# create dendrogram
dend<-as.dendrogram(hclust(dist(px_snps)))

#### add host color bar 
# replicate a vector of host names with Unknown, replace it with host name one by one
host_ab <- rep("Unknown", length(rownames(px_snps)))
# Create a vector giving a color for each isolate to which host it belongs to
host_list <- c("Unknown", "Cucurbita pepo", "Cucurbita maxima", "Cucumis melo", "Cucumis sativus", "Cucurbita moschata", "Lagenaria siceraria")

for (host in host_list) {
  # find index with this host
  hos_X <- grepl(chartr(old = " ", new = "_", host), rownames(px_snps))
  # replace it with this host name
  host_ab[hos_X] <- host
}

host_ab <- factor(host_ab)
n_host_ab <- length(unique(host_ab))
cols_7<-brewer.pal(7, "Dark2")
col_host_ab <- cols_7[host_ab]

####  add location color bar 
# Create a vector giving a color for each isolate to which location it belongs to
state_ab <- rep("", length(rownames(px_snps)))

loc_list<-c('CA', 'IL', 'IN', 'Italy', 'MI', 'NY', 'TX', 'WA', 'WI')
for (loc in loc_list) {
  # find index with this loc
  loc_X <- grepl(loc, rownames(px_snps))
  state_ab[loc_X] <- loc
}

state_ab <- factor(state_ab)
n_state_ab <- length(unique(state_ab))
cols_11<-brewer.pal(11, "Set3")
col_state_ab <- cols_11[state_ab]

# color labels by state:
labels_colors(dend) <- col_state_ab[order.dendrogram(dend)]

# start building tree
# showing the various clusters cuts 
k234 <- cutree(hclust(dist(px_snps)), k = 3:7)
# color branches based on cutting the tree into 7 clusters:
dend1 <- color_branches(dend, k = 7)

# set the edge very large c(bottom, left, top, right) 
# xpd=TRUE allow legend outside of the plot
par(mar = c(4,4,10,0), xpd=TRUE)

# show plot
plot(dend1, ylab = "Distance", leaflab="none", bty='L')

#customize color bar
#col_7 <- brewer.pal(7, "Pastel1")
col_7<-rainbow(7, s = 0.5)
col_cluster7<-col_7[k234[,5]]
col_cluster6<-col_7[k234[,4]]
col_cluster5<-col_7[k234[,3]]
col_cluster4<-col_7[k234[,2]]
col_cluster3<-col_7[k234[,1]]

colored_bars(cbind(col_cluster7, col_cluster6,col_cluster5,
                   col_cluster4,col_cluster3, 
                   col_host_ab, col_state_ab), 
             dend1, 
             rowLabels = c(paste0("k=", 7:3), "host", 'state'),
             y_shift=2, cex.rowLabels=1, text_shift=0.8)

# customize legend pos
legend(title="state",
       title.adj =0.01, 
       x = 1, y = 68, 
       cex =0.6, ncol=6,
       legend = as.character(unique(state_ab)),
       fill = unique(col_state_ab))

legend(title = "host", 
       text.font = 3, 
       title.adj =0.01, 
       x = 1, y = 78,
       cex =0.6, ncol=3,
       legend = as.character(unique(host_ab)),
       fill = unique(col_host_ab)
)

