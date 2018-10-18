##promoter SIDDs for 3 bacteria


library(ape)
library(Biostrings)


#promoter data
mycoplasma_gallisepticum_proms <- read.csv('/home/mikhail/Documents/mol_a_2017/mycoplasma_gallisepricum_fisunov_2016.csv')
colnames(mycoplasma_gallisepticum_proms)[1:2] <- c('strand', 'TSS')

mycoplasma_gallisepticum_genome <- paste0(read.GenBank(access.nb = 'CP006916.3', as.character = T)[[1]], collapse = '')[[1]]
str(mycoplasma_gallisepticum_genome)

#one promoter is at the right flank, need to extend
extended_mycoplasma_gallisepticum_genome <- paste0(mycoplasma_gallisepticum_genome,
                                                       substr(mycoplasma_gallisepticum_genome, 1, 1000))

#sequences 

mycoplasma_gallisepticum_proms_seqs <- c()
for (i in seq_along(mycoplasma_gallisepticum_proms$strand)){
  if (mycoplasma_gallisepticum_proms$strand[i] == '1') {
    mycoplasma_gallisepticum_proms_seqs <- c(mycoplasma_gallisepticum_proms_seqs, substr(extended_mycoplasma_gallisepticum_genome, start = mycoplasma_gallisepticum_proms$TSS[i]-750, mycoplasma_gallisepticum_proms$TSS[i]+750 ))
  } else {
    mycoplasma_gallisepticum_proms_seqs <- c(mycoplasma_gallisepticum_proms_seqs, reverse(substr(extended_mycoplasma_gallisepticum_genome, start = mycoplasma_gallisepticum_proms$TSS[i]-750, mycoplasma_gallisepticum_proms$TSS[i]+750 )))
  }
}

table(sapply(mycoplasma_gallisepticum_proms_seqs, nchar))

#sequence statistics
par(mfrow = c(1,2))
mycoplasma_gallisepticum_proms_gc <- sapply(mycoplasma_gallisepticum_proms_seqs, FUN = function(x) {GC(strsplit(x, split = '')[[1]])})
boxplot(mycoplasma_gallisepticum_proms_gc)

mycoplasma_gallisepticum_proms_gc1 <- sapply(mycoplasma_gallisepticum_proms_seqs, FUN = function(x) {GC1(strsplit(x, split = '')[[1]])})
boxplot(mycoplasma_gallisepticum_proms_gc1)

mycoplasma_gallisepticum_proms_6mers <- oligonucleotideFrequency(DNAString(paste0(mycoplasma_gallisepticum_proms_seqs, collapse = '')[[1]]), 6, as.prob = T)
barplot(sort(mycoplasma_gallisepticum_proms_6mers, decreasing = T)[1:20], las = 2)




dir_sidd_mycoplasma <- ('/home/mikhail/Documents/mol_a_2017/mycoplasma_sidds')
dir.create(dir_sidd_mycoplasma)
setwd(dir_sidd_mycoplasma)


# SIDDr function
rSIDD <- function(seq, algorithm = algorithm, sigma = sigma, form = form, file.out = file.out){
  library(seqinr)
  sidd_dir <- paste0(getwd(), '/')
  seq <- paste0(seq[[1]], collapse = '')
  #  sidd_dir <- paste0(getwd(), "/TMPSIDD", as.numeric(Sys.time()))
  # dir.create(sidd_dir)
  writeLines(seq, con = substr(seq, 0, 25))
  print(paste('Processing', substr(seq, 0, 25)))
  system(paste0('cd ', '/home/mikhail/Documents/sist/
        ',
                'perl -X master.pl -a M -f ', paste0(sidd_dir, substr(seq, 0, 25)),' -s ', sigma, ' -o ' ,  paste0(sidd_dir, file.out,'sidd_output.tsv')))
  tmp <- read.csv(paste0(sidd_dir, file.out,'sidd_output.tsv'), sep = '\t')
  # #
  tmp <- cbind(tmp, strsplit(seq, split = '')[[1]])
  colnames(tmp)[4] <- 'Nucleotide'
  return(tmp)
  unlink(substr(seq, 0, 25))
}

#check
seq <- mycoplasma_gallisepticum_proms_seqs[1]
tmp <- rSIDD(seq, algorithm = 'M', sigma = 0.06, form = 'l', file.out = 'tmpjunk.csv')

#SIDD calculation

mycoplasma_gallisepticum_proms_sidds <- c()

for (i in mycoplasma_gallisepticum_proms_seqs) {
  print(i)
  tmp <- rSIDD(i, algorithm = 'M', sigma = 0.06, form = 'l', file.out = 'tmpjunk.csv')
  mycoplasma_gallisepticum_proms_sidds <- cbind(mycoplasma_gallisepticum_proms_sidds, tmp$P.x.)
}

str(mycoplasma_gallisepticum_proms_sidds)
#i = mycoplasma_gallisepticum_proms_seqs[351]
#reading in 2 other bacteria SIDDs

load('/home/mikhail/Documents/mol_a_2017/sidds_ecoli.rda')
load('/home/mikhail/Documents/mol_a_2017/sidds_synechocystis.rda')



#function to cut heirarchical clusterization results at given clusters numer and plot all pfrofiles from each

cuplcl <- function(data, method, k) {
  library(fastcluster)
  library(dendextend)
  hclust <- hclust.vector(data, method = method)
  cut_hclust<-cutree(hclust, k=k)
  m <- rbind(rep(1, k), c(2:(k+1)))
  layout(m)
  par(mar = c(1, 1, 0, 0))
  
  dend <- color_branches(as.dendrogram(hclust), col = cut_hclust[order.dendrogram(as.dendrogram(hclust))])
  plot(dend%>%set('labels_cex', NA))
  
  for (i in unique(cut_hclust[order.dendrogram(as.dendrogram(hclust))])) {
    matplot(t(data[which(cut_hclust==i),]), type='l', col=i)
    
    print(c(i, names(which(cut_hclust==i))))
  }
}

#analysis itself

cuplcl(t(sidds_ecoli[751:2251,]), method = 'ward', k = 4)
cuplcl(t(sidds_syn), method = 'ward', k = 4)
cuplcl(t(mycoplasma_gallisepticum_proms_sidds), method = 'ward', k = 4)

cuplcl(t(cbind(sidds_ecoli[751:2251,], sidds_syn, mycoplasma_gallisepticum_proms_sidds)), method = 'ward', k = 4)
       

###vlhA non-complete list 1000 nts
library(seqinr)
library(ape)
vlhA_1000 <- read.fasta('/home/mikhail/Documents/mycoplasma/gaa_for_biophysics1000(1).fasta', as.string = T, seqonly = T)
str(vlhA_1000)
image.DNAbin(as.DNAbin(t(sapply(vlhA_1000, FUN = function(x) {strsplit(x, split = '')[[1]]}))))

dir_sidd_vlhA <- ('/home/mikhail/Documents/mycoplasma/vlhA_sidd')
unlink(dir_sidd_vlhA, recursive = T)
dir.create(dir_sidd_vlhA)
setwd(dir_sidd_vlhA)


# SIDDr function
rSIDD <- function(seq, algorithm = algorithm, temperature = temperature, sigma = sigma, form = form, file.out = file.out){
  library(seqinr)
  sidd_dir <- paste0(getwd(), '/')
  seq <- paste0(seq[[1]], collapse = '')
  #  sidd_dir <- paste0(getwd(), "/TMPSIDD", as.numeric(Sys.time()))
  # dir.create(sidd_dir)
  writeLines(seq, con = substr(seq, 0, 25))
  print(paste('Processing', substr(seq, 0, 25)))
  system(paste0('cd ', '/home/mikhail/Documents/sist/
                ',
                'perl -X master.pl -a M -f ', paste0(sidd_dir, substr(seq, 0, 25)), ' -t ', temperature, ' -s ', sigma, ' -o ' ,  paste0(sidd_dir, file.out,'sidd_output.tsv')))
  tmp <- read.csv(paste0(sidd_dir, file.out,'sidd_output.tsv'), sep = '\t', skip = 1)
  # #
  tmp <- cbind(tmp, strsplit(seq, split = '')[[1]])
  colnames(tmp)[4] <- 'Nucleotide'
  return(tmp)
  unlink(substr(seq, 0, 25))
}


vlhA_sidds <- c()

for (i in vlhA_1000) {
  print(i)
  tmp <- rSIDD(i, algorithm = 'M', temperature = 314, sigma = 0.06, form = 'l', file.out = 'tmpjunk.csv')
  vlhA_sidds <- cbind(vlhA_sidds, tmp$P.x.)
}

str(vlhA_sidds)
matplot(vlhA_sidds, type = 'l', col = 1, ylim = 0:1, lwd = 0.3)

#averaged values
plot(apply(vlhA_sidds, MARGIN = 1, FUN = mean), type = 'l')
lines(apply(vlhA_sidds, MARGIN = 1, FUN = median), type = 'l', col = 'red')
source('https://raw.githubusercontent.com/FVortex/different_bacteria_physics_and_text/master/cuplcl.R')

cuplcl(t(vlhA_sidds), method = 'ward', k =2)
abline(v = 1:10*100, lty =2)


#############non-vlhA
###vlhA non-complete list 1000 nts
library(seqinr)
library(ape)
nonvlhA_1000 <- read.fasta('/home/mikhail/Documents/mycoplasma/S6_promoters_all_genes.fasta', as.string = T, seqonly = T)

nonvlhA_1000 <- sapply(nonvlhA_1000, FUN = function(x) {substr(x, start = 500, stop = 1500)})
str(nonvlhA_1000)
image.DNAbin(as.DNAbin(t(sapply(nonvlhA_1000, FUN = function(x) {strsplit(x, split = '')[[1]]}))))

dir_sidd_nonvlhA <- ('/home/mikhail/Documents/mycoplasma/nonvlhA_sidd')
unlink(dir_sidd_nonvlhA, recursive = T)
dir.create(dir_sidd_nonvlhA)
setwd(dir_sidd_nonvlhA)


nonvlhA_sidds <- c()

for (i in nonvlhA_1000) {
  print(i)
  tmp <- rSIDD(i, algorithm = 'M', temperature = 314, sigma = 0.06, form = 'l', file.out = 'tmpjunk.csv')
  nonvlhA_sidds <- cbind(nonvlhA_sidds, tmp$P.x.)
}

str(nonvlhA_sidds)
matplot(nonvlhA_sidds, type = 'l', col = 1, ylim = 0:1, lwd = 0.3)

#averaged values
plot(apply(nonvlhA_sidds, MARGIN = 1, FUN = mean), type = 'l')
lines(apply(nonvlhA_sidds, MARGIN = 1, FUN = median), type = 'l', col = 'red')
source('https://raw.githubusercontent.com/FVortex/different_bacteria_physics_and_text/master/cuplcl.R')

cuplcl(t(nonvlhA_sidds), method = 'ward', k =2)
abline(v = 1:10*100, lty =2)



###normality test for both
shapiro.test(apply(vlhA_sidds, MARGIN = 1, FUN = max))
shapiro.test(apply(nonvlhA_sidds, MARGIN = 1, FUN = max))


vlhA_maxs <- (apply(vlhA_sidds, MARGIN = 1, FUN = max))
nonvlhA_maxs <- (apply(nonvlhA_sidds, MARGIN = 1, FUN = max))


wilcox.test(c(vlhA_maxs, nonvlhA_maxs) ~ c(rep('v', length(vlhA_maxs)), rep('n', length(nonvlhA_maxs))))

###

par(mfrow = c(1,2))
boxplot(vlhA_maxs, main = 'vlhA promoters from various strains', ylim = c(0,1.1), frame = F, xaxt = 'n', ylab = 'SIDD profile maximum for promoters')
text(x = c(1,2), y = c(1, 1), labels = c('***','***'),cex=2)

boxplot(nonvlhA_maxs, main = 'S6 strain promoters', ylim = c(0,1.1), frame = F, xaxt = 'n', yaxt = 'n')
#text(x = c(1,2), y = c(1, 1), labels = c('***','***'),cex=2)

###matplots
par(mfrow = c(2,1))
`Sequence (nts)` <- -300:300
matplot(`Sequence (nts)`, vlhA_sidds[260:860,], type = 'l', col = 1, ylim = 0:1, lwd = 0.1, ylab = 'Opening probability', main = 'vlhA promoters from various strains')
matplot(`Sequence (nts)`, nonvlhA_sidds[260:860,], type = 'l', col = 1, ylim = 0:1, lwd = 0.1, ylab = 'Opening probability', main = 'S6 strain promoters')
