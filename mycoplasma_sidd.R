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
tmp <- rSIDD(seq, algorithm = 'M', sigma = 0.09, form = 'l', file.out = 'tmpjunk.csv')

#SIDD calculation

mycoplasma_gallisepticum_proms_sidds <- c()

for (i in mycoplasma_gallisepticum_proms_seqs) {
  print(i)
  tmp <- rSIDD(i, algorithm = 'M', sigma = 0.09, form = 'l', file.out = 'tmpjunk.csv')
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
library(stringr)
library(Biostrings)
vlhA_1000 <- read.fasta('/home/mikhail/Documents/mycoplasma/gaa_for_biophysics1000.fasta', as.string = T, seqonly = T)
str(vlhA_1000)
image.DNAbin(as.DNAbin(t(sapply(vlhA_1000, FUN = function(x) {strsplit(x, split = '')[[1]]}))))


str_count(vlhA_1000[[1]], pattern = 'GAAGAA')
counts <- sapply(vlhA_1000, FUN = function(x) {str_count(x, pattern = 'GAAGAA')})
table(counts)

ord <- order(counts)

image.DNAbin(as.DNAbin(t(sapply(vlhA_1000[ord], FUN = function(x) {strsplit(x, split = '')[[1]]}))))
#!half of promoters is reverse complement

vlhA_1000_reverse <- sapply(vlhA_1000[which(counts ==0)], FUN = function(x) {as.character(reverseComplement(DNAString(x)))})
vlhA_1000_forward <- sapply(vlhA_1000[which(counts >0)], FUN = function(x) {x})
table(sapply(vlhA_1000_reverse, nchar))
table(sapply(vlhA_1000_forward, nchar))


vlhA_1000_fixed <- c(vlhA_1000_forward, vlhA_1000_reverse)
image.DNAbin(as.DNAbin(t(sapply(vlhA_1000_fixed, FUN = function(x) {strsplit(x, split = '')[[1]]}))))

str_count(vlhA_1000_fixed[[1]], pattern = 'GAAGAA')
counts <- sapply(vlhA_1000_fixed, FUN = function(x) {str_count(x, pattern = 'GAAGAA')})
table(counts)

ord <- order(counts)

dir_sidd_vlhA <- ('/home/mikhail/Documents/mycoplasma/vlhA_sidd')
dir.create(dir_sidd_vlhA)
setwd(dir_sidd_vlhA)


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
  tmp <- rSIDD(i, algorithm = 'M', sigma = 0.09, form = 'l', file.out = 'tmpjunk.csv')
  vlhA_sidds <- cbind(vlhA_sidds, tmp$P.x.)
}

str(vlhA_sidds)
matplot(vlhA_sidds, type = 'l', col = 1, ylim = 0:1)

#averaged values
plot(apply(vlhA_sidds, MARGIN = 1, FUN = mean), type = 'l')
lines(apply(vlhA_sidds, MARGIN = 1, FUN = median), type = 'l', col = 'red')
source('https://raw.githubusercontent.com/FVortex/different_bacteria_physics_and_text/master/cuplcl.R')

cuplcl(t(vlhA_sidds), method = 'ward', k =4)


### non-vlhA list 1000 nts
library(seqinr)
library(ape)
nonvlhA_1000 <- read.fasta('/home/mikhail/Documents/mycoplasma/S6_promoters_all_genes.fasta', as.string = T, seqonly = T)
nonvlhA_1000 <- sapply(nonvlhA_1000, FUN = function(x) {substr(x, start = 500, stop = 1500)})

str(nonvlhA_1000)
image.DNAbin(as.DNAbin(t(sapply(nonvlhA_1000, FUN = function(x) {strsplit(x, split = '')[[1]]}))))

dir_sidd_nonvlhA <- ('/home/mikhail/Documents/mycoplasma/nonvlhA_sidd')
dir.create(dir_sidd_nonvlhA)
setwd(dir_sidd_nonvlhA)


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
  tmp <- read.csv(paste0(sidd_dir, file.out,'sidd_output.tsv'), sep = '\t', skip = 1)
  # #
  tmp <- cbind(tmp, strsplit(seq, split = '')[[1]])
  colnames(tmp)[4] <- 'Nucleotide'
  return(tmp)
  unlink(substr(seq, 0, 25))
}


nonvlhA_sidds <- c()

for (i in nonvlhA_1000) {
  print(i)
  tmp <- rSIDD(i, algorithm = 'M', sigma = 0.09, form = 'l', file.out = 'tmpjunk.csv')
  nonvlhA_sidds <- cbind(nonvlhA_sidds, tmp$P.x.)
}

str(nonvlhA_sidds)
matplot(nonvlhA_sidds, type = 'l', col = 1, ylim = 0:1)

#averaged values
plot(apply(nonvlhA_sidds, MARGIN = 1, FUN = mean), type = 'l')
lines(apply(nonvlhA_sidds, MARGIN = 1, FUN = median), type = 'l', col = 'red')


source('https://raw.githubusercontent.com/FVortex/different_bacteria_physics_and_text/master/cuplcl.R')

cuplcl(t(nonvlhA_sidds), method = 'ward', k =4)


#vlhA vs non-vlhA statistics

par(mfrow = c(1,2))
boxplot(apply(vlhA_sidds, MARGIN = 1, FUN = max), main = 'vlhA promoters, various strains', ylab = 'SIDD profile maximum for a sequence', ylim = 0:1)
boxplot(apply(onvlhA_sidds, MARGIN = 1, FUN = max), main = 'S strain promoters', ylim = 0:1)


###Clearing workspace, setting working directory and loading R packages{r, message=FALSE}
rm(list=ls())

#setwd('/home/mikhail/Documents/Script_2016_all/PCA/')

library(R.matlab)
library(seqinr)
library(factoextra)
#library(data.table)
library(caret)
#library(doMC)
library(Biostrings)
library(ape)
library(reldna)
library(fastcluster)
library(dendextend)
library(ggplot2)
library(ggseqlogo)



#r strReverse funtcion, message = FALSE}
## a useful function: rev() for strings
strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

if(!require(ape)){stop('Code requires library "ape".\n')}
if(!require(reldna)){stop('Code requires library "reldna".\n')}
if(!require(parallel)){stop('Code requires library "parallel".\n')}
if(!require(Biostrings)){stop('Code requires library "Biostrings".\n')}

get_forward_subseq <- function(seq, tss, boundaries) {
  return(as.character(
    DNAString(seq,
              start = tss-boundaries[1],
              nchar = sum(boundaries)+1)))
}

get_reverse_subseq <- function(seq, tss, boundaries) {
  return(as.character(
    reverseComplement(
      DNAString(seq,
                start = tss-boundaries[2],
                nchar = sum(boundaries)+1))))
}

get_forward_substring <- function(seq, tss, boundaries) {
  return(seq[(tss-boundaries[1]):(tss+boundaries[2])])
}

get_reverse_substring <- function(seq, tss, boundaries) {
  res <- rev(chartr("ATGC", "TACG", seq[(tss-boundaries[2]):(tss+boundaries[1])]))
  return(res)
}

cut_char_vec <- function(char_vec, tss, boundaries, strand) {
  res <- switch(strand,
                forward = char_vec[(tss-boundaries[1]):(tss+boundaries[2])],
                reverse = char_vec[(tss-boundaries[2]):(tss+boundaries[1])]
  )
  return(res)
}

#' Title
#'
#' @param seq full sequence (chromosome)
#' @param interval requested interval
#' @param tss
#' @param strand DNA strand
#'
#' @return
#' @export
#'
#' @examples
dynchars <- function(seq, average_interval_size, tss, boundaries, strand=c('forward','reverse')) {
  
  if (missing(average_interval_size))
    stop("Need to specify average_interval_size")
  
  if(!is.character(seq))
    stop("Sequence must be a character vector containing A, C, G, T letters only")
  #make sure no small letters are present
  seq <- toupper(seq)
  # if((tss - average_interval_size < 1) || (tss + average_interval_size > length(seq)))
  #   stop("Considered average_interval_size exceeds the size of the sequence")
  
  # seq <- switch(strand, 
  #               forward = seq,
  #               reverse = as.character(reverseComplement(DNAString(seq))))
  # seq <- unlist(strsplit(seq, ''))
  # seq <- toupper(seq)
  
  boundaries <- boundaries + average_interval_size %/% 2
  
  seq <- cut_char_vec(seq, tss, boundaries, strand)
  
  # strand <- match.arg(strand)
  # seq <- switch(strand,
  #             forward=get_forward_substring(seq, tss, boundaries),
  #             reverse=get_forward_substring(seq, tss, boundaries))#get_reverse_substring(seq, tss, boundaries))
  
  a<-3.4*10^(-10)
  #I<-c(7.6, 4.8, 8.2, 4.1)*10^(-44)
  K<-c(227, 155, 220, 149)*10^(-20)
  V<-c(2.09, 1.43, 3.12, 2.12)*10^(-20)
  
  csA<-cumsum(seq=='A')
  csT<-cumsum(seq=='T')
  csG<-cumsum(seq=='G')
  csC<-cumsum(seq=='C')
  
  countA = csA[average_interval_size:length(csA)]-c(0, csA[1:(length(csA)-average_interval_size)])
  countT = csT[average_interval_size:length(csT)]-c(0, csT[1:(length(csT)-average_interval_size)])
  countG = csG[average_interval_size:length(csG)]-c(0, csG[1:(length(csG)-average_interval_size)])
  countC = csC[average_interval_size:length(csC)]-c(0, csC[1:(length(csC)-average_interval_size)])
  
  M <- cbind(countA, countT, countG, countC)/average_interval_size
  
  #Is <- as.numeric(M%*%I)#! numeric conversion
  Ks <- as.numeric(M %*% K)
  Vs <- as.numeric(M %*% V)
  
  E0 <- (8*(Ks*Vs)^0.5)* 6E23 / 4184
  d <- ((Ks*a^2)/Vs)^(0.5)/a;
  gc = M[,3] + M[,4]
  
  dynchars_return <- list(E0=E0, d=d, gc=gc)
  return(dynchars_return)
}

#' Title
#'
#' @param tss_position
#' @param extended_string -- vector of characters
#' @param ep_interval
#' @param zout
#' @param strand
#'
#' @return
#' @export
#'
#' @examples
calculate_EP_on_interval <- function(tss_position, extended_string, ep_interval=c(270, 130), zout=-540:179, strand=c('forward','reverse')) {
  #lseqspline1D(substr(e.coli_U00096.2, exp_tsss[i]-250, exp_tsss[i]+150), bound=c(50, 350), ref=251 )
  # strand<- match.arg(strand)
  #extended_string <- unlist(strsplit(extended_string, ''))
  # subseq <-switch(strand,
  #                 forward = get_forward_substring(extended_string, tss_position, ep_interval),
  #                 reverse = get_reverse_substring(extended_string, tss_position, ep_interval))
  # subseq <- paste0(subseq, collapse = "")
  
  # subseq <- switch(strand,
  #                  forward = get_forward_subseq(extended_string, tss_position, ep_interval),
  #                  reverse = get_reverse_subseq(extended_string, tss_position, ep_interval))
  # subseq <- switch(strand,
  #                  forward = get_forward_substring(extended_string, tss_position, ep_interval),
  #                  reverse = get_reverse_substring(extended_string, tss_position, ep_interval))
  # subseq <- as.character(subseq)
  subseq <- as.character(cut_char_vec(extended_string, tss_position, ep_interval, strand))
  subseq <- paste0(subseq, collapse = "")
  p <- lseqspline1D(
    subseq,
    bound=c(50, 350),
    ref=271)
  return(p$mpot[p$x %in% zout])
}

#' Title
#'
#' @param dnaSeq
#' @param dynamic_interval
#' @param tss_position
#' @param ep_interval
#' @param zout
#' @param strand
#'
#' @return
#' @export
#'
#' @examples
calculate_profile <- function(dnaSeq, average_interval_size, tss_position, dynamic_interval, ep_interval=c(270, 130), zout=-540:179, strand=c('forward','reverse')) {
  #ep_interval = c(163, 213) for reverse, for forward c(267, 217), zout - the same
  #seq, average_interval_size, tss, boundaries, strand=c('forward','reverse')
  #dnaSeq <- unlist(strsplit(dnaSeq, ''))
  dynRes <- dynchars(dnaSeq, average_interval_size, tss_position,  dynamic_interval, strand)
  epRes <- calculate_EP_on_interval(tss_position, dnaSeq, ep_interval, zout, strand)
  return(c(dynRes$E0, dynRes$d, epRes, dynRes$gc))
}
```

Ldyn <- lapply(mycoplasma_strings_major_cluster, FUN = function(x) {dynchars(seq = strsplit(x, split = '')[[1]], average_interval_size = dyn_int, tss = 600-80, boundaries = c(160, 70), strand = 'forward')})
#E0
mat_dyn <- sapply(Ldyn, function(x) {x$E0})
means_to_mat <- c( apply((mat_dyn),  MARGIN = 1, FUN = median))
`Sequence (nts)` <- -rev(seq_along(means_to_mat))

matplot(`Sequence (nts)`, mat_dyn, type = 'l', lwd = 0.1,  col = 1, ylab = 'Open states activation energy (kcal/mole)', main = paste('Open states activation energy (sliding window', dyn_int, 'nts)\n for 370 non-vlhA promoter sequences'))
lines(`Sequence (nts)`, means_to_mat, col = 'magenta', lwd = 2)

gaa_repeats <- -155:-130
approx_tsss <- max(gaa_repeats)+50

segments(x0 = min(gaa_repeats), y0 = max(mat_dyn, na.rm = T), x1 = max(gaa_repeats), y1 = max(mat_dyn, na.rm = T), lwd = 3, col = 'darkred')
# #arrows(x0 = approx_tsss, y0 = 280, x1 = approx_tsss, y1 = 255, col = 'darkred', lwd = 3)
abline(v = approx_tsss, lty=2, lwd = 2)
legend('topright', legend = c('GAA-repeats', 'Approximate TSS', 'Median values'), lty = c(1,2,1), col = c('darkred', 'black', 'magenta'), lwd = 2)
