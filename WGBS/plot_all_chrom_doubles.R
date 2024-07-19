#Delfi Dorussen
#Aim to produce figures showing the methylation profile of the met1 mutants 
#across all chromosomes, and add vertical lines to indicate partitions between
#the different chromosome regions

library(DMRcaller)

setwd("/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/")

chromosomes <- c("Chr1A","Chr1B","Chr1D","Chr2A","Chr2B","Chr2D","Chr3A","Chr3B","Chr3D",
                 "Chr4A","Chr4B","Chr4D","Chr5A","Chr5B","Chr5D","Chr6A","Chr6B","Chr6D",
                 "Chr7A","Chr7B","Chr7D")

chromosomes_lc <- c("chr1A","chr1B","chr1D","chr2A","chr2B","chr2D","chr3A","chr3B","chr3D",
                 "chr4A","chr4B","chr4D","chr5A","chr5B","chr5D","chr6A","chr6B","chr6D",
                 "chr7A","chr7B","chr7D")

lengths <- c(594102053, 689851867, 495453184, 780798554, 801256712, 651852608, 750843637, 830829763, 615552420,
             744588156, 673617493, 509857065, 709773740, 713149756, 566080676, 618079258, 720988477, 473592717,
             736706231, 750620384, 638686045)

R1_R2as <- c(59000000, 62000000, 44000000, 88000000, 59000000, 63000000, 62000000, 66000000, 49000000,
            41000000, 42000000, 10000000, 39000000, 52000000, 46000000, 60000000, 56000000, 98000000,
            89000000, 145000000, 84000000)

R2a_Cs <- c(151000000, 172000000, 81000000, 206000000, 248000000, 192000000, 249000000, 257000000, 167000000,
           180000000, 186000000, 135000000, 140000000, 140000000, 128000000, 216000000, 221000000, 164000000,
           239000000, 229000000, 266000000)

C_R2bs <- c(231000000, 277000000, 143000000, 367000000, 433000000, 364000000, 414000000, 407000000, 287000000,
           414000000, 360000000, 288000000, 260000000, 221000000, 207000000, 409000000, 429000000, 258000000,
           428000000, 418000000, 373000000)

R2b_R3s <- c(480000000, 534000000, 385000000, 662000000, 644000000, 520000000, 670000000, 728000000, 543000000,
            594000000, 537000000, 432000000, 427000000, 430000000, 345000000, 545000000, 651000000, 410000000,
            633000000, 578000000, 536000000)
  
  
for(i in 1:21){
  setwd("/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/")
  
  chromosome <- chromosomes[i]
  chromosome_lc <- chromosomes_lc[i]
  length <- lengths[i]
  R1_R2a <- R1_R2as[i]
  R2a_C <- R2a_Cs[i]
  C_R2b <- C_R2bs[i]
  R2b_R3 <- R2b_R3s[i]
  
  WT <- readBismark(paste("Unknown_BU703-001Z0001_1.",chromosome,".CX_report.txt", sep=""))
  aabbDD <- readBismark(paste("Unknown_BU703-001Z0005_1.",chromosome,".CX_report.txt", sep=""))
  aaBBdd <- readBismark(paste("Unknown_BU703-001Z0006_1.",chromosome,".CX_report.txt", sep=""))
  AAbbdd <- readBismark(paste("Unknown_BU703-001Z0007_1.",chromosome,".CX_report.txt", sep=""))
  
  chr <- GRanges(seqnames=chromosome_lc,ranges=IRanges(1,as.numeric(length)))
  
  chr_wt_CG <- computeMethylationProfile(WT, chr,
                                           windowSize=1000000, context="CG")
  chr_aabbDD_CG <- computeMethylationProfile(aabbDD, chr,
                                               windowSize=1000000, context="CG")
  chr_aaBBdd_CG <- computeMethylationProfile(aaBBdd, chr,
                                               windowSize=1000000, context="CG")
  chr_AAbbdd_CG <- computeMethylationProfile(AAbbdd, chr,
                                               windowSize=1000000, context="CG")
  
  chr_wt_CHG <- computeMethylationProfile(WT, chr,
                                            windowSize=1000000, context="CHG")
  chr_aabbDD_CHG <- computeMethylationProfile(aabbDD, chr,
                                                windowSize=1000000, context="CHG")
  chr_aaBBdd_CHG <- computeMethylationProfile(aaBBdd, chr,
                                                windowSize=1000000, context="CHG")
  chr_AAbbdd_CHG <- computeMethylationProfile(AAbbdd, chr,
                                                windowSize=1000000, context="CHG")
  
  chr_wt_CHH <- computeMethylationProfile(WT, chr,
                                            windowSize=1000000, context="CHH")
  chr_aabbDD_CHH <- computeMethylationProfile(aabbDD, chr,
                                                windowSize=1000000, context="CHH")
  chr_aaBBdd_CHH <- computeMethylationProfile(aaBBdd, chr,
                                                windowSize=1000000, context="CHH")
  chr_AAbbdd_CHH <- computeMethylationProfile(AAbbdd, chr,
                                                windowSize=1000000, context="CHH")
  
  setwd("/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/chr_plots_doubles/")
  
  pdf(paste("WT_all",chromosome,"_chr_CG.pdf", sep=""), width=8, height=6)
  plot((start(chr_wt_CG) + 500000), 100*chr_wt_CG$Proportion, type="l", lty=1,
       lwd=1, ylim=c(0,100), col="black", xlab="Position on Chromosome (bp)", ylab="Methylation Percentage")
  lines((start(chr_aaBBdd_CG) + 500000), 100*chr_aaBBdd_CG$Proportion, lty=1,
        lwd=1, col="#263392")
  lines((start(chr_AAbbdd_CG) + 500000), 100*chr_AAbbdd_CG$Proportion, lty=1,
        lwd=1, col="#4056f4")
  lines((start(chr_aabbDD_CG) + 500000), 100*chr_aabbDD_CG$Proportion, lty=1,
        lwd=1, col="#8c99f8")
  abline(v=as.numeric(R1_R2a), col="black")
  abline(v=as.numeric(R2a_C), col="black")
  abline(v=as.numeric(C_R2b), col="black")
  abline(v=as.numeric(R2b_R3), col="black")
  dev.off()
  
  pdf(paste("WT_all",chromosome,"_chr_CHG.pdf", sep=""), width=8, height=6)
  plot((start(chr_wt_CHG) + 500000), 100*chr_wt_CHG$Proportion, type="l", lty=1,
       lwd=1, ylim=c(0,100), col="black", xlab="Position on Chromosome (bp)", ylab="Methylation Percentage")
  lines((start(chr_aaBBdd_CHG) + 500000), 100*chr_aaBBdd_CHG$Proportion, lty=1,
        lwd=1, col="#263392")
  lines((start(chr_AAbbdd_CHG) + 500000), 100*chr_AAbbdd_CHG$Proportion, lty=1,
        lwd=1, col="#4056f4")
  lines((start(chr_aabbDD_CHG) + 500000), 100*chr_aabbDD_CHG$Proportion, lty=1,
        lwd=1, col="#8c99f8")
  abline(v=as.numeric(R1_R2a), col="black")
  abline(v=as.numeric(R2a_C), col="black")
  abline(v=as.numeric(C_R2b), col="black")
  abline(v=as.numeric(R2b_R3), col="black")
  dev.off()
  
  pdf(paste("WT_all",chromosome,"_chr_CHH.pdf", sep=""), width=8, height=6)
  plot((start(chr_wt_CHH) + 500000), 100*chr_wt_CHH$Proportion, type="l", lty=1,
       lwd=1, ylim=c(0,15), col="black", xlab="Position on Chromosome (bp)", ylab="Methylation Percentage")
  lines((start(chr_aaBBdd_CHH) + 500000), 100*chr_aaBBdd_CHH$Proportion, lty=1,
        lwd=1, col="#263392")
  lines((start(chr_AAbbdd_CHH) + 500000), 100*chr_AAbbdd_CHH$Proportion, lty=1,
        lwd=1, col="#4056f4")
  lines((start(chr_aabbDD_CHH) + 500000), 100*chr_aabbDD_CHH$Proportion, lty=1,
        lwd=1, col="#8c99f8")
  abline(v=as.numeric(R1_R2a), col="black")
  abline(v=as.numeric(R2a_C), col="black")
  abline(v=as.numeric(C_R2b), col="black")
  abline(v=as.numeric(R2b_R3), col="black")
  dev.off()
}
