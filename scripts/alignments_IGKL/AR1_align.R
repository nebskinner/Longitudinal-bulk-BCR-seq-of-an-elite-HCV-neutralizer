.libPaths("/home/gdskinnerlab/nes002/software/R_packages")
  
## load packages
install.packages("utf8", repos = "http://cran.us.r-project.org")
library(utf8)

install.packages("dplyr", repos = "http://cran.us.r-project.org")
library(dplyr)

install.packages("stringi", repos = "http://cran.us.r-project.org")
library(stringi)

install.packages("ff", repos = "http://cran.us.r-project.org")
library(ff)

install.packages("RecordLinkage", repos = "http://cran.us.r-project.org")
library(RecordLinkage)

install.packages("RcppAlgos", repos = "http://cran.us.r-project.org")
library(RcppAlgos)

install.packages("withr", repos = "http://cran.us.r-project.org")
library(withr)

install.packages("doSNOW", repos = "http://cran.us.r-project.org")
library(doSNOW)

install.packages("foreach", repos = "http://cran.us.r-project.org")
library(foreach)

#### bnabs
## load sequences to be aligned
seq_df  <- read.table("/home/gdskinnerlab/nes002/serial_timepoints_C110/CDR3_alignment/CDR3_sequences_IGKL/AR1_bnab_IGKL.txt", sep = "\t", header=TRUE)

CDR3_sequences <- seq_df$CDR3_aa

## quantify sequence homologies
simil_quant <- data.frame()
  for(i in 1:length(CDR3_sequences)){
    Ldist <- levenshteinSim(CDR3_sequences[i], CDR3_sequences[-i])
    Ldist_df <- cbind.data.frame(CDR3.aa = CDR3_sequences[-i], Ldist, "Comp_seq" = rep(CDR3_sequences[i], length(Ldist)), "Mean_homol" =  rep(mean(Ldist), length(Ldist)), stringsAsFactors=FALSE)
    simil_quant <- rbind.data.frame(simil_quant, Ldist_df)
  }
  
  homol_order <- simil_quant %>%
    distinct(Comp_seq, Mean_homol) %>%
    mutate(AA_length = nchar(Comp_seq)) %>%
    group_by(AA_length) %>%
    arrange(., desc(Mean_homol), by_group = TRUE) %>%
    arrange(., desc(AA_length)) %>%
    ungroup() %>%
    mutate(blanks = max(AA_length)-AA_length)
  
  
## align sequences, optimizing homology
cluster = makeCluster(5, type = "SOCK")
registerDoSNOW(cluster)  

homol_seq <- data.frame(Seq = homol_order$Comp_seq[1], Avg_homol = NA)
  
homol_seq <-
  foreach(j=2:nrow(homol_order), .combine = "rbind", .packages = c("doSNOW", "dplyr", "ff", "RecordLinkage", "stringi", "RcppAlgos")) %dopar% {
    num_gaps <- homol_order$blanks[j]
    perm <- permuteGeneral(v = c(0:num_gaps), m= homol_order$AA_length[j], repetition = TRUE, constraintFun = "sum", comparisonFun = c(">","<"), limitConstraints = c((num_gaps)-1,(num_gaps)+1))
    
    gap_perm_pre <- c(rep(homol_order$Comp_seq[j], nrow(perm)))
    gap_perm <-
      foreach(k=1:nrow(perm), .combine="c") %do% {
      stri_sub_replace_all(gap_perm_pre[k], from = c(1:ncol(perm)), to = c(1:ncol(perm))-1, replacement = strrep("-",perm[k,]))
      }
    
    max_homol <-
      foreach(m=1:nrow(homol_seq), .combine = "rbind") %do% {
        leven_gap <- levenshteinSim(gap_perm, homol_seq$Seq[m])
        leven_max <- cbind.data.frame(Seq = gap_perm, Homol = leven_gap) 
      }
    
    max_avg_homol <- max_homol %>%
      group_by(Seq) %>%
      mutate(Avg_homol = mean(Homol)) %>%
      ungroup() %>%
      filter(Avg_homol == max(Avg_homol)) %>%
      distinct(Seq, .keep_all = TRUE) %>%
      select(Seq, Avg_homol)
    
    max_avg_homol <- max_avg_homol[1,]
  }

stopCluster(cluster)   
  
homol_seq <- rbind.data.frame(data.frame(Seq = homol_order$Comp_seq[1], Avg_homol = NA), homol_seq, stringsAsFactors = FALSE)
write.csv(homol_seq, "/home/gdskinnerlab/nes002/serial_timepoints_C110/CDR3_alignment/CDR3_alignments_IGKL/AR1_bnab_align_IGKL.csv")  

