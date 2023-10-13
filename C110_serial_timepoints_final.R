###### TABLE OF CONTENTS
  ##### PRIOR TO ANALYSIS
    #### Install necessary packages
    #### Establish novel functions used in the analysis
    #### Read in data
      ### Heavy chain
      ### Light chain
      ### mAb data
  ##### ANALYSIS
    #### CDR3 length
      ### CDRH3
      ### CDRL3
    #### V gene usage
      ### IGHV
      ### IGKV/IGLV
    #### Somatic hypermutation of the V gene
      ### IGHV
      ### IGKV/IGLV
  #### Conserved front-layer and non-front-layer bnab mutations compared to E2 negative clonotypes
    ### IGH
      ## CDR1/CDR2
          # front-layer
          # non-front-layer
      ## CDR3
        # front-layer
        # non-front-layer
    ### IGK/IGL
      ## CDR1/CDR2
        # front-layer
        # non-front-layer
      ## CDR3
        # front-layer
        # non-front-layer

####### PRIOR TO ANALYSIS 
#### Install necessary packages
library(tidyverse)
library(reshape2)
library(pheatmap)
library(readxl)
library(vegan)
library(ggseqlogo)
library(dunn.test)
library(ggpubr)
library(rstatix)
library(gtools)
library(ggpattern)
library(RColorBrewer)


###### Establish novel functions used in the analysis
### function to italicize text
make_italics <- function(x) {
  as.expression(lapply(x, function(y) bquote(italic(.(y)))))
}

### functions for mutation rate analyses
## function quantifying mutation rates
mutation_quant <- function(df1, df2, Vgene, mutation_type, mutation_subset) {
  # df1 and df2 are the dataframes of mutation usage to be compared
  # Vgene can be the name of a Vgene (e.g., "IGHV1-69") or "all" if all V genes are being analyzed
  # mutation_type is "position" (e.g., 36>), "subsitution" (e.g., N36>S), or "delta" (e.g., 36>S)
  # mutation_subset is "all" or a vector containing positions to analyze (e.g c(27:38) for CDR1)
 
  ##extract mutations relevant for specific V gene of interest ("Vgene") or for all V genes
  if(Vgene != "all"){                     ##filtering for a specific V gene
    df1_V <- filter(df1, V_gene == Vgene) %>%
      filter(V_gene_IMGT == Vgene) %>%
      mutate(V_mutation_pos = gsub("[A-Z]","", V_mutation)) %>%   #format mutations to reveal IMGT position only
      mutate(V_mutation_delta = sapply(V_mutation, function(x) paste(sub(".","", strsplit(x," ")[[1]]), collapse = " ")))  #format mutations to reveal IMGT position + mutated amino acid (but not germline amino acid)
    df2_V <- filter(df2, V_gene == Vgene) %>%
      filter(V_gene_IMGT == Vgene) %>%
      mutate(V_mutation_pos = gsub("[A-Z]","", V_mutation)) %>%
      mutate(V_mutation_delta = sapply(V_mutation, function(x) paste(sub(".","", strsplit(x," ")[[1]]), collapse = " ")))
  } else {                                ##if analyzing all V genes
    df1_V <- df1 %>%
      mutate(V_mutation_pos = gsub("[A-Z]","", V_mutation)) %>%
      mutate(V_mutation_delta = sapply(V_mutation, function(x) paste(sub(".","", strsplit(x," ")[[1]]), collapse = " ")))
    df2_V <- df2 %>%
      mutate(V_mutation_pos = gsub("[A-Z]","", V_mutation)) %>%
      mutate(V_mutation_delta = sapply(V_mutation, function(x) paste(sub(".","", strsplit(x," ")[[1]]), collapse = " ")))
  }
  
  ###extract positions intended for analysis
  if(mutation_subset[1] != "all"){      ##if analyzing specified IMGT positions
    filt_Vmut1 <- c()
    filt_Vmut_pos1 <- c()
    filt_Vmut_delta1 <- c()
    for(i in 1:nrow(df1_V)){      # filtering IMGT positions for each df of V mutations
      positions = paste(mutation_subset, collapse="|")
      
      V_muts <- unlist(strsplit(df1_V$V_mutation[i]," "))
      filt <- paste(V_muts[grepl(positions, V_muts)], collapse = " ")
      filt_Vmut1 <- c(filt_Vmut1, filt)
      
      V_muts_pos <- unlist(strsplit(df1_V$V_mutation_pos[i]," "))
      filt_pos <- paste(V_muts_pos[grepl(positions, V_muts_pos)], collapse = " ")
      filt_Vmut_pos1 <- c(filt_Vmut_pos1, filt_pos)
      
      V_muts_delta <- unlist(strsplit(df1_V$V_mutation_delta[i]," "))
      filt_delta <- paste(V_muts_delta[grepl(positions, V_muts_delta)], collapse = " ")
      filt_Vmut_delta1 <- c(filt_Vmut_delta1, filt_delta)
    }
    
    df1_V <- df1_V %>%          #inputing filtered positions into dfs
      mutate(V_mutation = filt_Vmut1) %>%
      mutate(V_mutation_pos = filt_Vmut_pos1) %>%
      mutate(V_mutation_delta = filt_Vmut_delta1)
    df1_V[df1_V==""] <- "none"
    
    
    # repeating the process for df2
    filt_Vmut2 <- c()
    filt_Vmut_pos2 <- c()
    filt_Vmut_delta2 <- c()
    for(i in 1:nrow(df2_V)){
      positions = paste(mutation_subset, collapse="|")
      V_muts <- unlist(strsplit(df2_V$V_mutation[i]," "))
      filt <- paste(V_muts[grepl(positions, V_muts)], collapse = " ")
      filt_Vmut2 <- c(filt_Vmut2, filt)
      
      V_muts_pos <- unlist(strsplit(df2_V$V_mutation_pos[i]," "))
      filt_pos <- paste(V_muts_pos[grepl(positions, V_muts_pos)], collapse = " ")
      filt_Vmut_pos2 <- c(filt_Vmut_pos2, filt_pos)
      
      V_muts_delta <- unlist(strsplit(df2_V$V_mutation_delta[i]," "))
      filt_delta <- paste(V_muts_delta[grepl(positions, V_muts_delta)], collapse = " ")
      filt_Vmut_delta2 <- c(filt_Vmut_delta2, filt_delta)
    }
    df2_V <- df2_V %>%
      mutate(V_mutation = filt_Vmut2) %>%
      mutate(V_mutation_pos = filt_Vmut_pos2) %>%
      mutate(V_mutation_delta = filt_Vmut_delta2)
    df2_V[df2_V==""] <- "none"
  } else{}              ## no filtering if analyzing all positions
  
  ## generate list of mutations (modified by whether analyzing substitution, position, or just the delta substitution)
  df1_subs <- unique(unlist(str_split(df1_V$V_mutation, " "))) 
  df1_pos <- unique(gsub("[A-Z]","", df1_subs))
  df1_delta <- unique(sub(".","", df1_subs))
  
  df2_subs <- unique(unlist(str_split(df2_V$V_mutation, " "))) 
  df2_pos <- unique(gsub("[A-Z]","", df2_subs))
  df2_delta <- unique(sub(".","", df2_subs))
  
  ## quantify df1 and df2 mutations
  if(mutation_type == "substitution"){      # if interested in substitutions
    mut_list1 = df1_subs
    df_mut_list1 = df1_V$V_mutation
    mut_list2 = df2_subs
    df_mut_list2 = df2_V$V_mutation
  } else{if(mutation_type == "position"){   # if interested in position
    mut_list1 = df1_pos
    df_mut_list1 = df1_V$V_mutation_pos
    mut_list2 = df2_pos
    df_mut_list2 = df2_V$V_mutation_pos
  } else{if(mutation_type == "delta"){     # if interested in delta substitution
    mut_list1 = df1_delta
    df_mut_list1 = df1_V$V_mutation_delta
    mut_list2 = df2_delta
    df_mut_list2 = df2_V$V_mutation_delta
  } else{print("ERROR: You must enter an accepted mutation type: substitution (e.g., A36>G), position (e.g., 36>), or delta (e.g., 36>G)")}
  }}
  
  # generate df with mutation counts and proportions for df1 and df2
  df1_pos_counts <- data.frame()
  for(i in 1:length(mut_list1)){
    subs_count <- sum(str_count(df_mut_list1, mut_list1[i]))
    subs_prop <- subs_count/nrow(df1_V)
    df_count <- cbind.data.frame(Substitution = mut_list1[i], Count = subs_count, Total = nrow(df1_V), Prop = subs_prop, Group = "df1")
    df1_pos_counts <- rbind.data.frame(df1_pos_counts, df_count)
  }
  
  df2_pos_counts <- data.frame()
  for(i in 1:length(mut_list2)){
    subs_count <- sum(str_count(df_mut_list2, mut_list2[i]))
    subs_prop <- subs_count/nrow(df2_V)
    df_count <- cbind.data.frame(Substitution = mut_list2[i], Count = subs_count, Total = nrow(df2_V), Prop = subs_prop, Group = "df2")
    df2_pos_counts <- rbind.data.frame(df2_pos_counts, df_count)
  }
  
  # combine df1 and df2 into one df for plotting
  df_pos_counts <- rbind.data.frame(df1_pos_counts, df2_pos_counts) %>%
    complete(Substitution, Group, fill = list(Count = 0)) %>%
    select(Substitution, Group, Count, Prop) %>%
    replace(is.na(.), 0) %>%
    reshape2::dcast(., Group ~ Substitution, value.var = "Prop")
  row.names(df_pos_counts) <- c("df1", "df2")
  
  # format df for pheatmap
  df_pos_hm <- t(df_pos_counts[c("df1", "df2"),-1])
  
  return(list(df_pos_hm, df1_pos_counts, df2_pos_counts))   ## generate a list with the df structured for pheatmap [[1]], df1 mutation counts/props [[2]], df2 mutation counts/props [[3]] 
}

## function for statistical comparison of mutation usage; it relies on the mutation_quant function and therefore requires the same inputs (df1, df2, Vgene, mutation_type, mutation_subset)
stat_mut <- function(df1, df2, Vgene, mutation_type, mutation_subset){
  # generate one dataframe from df1 df2
  df_stat <- rbind.data.frame(mutation_quant(df1, df2, "all", "delta", cdr1cdr2)[[2]], 
                              mutation_quant(df1, df2, "all", "delta", cdr1cdr2)[[3]]) %>%
    complete(Substitution, Group, fill = list(Count = 0)) %>%   # for each mutation that is only present in one group (df1 or df2), add a row for the other group with count = 0
  group_by(Group) %>%
  fill(Total, .direction = "updown") %>%
  replace(is.na(.), 0)

## obtain uncorrected P values using Fisher test to compare proportion of mutation usage (# times used/total # cloonotypes) between the 2 groups
  df_pvalue <- data.frame()
  for(i in 1:length(unique(df_stat$Substitution))){
    subst <- filter(df_stat, Substitution == unique(df_stat$Substitution)[i])
    fisherP <- fisher.test(matrix(c(subst$Count[1], subst$Count[2], subst$Total[1]-subst$Count[1], subst$Total[2]-subst$Count[2]), ncol=2))$p.value
    subs_df <- cbind.data.frame(Substitution = unique(df_stat$Substitution)[i], 
                                df1_df2 = fisherP)
    df_pvalue  <- rbind.data.frame(df_pvalue, subs_df, stringsAsFactors = FALSE)
  }
  
  ## correct P values using Benjamini-Hochberg correction (aka FDR)
  df_pvalue_corr <- cbind.data.frame(df_pvalue, df_pvalue_corr = p.adjust(df_pvalue$df1_df2, method = "fdr")) 
  
  ## extract mutations with significant P values (p < 0.01 is used; for ~ 20 comparisons, this would mean 0.2 false positives; for 100 comparisons, it would mean 1 false positive) 
  df_sig_pvalue <- df_pvalue_corr %>%
    filter(df_pvalue_corr < 0.01) %>%
    select(Substitution, df_pvalue_corr)
  
  # generate a list with the df containing mutation counts/proportions [[1]], P values and corrected P values for each mutation [[2]], and just the mutations with significant P values [[3]] 
  return(list(df_stat, df_pvalue_corr, df_sig_pvalue))    
}


### function to find amino acid usage at each position of a sequence:
AA_usage <- function (AA_vec, AA_length) {
  AA_use_df <- data.frame()
  for(i in 1:AA_length){
    AAs <- paste0(substring(AA_vec,i,i), collapse="")
    A_count <- str_count(AAs, "A")
    R_count <- str_count(AAs, "R")
    N_count <- str_count(AAs, "N")
    D_count <- str_count(AAs, "D")
    C_count <- str_count(AAs, "C")
    E_count <- str_count(AAs, "E")
    Q_count <- str_count(AAs, "Q")
    G_count <- str_count(AAs, "G")
    H_count <- str_count(AAs, "H")
    I_count <- str_count(AAs, "I")
    L_count <- str_count(AAs, "L")
    K_count <- str_count(AAs, "K")
    M_count <- str_count(AAs, "M")
    F_count <- str_count(AAs, "F")
    P_count <- str_count(AAs, "P")
    S_count <- str_count(AAs, "S")
    T_count <- str_count(AAs, "T")
    W_count <- str_count(AAs, "W")
    Y_count <- str_count(AAs, "Y")
    V_count <- str_count(AAs, "V")
    use_freq<- c(A_count, R_count, N_count, D_count, C_count, E_count, Q_count, G_count, H_count, I_count, L_count, K_count, M_count, F_count, P_count, S_count, T_count, W_count, Y_count, V_count)
    AA_use_df <- rbind.data.frame(AA_use_df, use_freq)
  }
  
  AA_freq <- t(AA_use_df)
  rownames(AA_freq) <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  colnames(AA_freq) <- paste0("Pos", 1:AA_length)
  
  return(AA_freq)
}


####### Read in data
### Heavy chain
## E2 positive for each time point
C110_d279_pos_IGH <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d279_E2pos_IGH_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d279")  %>%
  mutate(Group = "E2pos")

C110_d540_pos_IGH <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d540_E2pos_IGH_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d540")  %>%
  mutate(Group = "E2pos")

C110_d1267_pos_IGH <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d1267_E2pos_IGH_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d1267")  %>%
  mutate(Group = "E2pos")

C110_d1842_pos_IGH <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d1842_E2pos_IGH_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d1842")  %>%
  mutate(Group = "E2pos") 

## E2 negative for each time point
C110_d279_neg_IGH <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d279_E2neg_IGH_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d279")  %>%
  mutate(Group = "E2neg")

C110_d540_neg_IGH <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d540_E2neg_IGH_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d540")  %>%
  mutate(Group = "E2neg")

C110_d1267_neg_IGH <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d1267_E2neg_IGH_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d1267") %>%
  mutate(Group = "E2neg")

C110_d1842_neg_IGH <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d1842_E2neg_IGH_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d1842")  %>%
  mutate(Group = "E2neg")

## combine all heavy chain data
C110_IGH <- rbind.data.frame(C110_d279_neg_IGH, C110_d279_pos_IGH, C110_d540_neg_IGH, C110_d540_pos_IGH, C110_d1267_neg_IGH, C110_d1267_pos_IGH, C110_d1842_neg_IGH, C110_d1842_pos_IGH, stringsAsFactors = FALSE) %>%
  replace(is.na(.), 0) %>%
  filter(Frame.Shift != "Frame Shift") %>%
  filter(Stop.Codon != "Stop Codon") %>%
  select(Timepoint, Group, Sequence = Clonal.Sequence, CDR3_nt = CDR3.Sequence, CDR3_aa = CDR3.Amino.Acid.Sequence, CDR3_len = Amino.Acid.Length, V_gene = V.segment, D_gene = D.segment, J_gene = J.segment, C_gene = C.segment) %>%
  filter(CDR3_aa != "")

### Light chain
## Kappa
# E2 positive for each time point
C110_d279_pos_IGK <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d279_E2pos_IGK_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d279")  %>%
  mutate(Group = "E2pos")

C110_d540_pos_IGK <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d540_E2pos_IGK_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d540")  %>%
  mutate(Group = "E2pos")

C110_d1267_pos_IGK <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d1267_E2pos_IGK_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d1267")  %>%
  mutate(Group = "E2pos")

C110_d1842_pos_IGK <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d1842_E2pos_IGK_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d1842")  %>%
  mutate(Group = "E2pos") 

# E2 negative for each time point
C110_d279_neg_IGK <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d279_E2neg_IGK_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d279")  %>%
  mutate(Group = "E2neg")

C110_d540_neg_IGK <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d540_E2neg_IGK_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d540")  %>%
  mutate(Group = "E2neg")

C110_d1267_neg_IGK <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d1267_E2neg_IGK_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d1267") %>%
  mutate(Group = "E2neg")

C110_d1842_neg_IGK <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d1842_E2neg_IGK_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d1842")  %>%
  mutate(Group = "E2neg")

## Lambda
# E2 positive for each time point
C110_d279_pos_IGL <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d279_E2pos_IGL_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d279")  %>%
  mutate(Group = "E2pos")

C110_d540_pos_IGL <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d540_E2pos_IGL_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d540")  %>%
  mutate(Group = "E2pos")

C110_d1267_pos_IGL <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d1267_E2pos_IGL_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d1267")  %>%
  mutate(Group = "E2pos")

C110_d1842_pos_IGL <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d1842_E2pos_IGL_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d1842")  %>%
  mutate(Group = "E2pos") 

# E2 negative for each time point
C110_d279_neg_IGL <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d279_E2neg_IGL_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d279")  %>%
  mutate(Group = "E2neg")

C110_d540_neg_IGL <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d540_E2neg_IGL_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d540")  %>%
  mutate(Group = "E2neg")

C110_d1267_neg_IGL <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d1267_E2neg_IGL_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d1267") %>%
  mutate(Group = "E2neg")

C110_d1842_neg_IGL <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Raw_data/C110_ImmuneProfiler_output/C110_d1842_E2neg_IGL_mig_fl_clones_result.csv") %>%
  mutate(Timepoint = "d1842")  %>%
  mutate(Group = "E2neg")

## combine all light chain data
C110_IGKL <- rbind.data.frame(C110_d279_neg_IGK, C110_d279_pos_IGK, C110_d540_neg_IGK, C110_d540_pos_IGK, C110_d1267_neg_IGK, C110_d1267_pos_IGK, C110_d1842_neg_IGK, C110_d1842_pos_IGK, C110_d279_neg_IGL, C110_d279_pos_IGL, C110_d540_neg_IGL, C110_d540_pos_IGL, C110_d1267_neg_IGL, C110_d1267_pos_IGL, C110_d1842_neg_IGL, C110_d1842_pos_IGL, stringsAsFactors = FALSE) %>%
  filter(!(grepl("IGHV", V.segment))) %>%
  replace(is.na(.), 0) %>%
  filter(Frame.Shift != "Frame Shift") %>%
  filter(Stop.Codon != "Stop Codon") %>%
  select(Timepoint, Group, Sequence = Clonal.Sequence, CDR3_nt = CDR3.Sequence, CDR3_aa = CDR3.Amino.Acid.Sequence, CDR3_len = Amino.Acid.Length, V_gene = V.segment, D_gene = D.segment, J_gene = J.segment, C_gene = C.segment) %>%
  filter(CDR3_aa != "")


### Somatic hypermutation data from IMGT (heavy and light chain combined)
C110_IMGT_all <- read_excel("/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Sequences/C110_IMGT_all_prod.xlsx") %>%
  mutate(Sequence = toupper(Sequence)) 

##format VJ genes 
C110_IMGT_all$V_gene_IMGT <- gsub("(\\*+[0-9]+\\s+[A-Z]*)|(\\*+[0-9]+\\s\\()","", C110_IMGT_all$V_gene_IMGT)
C110_IMGT_all$V_gene_IMGT <- gsub("[A-Z]+\\)","",C110_IMGT_all$V_gene_IMGT)
C110_IMGT_all$V_gene_IMGT <- gsub("or", "", C110_IMGT_all$V_gene_IMGT)
C110_IMGT_all$V_gene_IMGT <- gsub("\\,", "", C110_IMGT_all$V_gene_IMGT)
C110_IMGT_all$V_gene_IMGT <- gsub("D", "", C110_IMGT_all$V_gene_IMGT)
C110_IMGT_all$V_gene_IMGT <- gsub("\\[.|\\]", "", C110_IMGT_all$V_gene_IMGT)

C110_IMGT_all$V_gene_IMGT <- sub(" .*", "", C110_IMGT_all$V_gene_IMGT)

C110_IMGT_all$J_gene_IMGT <- gsub("\\*+.*","", C110_IMGT_all$J_gene_IMGT)

## change format of mutations
C110_IMGT_all$V_mutation <- sapply(C110_IMGT_all$V_mutation, function(x) paste(unique(str_extract_all(x,"[A-Z]+[0-9]+>+[A-Z]")[[1]]), collapse = " "))
C110_IMGT_all$V_mutation <- as.vector(C110_IMGT_all$V_mutation)

C110_IMGT_all$CDR1_mutation <- sapply(C110_IMGT_all$CDR1_mutation, function(x) paste(unique(str_extract_all(x,"[A-Z]+[0-9]+>+[A-Z]")[[1]]), collapse = " "))
C110_IMGT_all$CDR1_mutation <- as.vector(C110_IMGT_all$CDR1_mutation)

C110_IMGT_all$CDR2_mutation <- sapply(C110_IMGT_all$CDR2_mutation, function(x) paste(unique(str_extract_all(x,"[A-Z]+[0-9]+>+[A-Z]")[[1]]), collapse = " "))
C110_IMGT_all$CDR2_mutation <- as.vector(C110_IMGT_all$CDR2_mutation)

## replace blank entries
C110_IMGT_all[C110_IMGT_all=="NA"] <- "none"
C110_IMGT_all[C110_IMGT_all==""] <- "none"

## Separate IGH and IGKL
C110_IMGT_IGH <- C110_IMGT_all %>%
  filter(Chain == "IGH")

C110_IMGT_IGKL <- C110_IMGT_all %>%
  filter(Chain != "IGH")

### mAb data
## Heavy chain
mAbs_IGH <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Sequences/mAbs_C110_IGH.csv") 

# change format of mutations
mAbs_IGH$V_mutation <- sapply(mAbs_IGH$V_mutation, function(x) paste(unique(str_extract_all(x,"[A-Z]+[0-9]+>+[A-Z]")[[1]]), collapse = " "))
mAbs_IGH$V_mutation <- as.vector(mAbs_IGH$V_mutation)

mAbs_IGH$CDR1_mutation <- sapply(mAbs_IGH$CDR1_mutation, function(x) paste(unique(str_extract_all(x,"[A-Z]+[0-9]+>+[A-Z]")[[1]]), collapse = " "))
mAbs_IGH$CDR1_mutation <- as.vector(mAbs_IGH$CDR1_mutation)

mAbs_IGH$CDR2_mutation <- sapply(mAbs_IGH$CDR2_mutation, function(x) paste(unique(str_extract_all(x,"[A-Z]+[0-9]+>+[A-Z]")[[1]]), collapse = " "))
mAbs_IGH$CDR2_mutation <- as.vector(mAbs_IGH$CDR2_mutation)

# replace blank entries
mAbs_IGH[mAbs_IGH=="NA"] <- "none"
mAbs_IGH[mAbs_IGH==""] <- "none"


## Light chain
mAbs_IGKL <- read.csv(file = "/Volumes/nes002/Projects_Ongoing/Serial_Timepoints_C110/Sequences/mAbs_c110_IGKL.csv") 

# change format of mutations
mAbs_IGKL$V_mutation <- sapply(mAbs_IGKL$V_mutation, function(x) paste(unique(str_extract_all(x,"[A-Z]+[0-9]+>+[A-Z]")[[1]]), collapse = " "))
mAbs_IGKL$V_mutation <- as.vector(mAbs_IGKL$V_mutation)

mAbs_IGKL$CDR1_mutation <- sapply(mAbs_IGKL$CDR1_mutation, function(x) paste(unique(str_extract_all(x,"[A-Z]+[0-9]+>+[A-Z]")[[1]]), collapse = " "))
mAbs_IGKL$CDR1_mutation <- as.vector(mAbs_IGKL$CDR1_mutation)

mAbs_IGKL$CDR2_mutation <- sapply(mAbs_IGKL$CDR2_mutation, function(x) paste(unique(str_extract_all(x,"[A-Z]+[0-9]+>+[A-Z]")[[1]]), collapse = " "))
mAbs_IGKL$CDR2_mutation <- as.vector(mAbs_IGKL$CDR2_mutation)

# replace blank entries
mAbs_IGKL[mAbs_IGKL=="NA"] <- "none"
mAbs_IGKL[mAbs_IGKL==""] <- "none"



####### ANALYSIS
###### CDR3 length
#### CDRH3
### make dataframe
C110_CDR3_len_IGH <- rbind.data.frame(
  cbind.data.frame(SeqID = NA, Timepoint = C110_IGH$Timepoint, Group = C110_IGH$Group, CDR3_aa = C110_IGH$CDR3_aa),
  cbind.data.frame(SeqID = mAbs_IGH$SeqID, Timepoint = rep("mAb", nrow(mAbs_IGH)), Group=rep("E2pos", nrow(mAbs_IGH)), CDR3_aa = mAbs_IGH$CDR3_aa),
  stringsAsFactors = FALSE
) %>%
  mutate(CDR3_len = nchar(CDR3_aa)) %>%
  mutate(time_group = ifelse(Group == "E2pos", Timepoint, "E2neg"))

### statistics
## normality
shapiro.test(C110_CDR3_len_IGH[C110_CDR3_len_IGH$time_group == "d279",]$CDR3_len) # normal
shapiro.test(C110_CDR3_len_IGH[C110_CDR3_len_IGH$time_group == "d540",]$CDR3_len) # not normal
shapiro.test(C110_CDR3_len_IGH[C110_CDR3_len_IGH$time_group == "d1267",]$CDR3_len) # not normal
shapiro.test(C110_CDR3_len_IGH[C110_CDR3_len_IGH$time_group == "d1842",]$CDR3_len)  # normal
shapiro.test(C110_CDR3_len_IGH[C110_CDR3_len_IGH$time_group == "mAb",]$CDR3_len) # not normal

## hypothesis testing
C110_CDR3_len_IGH$time_group <- factor(C110_CDR3_len_IGH$time_group, levels = c("E2neg", "d279", "d540", "d1267", "d1842", "mAb"))

kruskal.test(CDR3_len ~ time_group, data = C110_CDR3_len_IGH)
CDR3_len_IGH_stat <- dunn_test(C110_CDR3_len_IGH, CDR3_len ~ time_group, p.adjust.method = "BH", detailed = FALSE)

### plot
med_IGH_cdr3 <- data.frame(x_pos=c(1, 2, 3, 4, 5, 6), 
                       y_pos=c(median(C110_CDR3_len_IGH[C110_CDR3_len_IGH$time_group == "E2neg",]$CDR3_len),
                               median(C110_CDR3_len_IGH[C110_CDR3_len_IGH$time_group== "d279",]$CDR3_len),
                               median(C110_CDR3_len_IGH[C110_CDR3_len_IGH$time_group == "d540",]$CDR3_len),
                               median(C110_CDR3_len_IGH[C110_CDR3_len_IGH$time_group == "d1267",]$CDR3_len),
                               median(C110_CDR3_len_IGH[C110_CDR3_len_IGH$time_group== "d1842",]$CDR3_len),
                               median(C110_CDR3_len_IGH[C110_CDR3_len_IGH$time_group == "mAb",]$CDR3_len)),
                       x_add = c(0.12, 0.12, 0.12, 0.12, 0.12, 0.12))

C110_CDR3_len_IGH_p <- ggplot(C110_CDR3_len_IGH, aes(x=time_group, y=CDR3_len)) +
  geom_point(aes(fill=time_group, color=time_group), shape =1, position = position_jitter(width=0.175), size=1, stroke=0.95)+
  geom_violin(color= "#808482", show.legend = FALSE, lwd=0.6, alpha=0) +
  geom_segment(data = med_IGH_cdr3 , aes(x = x_pos-x_add, xend = x_pos+x_add, y = y_pos, yend = y_pos), color = "#676767", linewidth=1.5)+
  scale_shape_manual(values=c(1, 1, 1, 1, 1, 1))+
  scale_color_manual(values=c("#BEBEBE", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#009E73")) +
  scale_x_discrete(breaks=c("E2neg", "d279", "d540", "d1267", "d1842", "mAb"),
                   labels=c("E2-", "d279", "d540", "d1267", "d1842", "mAb"))+
  stat_pvalue_manual(CDR3_len_IGH_stat, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(40.5), tip.length = 0.005, vjust= 0.7, bracket.size=0.45, size=5.5) +
  theme(axis.text.y=element_text(color="black", size=16),
        axis.text.x=element_text(color="black", size=16),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.55, color = "black"),
        axis.title=element_text(size=18),
        legend.position = "none",
        aspect.ratio = 1)+
  labs(x= "",
       y= "CDRH3 Amino Acid Length")
C110_CDR3_len_IGH_p 


#### CDRL3
C110_CDR3_len_IGKL <- rbind.data.frame(
  cbind.data.frame(SeqID = NA, Timepoint = C110_IGKL$Timepoint, Group = C110_IGKL$Group, CDR3_aa = C110_IGKL$CDR3_aa),
  cbind.data.frame(SeqID = mAbs_IGKL$SeqID, Timepoint = rep("mAb", nrow(mAbs_IGKL)), Group=rep("E2pos", nrow(mAbs_IGKL)), CDR3_aa = mAbs_IGKL$CDR3_aa),
  stringsAsFactors = FALSE
) %>%
  mutate(CDR3_len = nchar(CDR3_aa)) %>%
  mutate(time_group = ifelse(Group == "E2pos", Timepoint, "E2neg")) 

### statistics
## normality
shapiro.test(C110_CDR3_len_IGKL[C110_CDR3_len_IGKL$time_group == "d279",]$CDR3_len) # not normal
shapiro.test(C110_CDR3_len_IGKL[C110_CDR3_len_IGKL$time_group == "d540",]$CDR3_len) # not normal
shapiro.test(C110_CDR3_len_IGKL[C110_CDR3_len_IGKL$time_group == "d1267",]$CDR3_len) # not normal
shapiro.test(C110_CDR3_len_IGKL[C110_CDR3_len_IGKL$time_group == "d1842",]$CDR3_len)  # not normal
shapiro.test(C110_CDR3_len_IGKL[C110_CDR3_len_IGKL$time_group == "mAb",]$CDR3_len) # not normal

C110_CDR3_len_IGKL$time_group <- factor(C110_CDR3_len_IGKL$time_group, levels = c("E2neg", "d279", "d540", "d1267", "d1842", "mAb"))
kruskal.test(CDR3_len ~ time_group, data = C110_CDR3_len_IGKL) # not significant

### plot
med_IGKL_cdr3 <- data.frame(x_pos=c(1, 2, 3, 4, 5, 6), 
                       y_pos=c(median(C110_CDR3_len_IGKL[C110_CDR3_len_IGKL$time_group == "E2neg",]$CDR3_len),
                               median(C110_CDR3_len_IGKL[C110_CDR3_len_IGKL$time_group== "d279",]$CDR3_len),
                               median(C110_CDR3_len_IGKL[C110_CDR3_len_IGKL$time_group == "d540",]$CDR3_len),
                               median(C110_CDR3_len_IGKL[C110_CDR3_len_IGKL$time_group == "d1267",]$CDR3_len),
                               median(C110_CDR3_len_IGKL[C110_CDR3_len_IGKL$time_group== "d1842",]$CDR3_len),
                               median(C110_CDR3_len_IGKL[C110_CDR3_len_IGKL$time_group == "mAb",]$CDR3_len)),
                       x_add = c(0.12, 0.12, 0.12, 0.12, 0.12, 0.12))

C110_CDR3_len_IGKL_p <- ggplot(C110_CDR3_len_IGKL, aes(x=time_group, y=CDR3_len)) +
  geom_point(aes(fill=time_group, color=time_group), shape =1, position = position_jitter(width=0.175), size=1, stroke=0.95)+
  geom_violin(color= "#808482", show.legend = FALSE, lwd=0.6, alpha=0) +
  geom_segment(data = med_IGKL_cdr3 , aes(x = x_pos-x_add, xend = x_pos+x_add, y = y_pos, yend = y_pos), color = "#676767", linewidth=1.5)+
  scale_shape_manual(values=c(1, 1, 1, 1, 1, 1))+
  scale_color_manual(values=c("#BEBEBE", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#009E73")) +
  scale_x_discrete(breaks=c("E2neg", "d279", "d540", "d1267", "d1842", "mAb"),
                   labels=c("E2-", "d279", "d540", "d1267", "d1842", "mAb"))+
  theme(axis.text.y=element_text(color="black", size=16),
        axis.text.x=element_text(color="black", size=16),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.55, color = "black"),
        axis.title=element_text(size=18),
        legend.position = "none",
        aspect.ratio = 1)+
  labs(x= "",
       y= "CDRL3 Amino Acid Length")
C110_CDR3_len_IGKL_p 


###### V gene usage 
#### IGHV
### make dataframe
C110_IGHVJ <- rbind.data.frame(
  cbind.data.frame(SeqID = NA, Timepoint = C110_IGH$Timepoint, Group = C110_IGH$Group, CDR3 = C110_IGH$CDR3_aa, V_gene = C110_IGH$V_gene, J_gene = C110_IGH$J_gene),
  cbind.data.frame(SeqID = mAbs_IGH$SeqID, Timepoint = rep("mAb", nrow(mAbs_IGH)), Group=rep("E2pos", nrow(mAbs_IGH)), CDR3 = mAbs_IGH$CDR3_aa, V_gene = mAbs_IGH$V_gene, J_gene = mAbs_IGH$J_gene),
  stringsAsFactors = FALSE
)

C110_Vgene_IGH <- C110_IGHVJ %>%
  mutate(time_group = ifelse(Group == "E2pos", Timepoint, "E2neg")) %>%   #combine all E2 neg clonotypes into 1 group
  group_by(time_group, V_gene) %>%
  mutate(V_count = n()) %>%
  distinct(V_gene, time_group, V_count) %>%
  ungroup() %>%
  complete(V_gene, time_group, fill = list(V_count =0)) %>%
  group_by(time_group) %>%
  mutate(V_total = sum(V_count)) %>%
  ungroup() %>%
  mutate(V_prop = V_count/V_total)
  
### statistics
C110_IGHV_list <- unique(C110_Vgene_IGH$V_gene)

C110_IGHV_stat_pre <- data.frame()
for(i in 1:length(C110_IGHV_list)){
  Vgene_df <- filter(C110_Vgene_IGH, V_gene == C110_IGHV_list[i]) %>%
    mutate(nonV_count = V_total-V_count)
  df_d279 <- Vgene_df %>%
    filter(time_group == "d279")
  df_d540 <- Vgene_df %>%
    filter(time_group == "d540")
  df_d1267 <- Vgene_df %>%
    filter(time_group == "d1267")
  df_d1842 <- Vgene_df %>%
    filter(time_group == "d1842")
  df_mAb <- Vgene_df %>%
    filter(time_group == "mAb")
  df_E2n <- Vgene_df %>%
    filter(time_group == "E2neg")
  
  d279_E2n <- fisher.test(matrix(c(df_d279$V_count, df_E2n$V_count, df_d279$nonV_count, df_E2n$nonV_count), ncol=2))$p.value
  d540_E2n <- fisher.test(matrix(c(df_d540$V_count, df_E2n$V_count, df_d540$nonV_count, df_E2n$nonV_count), ncol=2))$p.value
  d1267_E2n <- fisher.test(matrix(c(df_d1267$V_count, df_E2n$V_count, df_d1267$nonV_count, df_E2n$nonV_count), ncol=2))$p.value
  d1842_E2n <- fisher.test(matrix(c(df_d1842$V_count, df_E2n$V_count, df_d1842$nonV_count, df_E2n$nonV_count), ncol=2))$p.value
  mAb_E2n <- fisher.test(matrix(c(df_mAb$V_count, df_E2n$V_count, df_mAb$nonV_count, df_E2n$nonV_count), ncol=2))$p.value
  
  d279_d540 <- fisher.test(matrix(c(df_d279$V_count, df_d540$V_count, df_d279$nonV_count, df_d540$nonV_count), ncol=2))$p.value
  d279_d1267 <- fisher.test(matrix(c(df_d279$V_count, df_d1267$V_count, df_d279$nonV_count, df_d1267$nonV_count), ncol=2))$p.value
  d279_d1842 <- fisher.test(matrix(c(df_d279$V_count, df_d1842$V_count, df_d279$nonV_count, df_d1842$nonV_count), ncol=2))$p.value
  d279_mAb <- fisher.test(matrix(c(df_d279$V_count, df_mAb$V_count, df_d279$nonV_count, df_mAb$nonV_count), ncol=2))$p.value
  
  d540_d1267 <- fisher.test(matrix(c(df_d540$V_count, df_d1267$V_count, df_d540$nonV_count, df_d1267$nonV_count), ncol=2))$p.value
  d540_d1842 <- fisher.test(matrix(c(df_d540$V_count, df_d1842$V_count, df_d540$nonV_count, df_d1842$nonV_count), ncol=2))$p.value
  d540_mAb <- fisher.test(matrix(c(df_d540$V_count, df_mAb$V_count, df_d540$nonV_count, df_mAb$nonV_count), ncol=2))$p.value
  
  d1267_d1842 <- fisher.test(matrix(c(df_d1267$V_count, df_d1842$V_count, df_d1267$nonV_count, df_d1842$nonV_count), ncol=2))$p.value
  d1267_mAb <- fisher.test(matrix(c(df_d1267$V_count, df_mAb$V_count, df_d1267$nonV_count, df_mAb$nonV_count), ncol=2))$p.value
  
  d1842_mAb <- fisher.test(matrix(c(df_d1842$V_count, df_mAb$V_count, df_d1842$nonV_count, df_mAb$nonV_count), ncol=2))$p.value
  
  stat_df <- cbind.data.frame(V_gene = C110_IGHV_list[i], d279_E2n = d279_E2n, d540_E2n = d540_E2n, d1267_E2n = d1267_E2n, d1842_E2n = d1842_E2n, mAb_E2n = mAb_E2n, d279_d540 = d279_d540, d279_d1267 = d279_d1267, d279_d1842 = d279_d1842, d279_mAb = d279_mAb, d540_d1267 = d540_d1267, d540_d1842 = d540_d1842, d540_mAb = d540_mAb, d1267_d1842 = d1267_d1842, d1267_mAb = d1267_mAb, d1842_mAb = d1842_mAb)
  C110_IGHV_stat_pre <- rbind.data.frame(C110_IGHV_stat_pre, stat_df, stringsAsFactors = FALSE)
}

C110_IGHV_stat <- as.data.frame(matrix(p.adjust(as.vector(as.matrix(C110_IGHV_stat_pre[,-1])), method="fdr"),ncol= ncol(C110_IGHV_stat_pre)-1))
rownames(C110_IGHV_stat) <- C110_IGHV_stat_pre[,1]
colnames(C110_IGHV_stat) <- colnames(C110_IGHV_stat_pre[,-1])

C110_IGHV_stat_sig <- C110_IGHV_stat %>%
  filter(if_any(d279_E2n:d1842_mAb, ~ .x < 0.05))

### heatmaps
## all genes
C110_Vgene_IGH_hm_pre <- C110_Vgene_IGH %>%
  select(time_group, V_gene, V_prop) %>%
  spread(key=V_gene, value = V_prop) %>%
  arrange(factor(time_group, levels = c("E2neg", "d279", "d540", "d1267", "d1842", "mAb")))

Vgene_ord_IGH <- data.frame(Vgene = colnames(C110_Vgene_IGH_hm_pre[-1]), 
                            Vgene_edit = as.character(gsub("[A-Z]", "", colnames(C110_Vgene_IGH_hm_pre[-1]))),
                             group = sub("-.*", "", colnames(C110_Vgene_IGH_hm_pre)[-1]),
                             num = sub(".*?-", "", colnames(C110_Vgene_IGH_hm_pre)[-1])) %>%
  mutate(num = as.numeric(sub("-", "\\.", num))) %>%
  group_by(group) %>%
  arrange(group, num)

C110_Vgene_IGH_hm <- as.data.frame(C110_Vgene_IGH_hm_pre[,-1][,Vgene_ord_IGH$Vgene])
rownames(C110_Vgene_IGH_hm) <- c("E2neg", "d279", "d540", "d1267", "d1842", "mAb")

pheatmap(t(C110_Vgene_IGH_hm), scale="none", cluster_rows = FALSE, cluster_cols=FALSE, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 14, cellwidth = 14, fontsize = 15, fontsize_row = 15, fontsize_col = 16, labels_row = make_italics(Vgene_ord_IGH$Vgene_edit), labels_col = c("E2-", "d279", "d540", "d1267", "d1842", "mAb"), gaps_col=1, border_color=NA, color = colorRampPalette(brewer.pal(9, "Blues"))(50))

## just genes present in mAbs
mAb_Vgenes_IGH <- C110_Vgene_IGH %>%
  filter(time_group == "mAb") %>%
  filter(V_count > 0) 

C110_Vgene_filt_IGH_hm_pre <- C110_Vgene_IGH %>%
  filter(V_gene %in% mAb_Vgenes_IGH$V_gene) %>%
  select(time_group, V_gene, V_prop) %>%
  spread(key=V_gene, value = V_prop) %>%
  arrange(factor(time_group, levels = c("E2neg", "d279", "d540", "d1267", "d1842", "mAb")))

Vgene_ord_IGH_filt <- data.frame(Vgene = colnames(C110_Vgene_filt_IGH_hm_pre[-1]), 
                            Vgene_edit = as.character(gsub("^.{0,4}", "", colnames(C110_Vgene_filt_IGH_hm_pre[-1]))),
                             group = sub("-.*", "", colnames(C110_Vgene_filt_IGH_hm_pre)[-1]),
                             num = sub(".*?-", "", colnames(C110_Vgene_filt_IGH_hm_pre)[-1])) %>%
  mutate(num = as.numeric(sub("-", "\\.", num))) %>%
  group_by(group) %>%
  arrange(group, num)

C110_Vgene_filt_IGH_hm <- as.data.frame(C110_Vgene_filt_IGH_hm_pre[,-1][,Vgene_ord_IGH_filt$Vgene])
rownames(C110_Vgene_filt_IGH_hm) <- c("E2neg", "d279", "d540", "d1267", "d1842", "mAb")

pheatmap(t(C110_Vgene_filt_IGH_hm), scale="none", cluster_rows = FALSE, cluster_cols=FALSE, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 14, cellwidth = 14, fontsize = 15, fontsize_row = 15, fontsize_col = 16, labels_row = make_italics(Vgene_ord_IGH_filt$Vgene_edit), labels_col = c("E2-", "d279", "d540", "d1267", "d1842", "mAb"), gaps_col=1, border_color=NA, color = colorRampPalette(brewer.pal(9, "Blues"))(50))


#### IGKV/IGLV
### make dataframe
C110_IGKLVJ <- rbind.data.frame(
  cbind.data.frame(SeqID = NA, Timepoint = C110_IGKL$Timepoint, Group = C110_IGKL$Group, CDR3 = C110_IGKL$CDR3_aa, V_gene = C110_IGKL$V_gene, J_gene = C110_IGKL$J_gene),
  cbind.data.frame(SeqID = mAbs_IGKL$SeqID, Timepoint = rep("mAb", nrow(mAbs_IGKL)), Group=rep("E2pos", nrow(mAbs_IGKL)), CDR3 = mAbs_IGKL$CDR3_aa, V_gene = mAbs_IGKL$V_gene, J_gene = mAbs_IGKL$J_gene),
  stringsAsFactors = FALSE
)

C110_Vgene_IGKL <- C110_IGKLVJ %>%
  mutate(time_group = ifelse(Group == "E2pos", Timepoint, "E2neg")) %>%   #combine all E2 neg clonotypes into 1 group
  group_by(time_group, V_gene) %>%
  mutate(V_count = n()) %>%
  distinct(V_gene, time_group, V_count) %>%
  ungroup() %>%
  complete(V_gene, time_group, fill = list(V_count =0)) %>%
  group_by(time_group) %>%
  mutate(V_total = sum(V_count)) %>%
  ungroup() %>%
  mutate(V_prop = V_count/V_total)
  
### statistics
C110_IGKLV_list <- unique(C110_Vgene_IGKL$V_gene)

C110_IGKLV_stat_pre <- data.frame()
for(i in 1:length(C110_IGKLV_list)){
  Vgene_df <- filter(C110_Vgene_IGKL, V_gene == C110_IGKLV_list[i]) %>%
    mutate(nonV_count = V_total-V_count)
  df_d279 <- Vgene_df %>%
    filter(time_group == "d279")
  df_d540 <- Vgene_df %>%
    filter(time_group == "d540")
  df_d1267 <- Vgene_df %>%
    filter(time_group == "d1267")
  df_d1842 <- Vgene_df %>%
    filter(time_group == "d1842")
  df_mAb <- Vgene_df %>%
    filter(time_group == "mAb")
  df_E2n <- Vgene_df %>%
    filter(time_group == "E2neg")
  
  d279_E2n <- fisher.test(matrix(c(df_d279$V_count, df_E2n$V_count, df_d279$nonV_count, df_E2n$nonV_count), ncol=2))$p.value
  d540_E2n <- fisher.test(matrix(c(df_d540$V_count, df_E2n$V_count, df_d540$nonV_count, df_E2n$nonV_count), ncol=2))$p.value
  d1267_E2n <- fisher.test(matrix(c(df_d1267$V_count, df_E2n$V_count, df_d1267$nonV_count, df_E2n$nonV_count), ncol=2))$p.value
  d1842_E2n <- fisher.test(matrix(c(df_d1842$V_count, df_E2n$V_count, df_d1842$nonV_count, df_E2n$nonV_count), ncol=2))$p.value
  mAb_E2n <- fisher.test(matrix(c(df_mAb$V_count, df_E2n$V_count, df_mAb$nonV_count, df_E2n$nonV_count), ncol=2))$p.value
  
  d279_d540 <- fisher.test(matrix(c(df_d279$V_count, df_d540$V_count, df_d279$nonV_count, df_d540$nonV_count), ncol=2))$p.value
  d279_d1267 <- fisher.test(matrix(c(df_d279$V_count, df_d1267$V_count, df_d279$nonV_count, df_d1267$nonV_count), ncol=2))$p.value
  d279_d1842 <- fisher.test(matrix(c(df_d279$V_count, df_d1842$V_count, df_d279$nonV_count, df_d1842$nonV_count), ncol=2))$p.value
  d279_mAb <- fisher.test(matrix(c(df_d279$V_count, df_mAb$V_count, df_d279$nonV_count, df_mAb$nonV_count), ncol=2))$p.value
  
  d540_d1267 <- fisher.test(matrix(c(df_d540$V_count, df_d1267$V_count, df_d540$nonV_count, df_d1267$nonV_count), ncol=2))$p.value
  d540_d1842 <- fisher.test(matrix(c(df_d540$V_count, df_d1842$V_count, df_d540$nonV_count, df_d1842$nonV_count), ncol=2))$p.value
  d540_mAb <- fisher.test(matrix(c(df_d540$V_count, df_mAb$V_count, df_d540$nonV_count, df_mAb$nonV_count), ncol=2))$p.value
  
  d1267_d1842 <- fisher.test(matrix(c(df_d1267$V_count, df_d1842$V_count, df_d1267$nonV_count, df_d1842$nonV_count), ncol=2))$p.value
  d1267_mAb <- fisher.test(matrix(c(df_d1267$V_count, df_mAb$V_count, df_d1267$nonV_count, df_mAb$nonV_count), ncol=2))$p.value
  
  d1842_mAb <- fisher.test(matrix(c(df_d1842$V_count, df_mAb$V_count, df_d1842$nonV_count, df_mAb$nonV_count), ncol=2))$p.value
  
  stat_df <- cbind.data.frame(V_gene = C110_IGKLV_list[i], d279_E2n = d279_E2n, d540_E2n = d540_E2n, d1267_E2n = d1267_E2n, d1842_E2n = d1842_E2n, mAb_E2n = mAb_E2n, d279_d540 = d279_d540, d279_d1267 = d279_d1267, d279_d1842 = d279_d1842, d279_mAb = d279_mAb, d540_d1267 = d540_d1267, d540_d1842 = d540_d1842, d540_mAb = d540_mAb, d1267_d1842 = d1267_d1842, d1267_mAb = d1267_mAb, d1842_mAb = d1842_mAb)
  C110_IGKLV_stat_pre <- rbind.data.frame(C110_IGKLV_stat_pre, stat_df, stringsAsFactors = FALSE)
}

C110_IGKLV_stat <- as.data.frame(matrix(p.adjust(as.vector(as.matrix(C110_IGKLV_stat_pre[,-1])), method="fdr"),ncol= ncol(C110_IGKLV_stat_pre)-1))
rownames(C110_IGKLV_stat) <- C110_IGKLV_stat_pre[,1]
colnames(C110_IGKLV_stat) <- colnames(C110_IGKLV_stat_pre[,-1])

C110_IGKLV_stat_sig <- C110_IGKLV_stat %>%
  filter(if_any(d279_E2n:d1842_mAb, ~ .x < 0.05))

### heatmaps
## all genes
C110_Vgene_IGKL_hm_pre <- C110_Vgene_IGKL %>%
  select(time_group, V_gene, V_prop) %>%
  spread(key=V_gene, value = V_prop) %>%
  arrange(factor(time_group, levels = c("E2neg", "d279", "d540", "d1267", "d1842", "mAb")))

Vgene_ord_IGKL <- data.frame(Vgene = colnames(C110_Vgene_IGKL_hm_pre[-1]), 
                            Vgene_edit = as.character(gsub("^.{0,4}", "", colnames(C110_Vgene_IGKL_hm_pre[-1]))),
                             group = sub("-.*", "", colnames(C110_Vgene_IGKL_hm_pre)[-1]),
                             num = sub(".*?-", "", colnames(C110_Vgene_IGKL_hm_pre)[-1])) %>%
  mutate(num = as.numeric(sub("-", "\\.", num))) %>%
  group_by(group) %>%
  arrange(group, num)

C110_Vgene_IGKL_hm <- as.data.frame(C110_Vgene_IGKL_hm_pre[,-1][,Vgene_ord_IGKL$Vgene])
rownames(C110_Vgene_IGKL_hm) <- c("E2neg", "d279", "d540", "d1267", "d1842", "mAb")

pheatmap(t(C110_Vgene_IGKL_hm), scale="none", cluster_rows = FALSE, cluster_cols=FALSE, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 14, cellwidth = 14, fontsize = 15, fontsize_row = 15, fontsize_col = 16, labels_row = make_italics(Vgene_ord_IGKL$Vgene_edit), labels_col = c("E2-", "d279", "d540", "d1267", "d1842", "mAb"), gaps_col=1, gaps_row=29, border_color=NA, color = colorRampPalette(brewer.pal(9, "Blues"))(50))

## just genes present in mAbs
mAb_Vgenes_IGKL <- C110_Vgene_IGKL %>%
  filter(time_group == "mAb") %>%
  filter(V_count > 0) 

C110_Vgene_filt_IGKL_hm_pre <- C110_Vgene_IGKL %>%
  filter(V_gene %in% mAb_Vgenes_IGKL$V_gene) %>%
  select(time_group, V_gene, V_prop) %>%
  spread(key=V_gene, value = V_prop) %>%
  arrange(factor(time_group, levels = c("E2neg", "d279", "d540", "d1267", "d1842", "mAb")))

Vgene_ord_IGKL_filt <- data.frame(Vgene = colnames(C110_Vgene_filt_IGKL_hm_pre[-1]), 
                            Vgene_edit = as.character(sub("^.{0,4}", "", colnames(C110_Vgene_filt_IGKL_hm_pre[-1]))),
                             group = sub("-.*", "", colnames(C110_Vgene_filt_IGKL_hm_pre)[-1]),
                             num = sub(".*?-", "", colnames(C110_Vgene_filt_IGKL_hm_pre)[-1])) %>%
  mutate(num = as.numeric(sub("-", "\\.", num))) %>%
  group_by(group) %>%
  arrange(group, num)

C110_Vgene_filt_IGKL_hm <- as.data.frame(C110_Vgene_filt_IGKL_hm_pre[,-1][,Vgene_ord_IGKL_filt$Vgene])
rownames(C110_Vgene_filt_IGKL_hm) <- c("E2neg", "d279", "d540", "d1267", "d1842", "mAb")

pheatmap(t(C110_Vgene_filt_IGKL_hm), scale="none", cluster_rows = FALSE, cluster_cols=FALSE, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 14, cellwidth = 14, fontsize = 15, fontsize_row = 15, fontsize_col = 16, labels_row = make_italics(Vgene_ord_IGKL_filt$Vgene_edit), labels_col = c("E2-", "d279", "d540", "d1267", "d1842", "mAb"), gaps_col=1, gaps_row=11, border_color=NA, color = colorRampPalette(brewer.pal(9, "Blues"))(50))


###### Somatic hypermutation of the V gene
#### IGHV
### make dataframe
C110_IGH_SHM <- rbind.data.frame(
  cbind.data.frame(SeqID = NA, CDR3_aa = C110_IMGT_IGH$CDR3_aa_IMGT, V_identity = C110_IMGT_IGH$V_identity, Timepoint = C110_IMGT_IGH$Timepoint, Group = C110_IMGT_IGH$Group),
  cbind.data.frame(SeqID = mAbs_IGH$SeqID, CDR3_aa = mAbs_IGH$CDR3_aa, V_identity = mAbs_IGH$V_identity, Timepoint = rep("mAb", nrow(mAbs_IGH)), Group=rep("E2pos")),
  stringsAsFactors = FALSE
) %>%
  filter(V_identity >=85) %>%
  mutate(time_group = ifelse(Group == "E2pos", Timepoint, "E2neg"))

### statistics
shapiro.test(C110_IGH_SHM[C110_IGH_SHM$time_group == "d279",]$V_identity) # normal
shapiro.test(C110_IGH_SHM[C110_IGH_SHM$time_group == "d540",]$V_identity) # normal
shapiro.test(C110_IGH_SHM[C110_IGH_SHM$time_group == "d1267",]$V_identity) # not normal
shapiro.test(C110_IGH_SHM[C110_IGH_SHM$time_group == "d1842",]$V_identity) # normal
shapiro.test(C110_IGH_SHM[C110_IGH_SHM$time_group == "mAb",]$V_identity) # normal

kruskal.test(V_identity ~ time_group, data = C110_IGH_SHM)
C110_IGH_SHM_stat <- dunn_test(C110_IGH_SHM, V_identity ~ time_group, p.adjust.method = "BH") 
  
### plot
C110_IGH_SHM$time_group <- factor(C110_IGH_SHM$time_group, levels = c("E2neg", "d279", "d540", "d1267", "d1842", "mAb"))

C110_IGH_SHM_p <- ggplot(C110_IGH_SHM, aes(x=time_group, y=V_identity, fill=Group, color=Group)) +
  scale_fill_manual(values=c("#cae4f1", "white")) +
  scale_color_manual(values=c("#528aae", "#528aae")) +
  geom_boxplot(width = 0.75, position = position_dodge(width = 0.9), outlier.shape=NA) +
  #stat_pvalue_manual(C110_IGH_SHM_stat, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = , tip.length = 0.005, vjust= 0.75, bracket.size=0.6, size=7) +
  scale_y_continuous(limits = c(85,105))+
  theme(axis.title.x = element_blank(), 
        axis.text=element_text(color="black"),
        axis.text.x=element_text(size=32),
        axis.text.y=element_text(size=32),
        axis.title=element_text(size=36),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        legend.position = "none"
  )+
  labs(x= "",
       y= "% IGHV Gene Identity")
C110_IGH_SHM_p




med_IGH_Vid <- data.frame(x_pos=c(1, 2, 3, 4, 5, 6), 
                       y_pos=c(median(C110_IGH_SHM[C110_IGH_SHM$time_group == "E2neg",]$V_identity),
                               median(C110_IGH_SHM[C110_IGH_SHM$time_group== "d279",]$V_identity),
                               median(C110_IGH_SHM[C110_IGH_SHM$time_group == "d540",]$V_identity),
                               median(C110_IGH_SHM[C110_IGH_SHM$time_group == "d1267",]$V_identity),
                               median(C110_IGH_SHM[C110_IGH_SHM$time_group== "d1842",]$V_identity),
                               median(C110_IGH_SHM[C110_IGH_SHM$time_group == "mAb",]$V_identity)),
                       x_add = c(0.12, 0.12, 0.12, 0.12, 0.12, 0.12))

C110_IGH_SHM_p <- ggplot(C110_IGH_SHM, aes(x=time_group, y=V_identity)) +
  geom_point(aes(fill=time_group, color=time_group), shape =1, position = position_jitter(width=0.175), size=1, stroke=0.95)+
  geom_violin(color= "#808482", show.legend = FALSE, lwd=0.6, alpha=0) +
  geom_segment(data = med_IGH_Vid , aes(x = x_pos-x_add, xend = x_pos+x_add, y = y_pos, yend = y_pos), color = "#676767", linewidth=1.5)+
  scale_shape_manual(values=c(1, 1, 1, 1, 1, 1))+
  scale_color_manual(values=c("#BEBEBE", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#009E73")) +
  scale_x_discrete(breaks=c("E2neg", "d279", "d540", "d1267", "d1842", "mAb"),
                   labels=c("E2-", "d279", "d540", "d1267", "d1842", "mAb"))+
  scale_y_continuous(breaks = seq(85,100, by=5), limits = c(85,104.5))+
  stat_pvalue_manual(C110_IGH_SHM_stat, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(100.5, 101, 101.5, 102.5, 103, 104, 104.5), tip.length = 0.005, vjust= 0.75, bracket.size=0.45, size=5.5) +
  theme(axis.text.y=element_text(color="black", size=16),
        axis.text.x=element_text(color="black", size=16),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.55, color = "black"),
        axis.title=element_text(size=18),
        legend.position = "none",
        aspect.ratio = 1.3)+
  labs(x= "",
       y= "%"~italic(IGHV)~ " Identity")
C110_IGH_SHM_p 


#### IGKV/IGLV
### make dataframe
C110_IGKL_SHM <- rbind.data.frame(
  cbind.data.frame(SeqID = NA, CDR3_aa = C110_IMGT_IGKL$CDR3_aa_IMGT, V_identity = C110_IMGT_IGKL$V_identity, Timepoint = C110_IMGT_IGKL$Timepoint, Group = C110_IMGT_IGKL$Group),
  cbind.data.frame(SeqID = mAbs_IGKL$SeqID, CDR3_aa = mAbs_IGKL$CDR3_aa, V_identity = mAbs_IGKL$V_identity, Timepoint = rep("mAb", nrow(mAbs_IGKL)), Group=rep("E2pos")),
  stringsAsFactors = FALSE
) %>%
  filter(V_identity >= 85) %>%
  mutate(time_group = ifelse(Group == "E2pos", Timepoint, "E2neg"))

### statistics
shapiro.test(C110_IGKL_SHM[C110_IGKL_SHM$time_group == "d279",]$V_identity) # not normal
shapiro.test(C110_IGKL_SHM[C110_IGKL_SHM$time_group == "d540",]$V_identity) # not normal
shapiro.test(C110_IGKL_SHM[C110_IGKL_SHM$time_group == "d1267",]$V_identity) # not normal
shapiro.test(C110_IGKL_SHM[C110_IGKL_SHM$time_group == "d1842",]$V_identity) # normal
shapiro.test(C110_IGKL_SHM[C110_IGKL_SHM$time_group == "mAb",]$V_identity) # not normal

kruskal.test(V_identity ~ time_group, data = C110_IGKL_SHM)
C110_SHM_IGKL_stat <- dunn_test(C110_IGKL_SHM, V_identity ~ time_group, p.adjust.method = "BH") 

### plot
C110_IGKL_SHM$time_group <- factor(C110_IGKL_SHM$time_group, levels = c("E2neg", "d279", "d540", "d1267", "d1842", "mAb"))

med_IGKL_Vid <- data.frame(x_pos=c(1, 2, 3, 4, 5, 6), 
                       y_pos=c(median(C110_IGKL_SHM[C110_IGKL_SHM$time_group == "E2neg",]$V_identity),
                               median(C110_IGKL_SHM[C110_IGKL_SHM$time_group== "d279",]$V_identity),
                               median(C110_IGKL_SHM[C110_IGKL_SHM$time_group == "d540",]$V_identity),
                               median(C110_IGKL_SHM[C110_IGKL_SHM$time_group == "d1267",]$V_identity),
                               median(C110_IGKL_SHM[C110_IGKL_SHM$time_group== "d1842",]$V_identity),
                               median(C110_IGKL_SHM[C110_IGKL_SHM$time_group == "mAb",]$V_identity)),
                       x_add = c(0.12, 0.12, 0.12, 0.12, 0.12, 0.12))

C110_IGKL_SHM_p <- ggplot(C110_IGKL_SHM, aes(x=time_group, y=V_identity)) +
  geom_point(aes(fill=time_group, color=time_group), shape =1, position = position_jitter(width=0.175), size=1, stroke=0.95)+
  geom_violin(color= "#808482", show.legend = FALSE, lwd=0.6, alpha=0) +
  geom_segment(data = med_IGKL_Vid , aes(x = x_pos-x_add, xend = x_pos+x_add, y = y_pos, yend = y_pos), color = "#676767", linewidth=1.5)+
  scale_shape_manual(values=c(1, 1, 1, 1, 1, 1))+
  scale_color_manual(values=c("#BEBEBE", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9", "#009E73")) +
  scale_x_discrete(breaks=c("E2neg", "d279", "d540", "d1267", "d1842", "mAb"),
                   labels=c("E2-", "d279", "d540", "d1267", "d1842", "mAb"))+
  scale_y_continuous(breaks = seq(85,100, by=5), limits = c(85,104.25))+
  stat_pvalue_manual(C110_SHM_IGKL_stat, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(100.5, 101, 102, 102.5, 103.5, 104.25), tip.length = 0.005, vjust= 0.75, bracket.size=0.45, size=5.5) +
  theme(axis.text.y=element_text(color="black", size=16),
        axis.text.x=element_text(color="black", size=16),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.55, color = "black"),
        axis.title=element_text(size=18),
        legend.position = "none",
        aspect.ratio = 1.3)+
  labs(x= "",
       y= "%"~italic(IGKV/IGLV)~ " Identity")
C110_IGKL_SHM_p 


###### Conserved front-layer and non-front-layer bnab mutations compared to E2 negative clonotypes
##### IGHV
#### CDR1/CDR2
### front-layer bnabs
## make df of front-layer bnabs
frly_bnab <- c("hcab3", "hcab15", "hcab21", "hcab30", "hcab43", "hcab55", "hcab64", "hcab71", "HEPC74")
frly_bnab_IGH <- mAbs_IGH %>%
  filter(SeqID %in% frly_bnab) %>%
  mutate(cdr1cdr2_mut = paste(CDR1_mutation, CDR2_mutation, sep = " ")) %>%
  mutate(cdr1_seq = substr(V_aa_seq, start=27, stop=38)) %>%
  mutate(cdr2_seq = substr(V_aa_seq, start=56, stop=65)) %>%
  mutate(cdr1cdr2_seq = paste0(cdr1_seq, cdr2_seq)) %>%
  mutate(cdr3_len = nchar(CDR3_aa))

## filter E2neg clonotypes to be representative of front-layer bnabs (same V genes, same SHM range) 
# df with all E2 neg heavy chain sequences 
E2neg_IMGT_IGH <- C110_IMGT_IGH %>%   
  filter(Group == "E2neg") %>%
  mutate(cdr1_seq = substr(V_aa_seq, start=27, stop=38)) %>%
  mutate(cdr2_seq = substr(V_aa_seq, start=56, stop=65)) %>%
  mutate(cdr1cdr2_seq = paste0(cdr1_seq, cdr2_seq)) %>%
  mutate(cdr1cdr2_seq = gsub("\\.", "x", cdr1cdr2_seq)) %>%
  mutate(gaps_cdr1 = gsub("[A-Z]", "", cdr1_seq)) %>%
  mutate(num_gaps_cdr1 = nchar(gaps_cdr1)) %>%
  mutate(gaps_cdr2 = gsub("[A-Z]", "", cdr2_seq)) %>%
  mutate(num_gaps_cdr2 = nchar(gaps_cdr2)) 

# determine front-layer V gene usage and SHM distribution
frly_feat_dist_IGH <- frly_bnab_IGH %>%
  group_by(V_gene) %>%
  mutate(Vcount = n()) %>%
  mutate(Vprop = Vcount/nrow(.)) %>%
  ungroup() %>%
  mutate(Vid_min = min(V_identity)) %>%
  mutate(Vid_max = max(V_identity)) %>%
  distinct(V_gene, Vcount, Vprop, Vid_min, Vid_max)

# filter E2 neg clonotypes for frly V gene usage and SHM distribution
E2neg_frly_filt_IGH <- E2neg_IMGT_IGH %>%
  mutate(V_gene = V_gene_IMGT) %>%   
  filter(V_gene %in% frly_feat_dist_IGH$V_gene) %>%
  filter(V_identity >= unique(frly_feat_dist_IGH$Vid_min) &
           V_identity <= unique(frly_feat_dist_IGH$Vid_max)) 

## look for mutations enriched in frly bnabs compared to representative E2 neg clonotypes
# define mutations (delta) present in frly bnabs
frly_subs_IGH <- unique(unlist(str_split(frly_bnab_IGH$cdr1cdr2_mut, " ")))
frly_subs_IGH <- frly_subs_IGH[frly_subs_IGH!="none"]
frly_delta_IGH <- unique(sub("[A-Z]", "", frly_subs_IGH))
  
# dataframe of usage
cdr1cdr2 = c(27:38, 56:65)   #IMGT numbering for CDR1 and CDR2 found in any sample
frly_E2n_del_IGH <- stat_mut(frly_bnab_IGH, E2neg_frly_filt_IGH, "all", "delta", cdr1cdr2)[[1]]

# filter E2 neg mutations to only include those found in front-layer bnabs
frly_E2n_filt_IGH <- frly_E2n_del_IGH %>%
  filter(Substitution %in% frly_delta_IGH)

# statistical analysis of mutation usage between front-layer bnabs and representative E2 neg clonotypes
frly_E2n_filt_sig_IGH <- data.frame()
for(i in 1:length(frly_delta_IGH)){
  sub <- frly_delta_IGH[i]
  df <- filter(frly_E2n_filt_IGH, Substitution == sub)
  fisher_pvalue <- fisher.test(matrix(c(df$Count[1], df$Count[2], df$Total[1]-df$Count[1], df$Total[2]-df$Count[2]), ncol=2))$p.value
  df_fish <- data.frame(Substitution = sub, Pvalue = fisher_pvalue, frly_prop = df$Prop[1], E2n_prop = df$Prop[2]) %>%
    mutate(frly_enrich = frly_prop > E2n_prop)
  frly_E2n_filt_sig_IGH <- rbind.data.frame(frly_E2n_filt_sig_IGH, df_fish, stringsAsFactors = FALSE)
}

# fdr correct p-values, filter by adjusted P < 0.05, include Kabat numbering for frly positions
kabat_frly_cdr12_IGH <- data.frame(Position = c(27:30, 35:38, 56:59, 62:65), Kabat = c("26", "27", "28", "29", "30", "31", "32", "33", "51", "52", "52A", "53", "54", "55", "56", "57"))

frly_E2n_filt_sig_corr_IGH <- frly_E2n_filt_sig_IGH %>%
  mutate(Pvalue_adj = p.adjust(Pvalue, "fdr")) %>%
  filter(Pvalue_adj < 0.05) %>%
  mutate(FC = frly_prop/E2n_prop) %>%
  mutate(Position = as.numeric(gsub(">[A-Z]", "", Substitution))) %>%
  mutate(mut = gsub("[0-9][0-9]>", "", Substitution)) %>%
  left_join(kabat_frly_cdr12_IGH, by="Position") %>%
  mutate(Kabat = paste0(Kabat, mut))


### non-front-layer bnabs
## make df of non non-front-layer bnabs (AR1 and AR2-4 epitope binding)
non_frly_bnab <- c("hcab27", "hcab17", "hcab35", "hcab7", "hcab44", "hcab22", "hcab5", "hcab23", "hcab31", "hcab40")

nonfrly_bnab_IGH <- mAbs_IGH %>%
  filter(SeqID %in% non_frly_bnab) %>%
  mutate(epitope = ifelse(SeqID %in% c("hcab27", "hcab17", "hcab35", "hcab7", "hcab44"), "AR1",
                          ifelse(SeqID %in% c("hcab5", "hcab23", "hcab31", "hcab40"), "AR2-4",
                                 ifelse(SeqID == "hcab22", "overlap", "none")))) %>%
  mutate(cdr1cdr2_mut = paste(CDR1_mutation, CDR2_mutation, sep = " ")) %>%
  mutate(cdr1_seq = substr(V_aa_seq, start=27, stop=38)) %>%
  mutate(cdr2_seq = substr(V_aa_seq, start=56, stop=65)) %>%
  mutate(cdr1cdr2_seq = paste0(cdr1_seq, cdr2_seq)) %>%
  mutate(cdr3_len = nchar(CDR3_aa))

# AR1
AR1_bnab_IGH <- nonfrly_bnab_IGH  %>%
  filter(epitope == "AR1" | epitope == "overlap")

# AR2-4
AR24_bnab_IGH <- nonfrly_bnab_IGH  %>%
  filter(epitope == "AR2-4" | epitope == "overlap")

## filter E2neg clonotypes to be representative of non-front-layer bnabs (same V genes, same SHM range) 
# determine nonfront-layer V gene usage and SHM distribution
AR1_feat_dist_IGH <- AR1_bnab_IGH %>%
  group_by(V_gene) %>%
  mutate(Vcount = n()) %>%
  mutate(Vprop = Vcount/nrow(.)) %>%
  ungroup() %>%
  mutate(Vid_min = min(V_identity)) %>%
  mutate(Vid_max = max(V_identity)) %>%
  distinct(V_gene, Vcount, Vprop, Vid_min, Vid_max)

AR24_feat_dist_IGH <- AR24_bnab_IGH %>%
  group_by(V_gene) %>%
  mutate(Vcount = n()) %>%
  mutate(Vprop = Vcount/nrow(.)) %>%
  ungroup() %>%
  mutate(Vid_min = min(V_identity)) %>%
  mutate(Vid_max = max(V_identity)) %>%
  distinct(V_gene, Vcount, Vprop, Vid_min, Vid_max)

# filter E2 neg clonotypes for non-frly V gene usage and SHM distribution
E2neg_AR1_filt_IGH <- E2neg_IMGT_IGH %>%
  mutate(V_gene = V_gene_IMGT) %>%   
  filter(V_gene %in% AR1_feat_dist_IGH$V_gene) %>%
  filter(V_identity >= unique(AR1_feat_dist_IGH$Vid_min) &
           V_identity <= unique(AR1_feat_dist_IGH$Vid_max)) 

E2neg_AR24_filt_IGH <- E2neg_IMGT_IGH %>%
  mutate(V_gene = V_gene_IMGT) %>%   
  filter(V_gene %in% AR24_feat_dist_IGH$V_gene) %>%
  filter(V_identity >= unique(AR24_feat_dist_IGH$Vid_min) &
           V_identity <= unique(AR24_feat_dist_IGH$Vid_max)) 

## define mutations (delta) present in non-frly bnabs
AR1_subs_IGH <- unique(unlist(str_split(AR1_bnab_IGH$cdr1cdr2_mut, " ")))
AR1_subs_IGH <- AR1_subs_IGH[AR1_subs_IGH!="none"]
AR1_delta_IGH <- unique(sub("[A-Z]", "", AR1_subs_IGH))

AR24_subs_IGH <- unique(unlist(str_split(AR24_bnab_IGH$cdr1cdr2_mut, " ")))
AR24_subs_IGH <- AR24_subs_IGH[AR24_subs_IGH!="none"]
AR24_delta_IGH <- unique(sub("[A-Z]", "", AR24_subs_IGH))

  
# dataframe of usage
AR1_E2n_del_IGH <- stat_mut(AR1_bnab_IGH, E2neg_AR1_filt_IGH, "all", "delta", cdr1cdr2)[[1]]

AR24_E2n_del_IGH <- stat_mut(AR24_bnab_IGH, E2neg_AR24_filt_IGH, "all", "delta", cdr1cdr2)[[1]]

# filter E2 neg mutations to only include those found in non-front-layer bnabs
AR1_E2n_filt_IGH <- AR1_E2n_del_IGH %>%
  filter(Substitution %in% AR1_delta_IGH)

AR24_E2n_filt_IGH <- AR24_E2n_del_IGH %>%
  filter(Substitution %in% AR24_delta_IGH)

# statistical analysis of mutation usage between non-front-layer bnabs and representative E2 neg clonotypes
AR1_E2n_filt_sig_IGH <- data.frame()
for(i in 1:length(AR1_delta_IGH)){
  sub <- AR1_delta_IGH[i]
  df <- filter(AR1_E2n_filt_IGH, Substitution == sub)
  fisher_pvalue <- fisher.test(matrix(c(df$Count[1], df$Count[2], df$Total[1]-df$Count[1], df$Total[2]-df$Count[2]), ncol=2))$p.value
  df_fish <- data.frame(Substitution = sub, Pvalue = fisher_pvalue, AR1_prop = df$Prop[1], E2n_prop = df$Prop[2]) %>%
    mutate(AR1_enrich = AR1_prop > E2n_prop)
  AR1_E2n_filt_sig_IGH <- rbind.data.frame(AR1_E2n_filt_sig_IGH, df_fish, stringsAsFactors = FALSE)
}

AR24_E2n_filt_sig_IGH <- data.frame()
for(i in 1:length(AR24_delta_IGH)){
  sub <- AR24_delta_IGH[i]
  df <- filter(AR24_E2n_filt_IGH, Substitution == sub)
  fisher_pvalue <- fisher.test(matrix(c(df$Count[1], df$Count[2], df$Total[1]-df$Count[1], df$Total[2]-df$Count[2]), ncol=2))$p.value
  df_fish <- data.frame(Substitution = sub, Pvalue = fisher_pvalue, AR24_prop = df$Prop[1], E2n_prop = df$Prop[2]) %>%
    mutate(AR24_enrich = AR24_prop > E2n_prop)
  AR24_E2n_filt_sig_IGH <- rbind.data.frame(AR24_E2n_filt_sig_IGH, df_fish, stringsAsFactors = FALSE)
}

# fdr correct p-values, filter by adjusted P < 0.05, include Kabat numbering
kabat_nonfrly_cdr12_IGH <- data.frame(Position = c(27:31, 34:38, 56:59, 62:65), Kabat = c("26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "51", "52", "52A", "53", "54", "55", "56", "57"))

AR1_E2n_filt_sig_corr_IGH <- AR1_E2n_filt_sig_IGH %>%
  mutate(Pvalue_adj = p.adjust(Pvalue, "fdr")) %>%
  filter(Pvalue_adj < 0.05) %>%
  mutate(FC = AR1_prop/E2n_prop) %>%
  mutate(Position = as.numeric(gsub(">[A-Z]", "", Substitution))) %>%
  mutate(mut = gsub("[0-9][0-9]>", "", Substitution)) %>%
  left_join(kabat_nonfrly_cdr12_IGH, by="Position") %>%
  mutate(Kabat = paste0(Kabat, mut))

AR24_E2n_filt_sig_corr_IGH <- AR24_E2n_filt_sig_IGH %>%
  mutate(Pvalue_adj = p.adjust(Pvalue, "fdr")) %>%
  filter(Pvalue_adj < 0.05) %>%
  mutate(FC = AR24_prop/E2n_prop) %>%
  mutate(Position = as.numeric(gsub(">[A-Z]", "", Substitution))) %>%
  mutate(mut = gsub("[0-9][0-9]>", "", Substitution)) %>%
  left_join(kabat_nonfrly_cdr12_IGH, by="Position") %>%
  mutate(Kabat = paste0(Kabat, mut))


#### CDR3
# align front-layer and non-front-layer bnab CDR3s with a representative subset of E2 negative clonotypes

### front-layer
## filter E2neg clonotypes to be representative of front-layer bnabs (same V genes, same SHM range, same CDR3 lengths) 
E2n_frly_filt_IGH <- C110_IMGT_IGH %>%
  filter(Group == "E2neg") %>%
  filter(V_gene_IMGT %in% unique(frly_bnab_IGH$V_gene)) %>%
  filter(J_gene_IMGT %in% unique(frly_bnab_IGH$J_gene)) %>%
  mutate(CDR3_len = nchar(CDR3_aa_IMGT)) %>%
  filter(CDR3_len %in% unique(frly_bnab_IGH$cdr3_len)) %>%
  filter(V_identity >= min(frly_bnab_IGH$V_identity) &
           V_identity <= max(frly_bnab_IGH$V_identity)) %>%
  select(Group, CDR3_aa = CDR3_aa_IMGT, V_gene_IMGT, J_gene_IMGT, V_identity)


## export table of E2neg clonotypes and front-layer bnabs and run alignment for each on the cluster
  # write.table(E2n_frly_filt_IGH , "~/Desktop/E2n_frly_filt_IGH.txt", quote=FALSE , sep = "\t", row.names = FALSE)
  # write.table(frly_bnab_IGH, "~/Desktop/frly_bnab_IGH.txt.txt", quote=FALSE , sep = "\t", row.names = FALSE)

## import alignments and combine to one dataframe
frly_cdr3_align_IGH <- read.csv("~/Desktop/Serial_C110/Final_CDR3_alignments/frly_bnab_alignment.csv") %>%
  mutate(Group = "frly_bnab") %>%
  select(Seq, Group)

E2n_frlyfilt_cdr3_align_IGH <- read.csv("~/Desktop/Serial_C110/Final_CDR3_alignments/E2n_frly_align_IGH.csv") %>%
  mutate(Group = "E2neg") %>%
  select(Seq, Group)

E2neg_frly_cdr3_align <- rbind.data.frame(frly_cdr3_align_IGH, E2n_frlyfilt_cdr3_align_IGH, stringsAsFactors = FALSE)


## Determine amino acid usage
# unable to look for substitutions since no consensus germline CDR3 exists; consequently will look at amino acids used by front-layer bnabs at each position and determine if usage is statistically significantly higher than for E2 neg clonotypes 
frly_cdr3_usage <- AA_usage(E2neg_frly_cdr3_align[E2neg_frly_cdr3_align$Group == "frly_bnab",]$Seq, nchar(E2neg_frly_cdr3_align$Seq[1]))
E2n_frlyfilt_cdr3_usage <- AA_usage(E2neg_frly_cdr3_align[E2neg_frly_cdr3_align$Group == "E2neg",]$Seq, nchar(E2neg_frly_cdr3_align$Seq[1]))

## filter to only include analysis of amino acids present in front-layer bnabs; then compare usage and obtain P values
E2neg_frly_cdr3_align_stat_pre <- data.frame()
for(i in 1:nchar(E2neg_frly_cdr3_align$Seq[1])){
  frly_filt <- frly_cdr3_usage[,i]
  frly_filt2 <- frly_filt[frly_filt!=0]
  E2n_filt <- E2n_frlyfilt_cdr3_usage[,i] 
  E2n_filt2 <- E2n_filt[names(E2n_filt) %in% names(frly_filt2)]

  amino_acids <- names(frly_filt2)
  
  frly_E2n_fish <- data.frame()
  for(j in 1:length(amino_acids)){
    AA <- amino_acids[j]
    AA_count_frly <- frly_filt[names(frly_filt) == AA]
    AA_count_E2n <- E2n_filt[names(E2n_filt) == AA]
    AA_total_frly <- nrow(frly_cdr3_align_IGH)
    AA_total_E2n <- nrow(E2n_frlyfilt_cdr3_align_IGH )
    fish_mat <-matrix(c(AA_count_frly, AA_count_E2n, AA_total_frly-AA_count_frly, AA_total_E2n-AA_count_E2n), ncol = 2)
    fish_test <- cbind.data.frame(AA = AA, 
                                  frly_prop = AA_count_frly/AA_total_frly, 
                                  E2n_prop = AA_count_E2n/AA_total_E2n,
                                  P_value = fisher.test(fish_mat)$p.value)  
    frly_E2n_fish <- rbind.data.frame(frly_E2n_fish, fish_test, stringsAsFactors = FALSE)
  }
  
  frly_E2n_fish2 <- frly_E2n_fish %>%
    mutate(Position = i)
  
  E2neg_frly_cdr3_align_stat_pre <- rbind.data.frame(E2neg_frly_cdr3_align_stat_pre, frly_E2n_fish2, stringsAsFactors = FALSE)
}

## fdr correct p-values, convert to Kabat numbering, find fold-change
# kabat numbering
kabat_frly_cdr3_IGH <- data.frame(Position = c(1:22), Kabat = c("92", "93", "94", "95", "96", "97", "98", "99", "100", "100a", "100b", "100c", "100d", "100e", "100f", "100g", "100h", "100i", "100j", "101", "102", "103"))

# for fold-change calculation, have to replace 0s (replace with  1/2 of the smallest non-zero number)
cdr3_frly_zero_corr_IGH <- min(  # change zeros to a number close to zero (smallest non-zero number divded by 2) to avoid dividing by infinity
  min(E2neg_frly_cdr3_align_stat_pre$frly_prop[E2neg_frly_cdr3_align_stat_pre$frly_prop != 0])/2,
  min(E2neg_frly_cdr3_align_stat_pre$E2n_prop[E2neg_frly_cdr3_align_stat_pre$E2n_prop != 0])/2)

# fdr correction
E2neg_frly_cdr3_align_stat <- E2neg_frly_cdr3_align_stat_pre %>%
  left_join(kabat_frly_cdr3_IGH, by= "Position") %>%
  mutate(frly_prop_adj = ifelse(frly_prop == 0, cdr3_frly_zero_corr_IGH, frly_prop)) %>%
  mutate(E2n_prop_adj = ifelse(E2n_prop == 0, cdr3_frly_zero_corr_IGH, E2n_prop)) %>%
  mutate(FC = frly_prop_adj/E2n_prop_adj) %>%
  mutate(P_value_adj = p.adjust(P_value, "BH"))

# filter for significantly differences
E2neg_frly_cdr3_align_sig <- E2neg_frly_cdr3_align_stat  %>%
  filter(P_value_adj < 0.05) 

### non-front-layer
# For AR1, need toemove HCAB 35 with CDR3 length 31 to help with alignment run time
AR1_bnab_no35_IGH <- AR1_bnab_IGH %>%
  filter(SeqID != "hcab35")

## filter E2neg clonotypes to be representative of non-front-layer bnabs (same V genes, same SHM range, same CDR3 lengths) 
E2n_AR1_filt_IGH <- C110_IMGT_IGH %>%
  filter(Group == "E2neg") %>%
  filter(V_gene_IMGT %in% unique(AR1_bnab_no35_IGH$V_gene)) %>%
  filter(J_gene_IMGT %in% unique(AR1_bnab_no35_IGH$J_gene)) %>%
  mutate(CDR3_len = nchar(CDR3_aa_IMGT)) %>%
  filter(CDR3_len %in% unique(AR1_bnab_no35_IGH$cdr3_len)) %>%
  filter(V_identity >= min(AR1_bnab_no35_IGH$V_identity) &
           V_identity <= max(AR1_bnab_no35_IGH$V_identity)) %>%
  select(Group, CDR3_aa = CDR3_aa_IMGT, V_gene_IMGT, J_gene_IMGT, V_identity)

E2n_AR24_filt_IGH <- C110_IMGT_IGH %>%
  filter(Group == "E2neg") %>%
  filter(V_gene_IMGT %in% unique(AR24_bnab_IGH$V_gene)) %>%
  filter(J_gene_IMGT %in% unique(AR24_bnab_IGH$J_gene)) %>%
  mutate(CDR3_len = nchar(CDR3_aa_IMGT)) %>%
  filter(CDR3_len %in% unique(AR24_bnab_IGH$cdr3_len)) %>%
  filter(V_identity >= min(AR24_bnab_IGH$V_identity) &
           V_identity <= max(AR24_bnab_IGH$V_identity)) %>%
  select(Group, CDR3_aa = CDR3_aa_IMGT, V_gene_IMGT, J_gene_IMGT, V_identity)


## export table of E2neg clonotypes and non-front-layer bnabs and run alignment for each on the cluster
# write.table(E2n_AR1_filt_IGH, "~/Desktop/E2n_AR1_filt_IGH.txt", quote=FALSE , sep = "\t", row.names = FALSE)
# write.table(E2n_AR24_filt_IGH, "~/Desktop/E2n_AR24_filt_IGH.txt", quote=FALSE , sep = "\t", row.names = FALSE)

## import alignments and combine to one dataframe
AR1_cdr3_align_IGH <- read.csv("~/Desktop/Serial_C110/Final_CDR3_alignments/AR1_bnab_align_IGH.csv") %>%
  mutate(Group = "bnab") %>%
  select(Seq, Group)

E2n_AR1_cdr3_align_IGH <- read.csv("~/Desktop/Serial_C110/Final_CDR3_alignments/E2n_AR1_align_IGH.csv") %>%
  mutate(Group = "E2neg") %>%
  select(Seq, Group)


AR24_cdr3_align_IGH <- read.csv("~/Desktop/Serial_C110/Final_CDR3_alignments/AR24_IGH_alignment.csv") %>%
  mutate(Group = "bnab") %>%
  select(Seq, Group)

E2n_AR24_cdr3_align_IGH <- read.csv("~/Desktop/Serial_C110/Final_CDR3_alignments/E2n_AR24_align_IGH.csv") %>%
  mutate(Group = "E2neg") %>%
  select(Seq, Group)


E2neg_AR1_cdr3_align_IGH <- rbind.data.frame(AR1_cdr3_align_IGH, E2n_AR1_cdr3_align_IGH, stringsAsFactors = FALSE)

E2neg_AR24_cdr3_align_IGH <- rbind.data.frame(AR24_cdr3_align_IGH, E2n_AR24_cdr3_align_IGH, stringsAsFactors = FALSE)

## Determine amino acid usage
# unable to look for substitutions since no consensus germline CDR3 exists; consequently will look at amino acids used by non-front-layer bnabs at each position and determine if usage is statistically significantly higher than for E2 neg clonotypes 
AR1_cdr3_usage <- AA_usage(E2neg_AR1_cdr3_align_IGH[E2neg_AR1_cdr3_align_IGH$Group == "bnab",]$Seq, nchar(E2neg_AR1_cdr3_align_IGH$Seq[1]))
E2n_AR1_cdr3_usage <- AA_usage(E2neg_AR1_cdr3_align_IGH[E2neg_AR1_cdr3_align_IGH$Group == "E2neg",]$Seq, nchar(E2neg_AR1_cdr3_align_IGH$Seq[1]))

AR24_cdr3_usage <- AA_usage(E2neg_AR24_cdr3_align_IGH[E2neg_AR24_cdr3_align_IGH$Group == "bnab",]$Seq, nchar(E2neg_AR24_cdr3_align_IGH$Seq[1]))
E2n_AR24_cdr3_usage <- AA_usage(E2neg_AR24_cdr3_align_IGH[E2neg_AR24_cdr3_align_IGH$Group == "E2neg",]$Seq, nchar(E2neg_AR24_cdr3_align_IGH$Seq[1]))

## filter to only include analysis of amino acids present in bnabs; then compare usage and obtain P values
## AR1
E2neg_AR1_cdr3_align_IGH_stat_pre <- data.frame()
for(i in 1:nchar(E2neg_AR1_cdr3_align_IGH$Seq[1])){
  AR1_filt <- AR1_cdr3_usage[,i]
  AR1_filt2 <- AR1_filt[AR1_filt!=0]
  E2n_filt <- E2n_AR1_cdr3_usage[,i] 
  E2n_filt2 <- E2n_filt[names(E2n_filt) %in% names(AR1_filt2)]

  amino_acids <- names(AR1_filt2)
  
  AR1_E2n_fish <- data.frame()
  for(j in 1:length(amino_acids)){
    AA <- amino_acids[j]
    AA_count_AR1 <- AR1_filt[names(AR1_filt) == AA]
    AA_count_E2n <- E2n_filt[names(E2n_filt) == AA]
    AA_total_AR1 <- nrow(AR1_cdr3_align_IGH)
    AA_total_E2n <- nrow(E2n_AR1_cdr3_align_IGH)
    fish_mat <-matrix(c(AA_count_AR1, AA_count_E2n, AA_total_AR1-AA_count_AR1, AA_total_E2n-AA_count_E2n), ncol = 2)
    fish_test <- cbind.data.frame(AA = AA, 
                                  AR1_prop = AA_count_AR1/AA_total_AR1, 
                                  E2n_prop = AA_count_E2n/AA_total_E2n,
                                  P_value = fisher.test(fish_mat)$p.value)  
    AR1_E2n_fish <- rbind.data.frame(AR1_E2n_fish, fish_test, stringsAsFactors = FALSE)
  }
  
  AR1_E2n_fish2 <- AR1_E2n_fish %>%
    mutate(Position = i)
  
  E2neg_AR1_cdr3_align_IGH_stat_pre <- rbind.data.frame(E2neg_AR1_cdr3_align_IGH_stat_pre, AR1_E2n_fish2, stringsAsFactors = FALSE)
}

## AR2-4
E2neg_AR24_cdr3_align_IGH_stat_pre <- data.frame()
for(i in 1:nchar(E2neg_AR24_cdr3_align_IGH$Seq[1])){
  AR24_filt <- AR24_cdr3_usage[,i]
  AR24_filt2 <- AR24_filt[AR24_filt!=0]
  E2n_filt <- E2n_AR24_cdr3_usage[,i] 
  E2n_filt2 <- E2n_filt[names(E2n_filt) %in% names(AR24_filt2)]

  amino_acids <- names(AR24_filt2)
  
  AR24_E2n_fish <- data.frame()
  for(j in 1:length(amino_acids)){
    AA <- amino_acids[j]
    AA_count_AR24 <- AR24_filt[names(AR24_filt) == AA]
    AA_count_E2n <- E2n_filt[names(E2n_filt) == AA]
    AA_total_AR24 <- nrow(AR24_cdr3_align_IGH)
    AA_total_E2n <- nrow(E2n_AR24_cdr3_align_IGH)
    fish_mat <-matrix(c(AA_count_AR24, AA_count_E2n, AA_total_AR24-AA_count_AR24, AA_total_E2n-AA_count_E2n), ncol = 2)
    fish_test <- cbind.data.frame(AA = AA, 
                                  AR24_prop = AA_count_AR24/AA_total_AR24, 
                                  E2n_prop = AA_count_E2n/AA_total_E2n,
                                  P_value = fisher.test(fish_mat)$p.value)  
    AR24_E2n_fish <- rbind.data.frame(AR24_E2n_fish, fish_test, stringsAsFactors = FALSE)
  }
  
  AR24_E2n_fish2 <- AR24_E2n_fish %>%
    mutate(Position = i)
  
  E2neg_AR24_cdr3_align_IGH_stat_pre <- rbind.data.frame(E2neg_AR24_cdr3_align_IGH_stat_pre, AR24_E2n_fish2, stringsAsFactors = FALSE)
}

## FDR correct p-values, convert to Kabat numbering, find fold-change

# kabat numbering
kabat_AR1_cdr3_IGH <- data.frame(Position = c(1:26), Kabat = c("92", "93", "94", "95", "96", "97", "98", "99", "100", "100a", "100b", "100c", "100d", "100e", "100f", "100g", "100h", "100i", "100j", "100k", "100l", "100m", "100n", "101", "102", "103"))

# for fold-change calculation, have to replace 0s (replace with  1/2 of the smallest non-zero number)
cdr3_AR1_zero_corr_IGH <- min(
  min(E2neg_AR1_cdr3_align_IGH_stat_pre$AR1_prop[E2neg_AR1_cdr3_align_IGH_stat_pre$AR1_prop != 0])/2,
  min(E2neg_AR1_cdr3_align_IGH_stat_pre$E2n_prop[E2neg_AR1_cdr3_align_IGH_stat_pre$E2n_prop != 0])/2)

# fdr correction
E2neg_AR1_cdr3_align_IGH_stat <- E2neg_AR1_cdr3_align_IGH_stat_pre %>%
  left_join(kabat_AR1_cdr3_IGH, by= "Position") %>%
  mutate(AR1_prop_adj = ifelse(AR1_prop == 0, cdr3_AR1_zero_corr_IGH, AR1_prop)) %>%
  mutate(E2n_prop_adj = ifelse(E2n_prop == 0, cdr3_AR1_zero_corr_IGH, E2n_prop)) %>%
  mutate(FC = AR1_prop_adj/E2n_prop_adj) %>%
  mutate(P_value_adj = p.adjust(P_value, "fdr"))

# filter for significantly differences
E2neg_AR1_cdr3_align_IGH_sig <- E2neg_AR1_cdr3_align_IGH_stat  %>%
  filter(P_value_adj < 0.05)


## AR2-4
# kabat numbering
kabat_AR24_cdr3_IGH <- data.frame(Position = c(1:25), Kabat = c("92", "93", "94", "95", "96", "97", "98", "99", "100", "100a", "100b", "100c", "100d", "100e", "100f", "100g", "100h", "100i", "100j", "100k", "100l", "100m", "101", "102", "103"))

# for fold-change calculation, have to replace 0s (replace with  1/2 of the smallest non-zero number)
cdr3_AR24_zero_corr_IGH <- min(
  min(E2neg_AR24_cdr3_align_IGH_stat_pre$AR24_prop[E2neg_AR24_cdr3_align_IGH_stat_pre$AR24_prop != 0])/2,
  min(E2neg_AR24_cdr3_align_IGH_stat_pre$E2n_prop[E2neg_AR24_cdr3_align_IGH_stat_pre$E2n_prop != 0])/2)

# fdr correction
E2neg_AR24_cdr3_align_IGH_stat <- E2neg_AR24_cdr3_align_IGH_stat_pre %>%
  left_join(kabat_AR24_cdr3_IGH, by= "Position") %>%
  mutate(AR24_prop_adj = ifelse(AR24_prop == 0, cdr3_AR24_zero_corr_IGH, AR24_prop)) %>%
  mutate(E2n_prop_adj = ifelse(E2n_prop == 0, cdr3_AR24_zero_corr_IGH, E2n_prop)) %>%
  mutate(FC = AR24_prop_adj/E2n_prop_adj) %>%
  mutate(P_value_adj = p.adjust(P_value, "fdr"))

# filter for significantly differences
E2neg_AR24_cdr3_align_IGH_sig <- E2neg_AR24_cdr3_align_IGH_stat  %>%
  filter(P_value_adj < 0.05)


##### IGK/IGL
#### CDR1/CDR2
### front-layer
frly_bnab_IGKL <- mAbs_IGKL %>%
  filter(SeqID %in% frly_bnab) %>%
  mutate(cdr1cdr2_mut = paste(CDR1_mutation, CDR2_mutation, sep = " ")) %>%
  mutate(cdr1_seq = substr(V_aa_seq, start=27, stop=38)) %>%
  mutate(cdr2_seq = substr(V_aa_seq, start=56, stop=65)) %>%
  mutate(cdr1cdr2_seq = paste0(cdr1_seq, cdr2_seq)) %>%
  mutate(cdr3_len = nchar(CDR3_aa))

## filter E2neg clonotypes to be representative of front-layer bnabs (same V genes, same SHM range) 
# df with all E2 neg heavy chain sequences 
E2neg_IMGT_IGKL <- C110_IMGT_IGKL %>%   
  filter(Group == "E2neg") %>%
  mutate(cdr1_seq = substr(V_aa_seq, start=27, stop=38)) %>%
  mutate(cdr2_seq = substr(V_aa_seq, start=56, stop=65)) %>%
  mutate(cdr1cdr2_seq = paste0(cdr1_seq, cdr2_seq)) %>%
  mutate(cdr1cdr2_seq = gsub("\\.", "x", cdr1cdr2_seq)) %>%
  mutate(gaps_cdr1 = gsub("[A-Z]", "", cdr1_seq)) %>%
  mutate(num_gaps_cdr1 = nchar(gaps_cdr1)) %>%
  mutate(gaps_cdr2 = gsub("[A-Z]", "", cdr2_seq)) %>%
  mutate(num_gaps_cdr2 = nchar(gaps_cdr2)) 

# determine front-layer V gene usage and SHM distribution
frly_feat_dist_IGKL <- frly_bnab_IGKL %>%
  group_by(V_gene) %>%
  mutate(Vcount = n()) %>%
  mutate(Vprop = Vcount/nrow(.)) %>%
  ungroup() %>%
  mutate(Vid_min = min(V_identity)) %>%
  mutate(Vid_max = max(V_identity)) %>%
  distinct(V_gene, Vcount, Vprop, Vid_min, Vid_max)

# filter E2 neg clonotypes for frly V gene usage and SHM distribution
E2neg_frly_filt_IGKL <- E2neg_IMGT_IGKL %>%
  mutate(V_gene = V_gene_IMGT) %>%   
  filter(V_gene %in% frly_feat_dist_IGKL$V_gene) %>%
  filter(V_identity >= unique(frly_feat_dist_IGKL$Vid_min) &
           V_identity <= unique(frly_feat_dist_IGKL$Vid_max)) 

## look for mutations enriched in frly bnabs compared to representative E2 neg clonotypes
# define mutations (delta) present in frly bnabs
frly_subs_IGKL <- unique(unlist(str_split(frly_bnab_IGKL$cdr1cdr2_mut, " ")))
frly_subs_IGKL <- frly_subs_IGKL[frly_subs_IGKL!="none"]
frly_delta_IGKL <- unique(sub("[A-Z]", "", frly_subs_IGKL))
  
# dataframe of usage
frly_E2n_del_IGKL <- stat_mut(frly_bnab_IGKL, E2neg_frly_filt_IGKL, "all", "delta", cdr1cdr2)[[1]]

# filter E2 neg mutations to only include those found in front-layer bnabs
frly_E2n_filt_IGKL <- frly_E2n_del_IGKL %>%
  filter(Substitution %in% frly_delta_IGKL)

# statistical analysis of mutation usage between front-layer bnabs and representative E2 neg clonotypes
frly_E2n_filt_sig_IGKL <- data.frame()
for(i in 1:length(frly_delta_IGKL)){
  sub <- frly_delta_IGKL[i]
  df <- filter(frly_E2n_filt_IGKL, Substitution == sub)
  fisher_pvalue <- fisher.test(matrix(c(df$Count[1], df$Count[2], df$Total[1]-df$Count[1], df$Total[2]-df$Count[2]), ncol=2))$p.value
  df_fish <- data.frame(Substitution = sub, Pvalue = fisher_pvalue, frly_prop = df$Prop[1], E2n_prop = df$Prop[2]) %>%
    mutate(frly_enrich = frly_prop > E2n_prop)
  frly_E2n_filt_sig_IGKL <- rbind.data.frame(frly_E2n_filt_sig_IGKL, df_fish, stringsAsFactors = FALSE)
}

# fdr correct p-values, filter by adjusted P < 0.05
frly_E2n_filt_sig_corr_IGKL <- frly_E2n_filt_sig_IGKL %>%
  mutate(Pvalue_adj = p.adjust(Pvalue, "fdr")) %>%
  filter(Pvalue_adj < 0.05) %>%
  mutate(FC = frly_prop/E2n_prop) %>%
  mutate(Position = as.numeric(gsub(">[A-Z]", "", Substitution))) %>%
  mutate(mut = gsub("[0-9][0-9]>", "", Substitution))
  
### non-front-layer
## AR1
AR1_bnab_IGKL <- mAbs_IGKL %>%
  filter(SeqID %in% AR1_bnab_IGH$SeqID) %>%
  mutate(cdr1cdr2_mut = paste(CDR1_mutation, CDR2_mutation, sep = " ")) %>%
  mutate(cdr1_seq = substr(V_aa_seq, start=27, stop=38)) %>%
  mutate(cdr2_seq = substr(V_aa_seq, start=56, stop=65)) %>%
  mutate(cdr1cdr2_seq = paste0(cdr1_seq, cdr2_seq)) %>%
  mutate(cdr3_len = nchar(CDR3_aa))

## determine AR1 V gene usage and SHM distribution
AR1_feat_dist_IGKL <- AR1_bnab_IGKL %>%
  group_by(V_gene) %>%
  mutate(Vcount = n()) %>%
  mutate(Vprop = Vcount/nrow(.)) %>%
  ungroup() %>%
  mutate(Vid_min = min(V_identity)) %>%
  mutate(Vid_max = max(V_identity)) %>%
  distinct(V_gene, Vcount, Vprop, Vid_min, Vid_max)

# filter E2 neg clonotypes for AR1 V gene usage and SHM distribution
E2neg_AR1_filt_IGKL <- E2neg_IMGT_IGKL %>%
  mutate(V_gene = V_gene_IMGT) %>%   
  filter(V_gene %in% AR1_feat_dist_IGKL$V_gene) %>%
  filter(V_identity >= unique(AR1_feat_dist_IGKL$Vid_min) &
           V_identity <= unique(AR1_feat_dist_IGKL$Vid_max)) 

## define mutations (delta) present in AR1 bnabs
AR1_subs_IGKL <- unique(unlist(str_split(AR1_bnab_IGKL$cdr1cdr2_mut, " ")))
AR1_subs_IGKL <- AR1_subs_IGKL[AR1_subs_IGKL!="none"]
AR1_delta_IGKL <- unique(sub("[A-Z]", "", AR1_subs_IGKL))
  
# dataframe of usage
AR1_E2n_del_IGKL <- stat_mut(AR1_bnab_IGKL, E2neg_AR1_filt_IGKL, "all", "delta", cdr1cdr2)[[1]]

# filter E2 neg mutations to only include those found in AR1 bnabs
AR1_E2n_filt_IGKL <- AR1_E2n_del_IGKL %>%
  filter(Substitution %in% AR1_delta_IGKL)

# statistical analysis of mutation usage between AR1 bnabs and representative E2 neg clonotypes
AR1_E2n_filt_sig_IGKL <- data.frame()
for(i in 1:length(AR1_delta_IGKL)){
  sub <- AR1_delta_IGKL[i]
  df <- filter(AR1_E2n_filt_IGKL, Substitution == sub)
  fisher_pvalue <- fisher.test(matrix(c(df$Count[1], df$Count[2], df$Total[1]-df$Count[1], df$Total[2]-df$Count[2]), ncol=2))$p.value
  df_fish <- data.frame(Substitution = sub, Pvalue = fisher_pvalue, AR1_prop = df$Prop[1], E2n_prop = df$Prop[2]) %>%
    mutate(AR1_enrich = AR1_prop > E2n_prop)
  AR1_E2n_filt_sig_IGKL <- rbind.data.frame(AR1_E2n_filt_sig_IGKL, df_fish, stringsAsFactors = FALSE)
}

# fdr correct p-values, filter by adjusted P < 0.05, include Kabat numbering
AR1_E2n_filt_sig_corr_IGKL <- AR1_E2n_filt_sig_IGKL %>%
  mutate(Pvalue_adj = p.adjust(Pvalue, "fdr")) %>%
  filter(Pvalue_adj < 0.05) %>%
  mutate(FC = AR1_prop/E2n_prop) %>%
  mutate(Position = as.numeric(gsub(">[A-Z]", "", Substitution))) %>%
  mutate(mut = gsub("[0-9][0-9]>", "", Substitution))
  

## AR2-4
AR24_bnab_IGKL <- mAbs_IGKL %>%
  filter(SeqID %in% AR24_bnab_IGH$SeqID) %>%
  mutate(cdr1cdr2_mut = paste(CDR1_mutation, CDR2_mutation, sep = " ")) %>%
  mutate(cdr1_seq = substr(V_aa_seq, start=27, stop=38)) %>%
  mutate(cdr2_seq = substr(V_aa_seq, start=56, stop=65)) %>%
  mutate(cdr1cdr2_seq = paste0(cdr1_seq, cdr2_seq)) %>%
  mutate(cdr3_len = nchar(CDR3_aa))

## determine AR2-4 V gene usage and SHM distribution
AR24_feat_dist_IGKL <- AR24_bnab_IGKL %>%
  group_by(V_gene) %>%
  mutate(Vcount = n()) %>%
  mutate(Vprop = Vcount/nrow(.)) %>%
  ungroup() %>%
  mutate(Vid_min = min(V_identity)) %>%
  mutate(Vid_max = max(V_identity)) %>%
  distinct(V_gene, Vcount, Vprop, Vid_min, Vid_max)

# filter E2 neg clonotypes for AR2-4 V gene usage and SHM distribution
E2neg_AR24_filt_IGKL <- E2neg_IMGT_IGKL %>%
  mutate(V_gene = V_gene_IMGT) %>%   
  filter(V_gene %in% AR24_feat_dist_IGKL$V_gene) %>%
  filter(V_identity >= unique(AR24_feat_dist_IGKL$Vid_min) &
           V_identity <= unique(AR24_feat_dist_IGKL$Vid_max)) 

## define mutations (delta) present in AR24 bnabs
AR24_subs_IGKL <- unique(unlist(str_split(AR24_bnab_IGKL$cdr1cdr2_mut, " ")))
AR24_subs_IGKL <- AR24_subs_IGKL[AR24_subs_IGKL!="none"]
AR24_delta_IGKL <- unique(sub("[A-Z]", "", AR24_subs_IGKL))
  
# dataframe of usage
AR24_E2n_del_IGKL <- stat_mut(AR24_bnab_IGKL, E2neg_AR24_filt_IGKL, "all", "delta", cdr1cdr2)[[1]]

# filter E2 neg mutations to only include those found in AR2-4 bnabs
AR24_E2n_filt_IGKL <- AR24_E2n_del_IGKL %>%
  filter(Substitution %in% AR24_delta_IGKL)

# statistical analysis of mutation usage between AR2-4 bnabs and representative E2 neg clonotypes
AR24_E2n_filt_sig_IGKL <- data.frame()
for(i in 1:length(AR24_delta_IGKL)){
  sub <- AR24_delta_IGKL[i]
  df <- filter(AR24_E2n_filt_IGKL, Substitution == sub)
  fisher_pvalue <- fisher.test(matrix(c(df$Count[1], df$Count[2], df$Total[1]-df$Count[1], df$Total[2]-df$Count[2]), ncol=2))$p.value
  df_fish <- data.frame(Substitution = sub, Pvalue = fisher_pvalue, AR24_prop = df$Prop[1], E2n_prop = df$Prop[2]) %>%
    mutate(AR24_enrich = AR24_prop > E2n_prop)
  AR24_E2n_filt_sig_IGKL <- rbind.data.frame(AR24_E2n_filt_sig_IGKL, df_fish, stringsAsFactors = FALSE)
}

# fdr correct p-values, filter by adjusted P < 0.05, include Kabat numbering
kabat_frly_cdr12_IGKL <- data.frame(Position = c(27:38, 56:59, 62:65), Kabat = c("27", "27a", "27b", "27c", "27d", "27e", "27f", "28", "29", "30", "31", "32", "51", "52", "52A", "52B", "52C", "52D", "52D", "52E"))

AR24_E2n_filt_sig_corr_IGKL <- AR24_E2n_filt_sig_IGKL %>%
  mutate(Position = as.numeric(gsub(">.", "", Substitution))) %>%
  left_join(kabat_frly_cdr12_IGKL, by = "Position") %>%
  mutate(Pvalue_adj = p.adjust(Pvalue, "fdr")) %>%
  filter(Pvalue_adj < 0.05) %>%
  mutate(FC = AR24_prop/E2n_prop) %>%
  mutate(Position = as.numeric(gsub(">[A-Z]", "", Substitution))) %>%
  mutate(mut = gsub("[0-9][0-9]>", "", Substitution))


#### CDR3
# align front-layer and non-front-layer bnab CDR3s with a representative subset of E2 negative clonotypes 

### front-layer
## filter E2neg clonotypes to be representative of front-layer bnabs (same VJ genes, same SHM range, same CDR3 lengths) 
E2n_frly_filt_IGKL <- C110_IMGT_IGKL %>%
  filter(Group == "E2neg") %>%
  filter(V_gene_IMGT %in% unique(frly_bnab_IGKL$V_gene)) %>%
  filter(J_gene_IMGT %in% unique(frly_bnab_IGKL$J_gene)) %>%
  mutate(CDR3_len = nchar(CDR3_aa_IMGT)) %>%
  filter(CDR3_len >= unique(min(frly_bnab_IGKL$cdr3_len)) & 
           CDR3_len <= unique(max(frly_bnab_IGKL$cdr3_len))) %>%
  filter(V_identity >= unique(min(frly_bnab_IGKL$V_identity)) &
           V_identity <= unique(max(frly_bnab_IGKL$V_identity))) %>%
  select(Group, CDR3_aa = CDR3_aa_IMGT, V_gene_IMGT, J_gene_IMGT, V_identity)

## export table of E2neg clonotypes and front-layer bnabs and run alignment for each on the cluster
  # write.table(E2n_frly_filt_IGKL , "~/Desktop/E2n_frly_filt_IGKL.txt", quote=FALSE , sep = "\t", row.names = FALSE)
  # write.table(frly_bnab_IGKL, "~/Desktop/frly_bnab_IGKL.txt", quote=FALSE , sep = "\t", row.names = FALSE)

## import alignments and combine to one dataframe
frly_cdr3_align_IGKL <- read.csv("~/Desktop/Serial_C110/Final_CDR3_alignments/frly_bnab_align_IGKL.csv") %>%
  mutate(Group = "bnab") %>%
  select(Seq, Group)

E2n_frly_align_IGKL <- read.csv("~/Desktop/Serial_C110/Final_CDR3_alignments/E2n_frly_align_IGKL.csv") %>%
  mutate(Group = "E2neg") %>%
  select(Seq, Group)

E2neg_frly_cdr3_align_IGKL <- rbind.data.frame(frly_cdr3_align_IGKL, E2n_frly_align_IGKL, stringsAsFactors = FALSE)


## Determine amino acid usage
# unable to look for substitutions since no consensus germline CDR3 exists; consequently will look at amino acids used by front-layer bnabs at each position and determine if usage is statistically significantly higher than for E2 neg clonotypes 
frly_cdr3_usage_IGKL <- AA_usage(E2neg_frly_cdr3_align_IGKL[E2neg_frly_cdr3_align_IGKL$Group == "bnab",]$Seq, nchar(E2neg_frly_cdr3_align_IGKL$Seq[1]))
E2n_frlyfilt_cdr3_usage_IGKL <- AA_usage(E2neg_frly_cdr3_align_IGKL[E2neg_frly_cdr3_align_IGKL$Group == "E2neg",]$Seq, nchar(E2neg_frly_cdr3_align_IGKL$Seq[1]))

## filter to only include analysis of amino acids present in front-layer bnabs; then compare usage and obtain P values
E2neg_frly_cdr3_align_stat_pre_IGKL <- data.frame()
for(i in 1:nchar(E2neg_frly_cdr3_align_IGKL$Seq[1])){
  frly_filt <- frly_cdr3_usage_IGKL[,i]
  frly_filt2 <- frly_filt[frly_filt!=0]
  E2n_filt <- E2n_frlyfilt_cdr3_usage_IGKL[,i] 
  E2n_filt2 <- E2n_filt[names(E2n_filt) %in% names(frly_filt2)]

  amino_acids <- names(frly_filt2)
  
  frly_E2n_fish <- data.frame()
  for(j in 1:length(amino_acids)){
    AA <- amino_acids[j]
    AA_count_frly <- frly_filt[names(frly_filt) == AA]
    AA_count_E2n <- E2n_filt[names(E2n_filt) == AA]
    AA_total_frly <- nrow(frly_cdr3_align_IGKL)
    AA_total_E2n <- nrow(E2n_frly_align_IGKL)
    fish_mat <-matrix(c(AA_count_frly, AA_count_E2n, AA_total_frly-AA_count_frly, AA_total_E2n-AA_count_E2n), ncol = 2)
    fish_test <- cbind.data.frame(AA = AA, 
                                  frly_prop = AA_count_frly/AA_total_frly, 
                                  E2n_prop = AA_count_E2n/AA_total_E2n,
                                  P_value = fisher.test(fish_mat)$p.value)  
    frly_E2n_fish <- rbind.data.frame(frly_E2n_fish, fish_test, stringsAsFactors = FALSE)
  }
  
  frly_E2n_fish2 <- frly_E2n_fish %>%
    mutate(Position = i)
  
  E2neg_frly_cdr3_align_stat_pre_IGKL <- rbind.data.frame(E2neg_frly_cdr3_align_stat_pre_IGKL, frly_E2n_fish2, stringsAsFactors = FALSE)
}

## fdr correct p-values, convert to Kabat numbering, find fold-change

# kabat numbering
kabat_frly_cdr3_IGKL <- data.frame(Position = c(1:14), Kabat = c("89", "90", "91", "92", "93", "94", "95", "95a", "95b","95c", "95d", "95e","96", "97"))

# for fold-change calculation, have to replace 0s (replace with  1/2 of the smallest non-zero number)
cdr3_frly_zero_corr_IGKL <- min(
  min(E2neg_frly_cdr3_align_stat_pre_IGKL$frly_prop[E2neg_frly_cdr3_align_stat_pre_IGKL$frly_prop != 0])/2,
  min(E2neg_frly_cdr3_align_stat_pre_IGKL$E2n_prop[E2neg_frly_cdr3_align_stat_pre_IGKL$E2n_prop != 0])/2)

# fdr correction
E2neg_frly_cdr3_align_stat_IGKL <- E2neg_frly_cdr3_align_stat_pre_IGKL %>%
  left_join(kabat_frly_cdr3_IGKL, by= "Position") %>%
  mutate(frly_prop_adj = ifelse(frly_prop == 0, cdr3_frly_zero_corr_IGKL, frly_prop)) %>%
  mutate(E2n_prop_adj = ifelse(E2n_prop == 0, cdr3_frly_zero_corr_IGKL, E2n_prop)) %>%
  mutate(FC = frly_prop_adj/E2n_prop_adj) %>%
  mutate(P_value_adj = p.adjust(P_value, "fdr"))

# filter for significantly differences
E2neg_frly_cdr3_align_sig_IGKL <- E2neg_frly_cdr3_align_stat_IGKL  %>%
  filter(P_value_adj < 0.05) 

### non-front-layer
## AR1
## filter E2neg clonotypes to be representative of AR1 bnabs (same VJ genes, same SHM range, same CDR3 lengths) 
E2n_AR1_filt_IGKL <- C110_IMGT_IGKL %>%
  filter(Group == "E2neg") %>%
  filter(V_gene_IMGT %in% unique(AR1_bnab_IGKL$V_gene)) %>%
  filter(J_gene_IMGT %in% unique(AR1_bnab_IGKL$J_gene)) %>%
  mutate(CDR3_len = nchar(CDR3_aa_IMGT)) %>%
  filter(CDR3_len >= unique(min(AR1_bnab_IGKL$cdr3_len)) & 
           CDR3_len <= unique(max(AR1_bnab_IGKL$cdr3_len))) %>%
  filter(V_identity >= unique(min(AR1_bnab_IGKL$V_identity)) &
           V_identity <= unique(max(AR1_bnab_IGKL$V_identity))) %>%
  select(Group, CDR3_aa = CDR3_aa_IMGT, V_gene_IMGT, J_gene_IMGT, V_identity)

## export table of E2neg clonotypes and AR1 bnabs and run alignment for each on the cluster
  # write.table(E2n_AR1_filt_IGKL , "~/Desktop/E2n_AR1_filt_IGKL.txt", quote=FALSE , sep = "\t", row.names = FALSE)
  # write.table(AR1_bnab_IGKL, "~/Desktop/AR1_bnab_IGKL.txt", quote=FALSE , sep = "\t", row.names = FALSE)

## import alignments and combine to one dataframe
AR1_cdr3_align_IGKL <- read.csv("~/Desktop/Serial_C110/Final_CDR3_alignments/AR1_bnab_align_IGKL.csv") %>%
  mutate(Group = "bnab") %>%
  select(Seq, Group)

E2n_AR1_align_IGKL <- read.csv("~/Desktop/Serial_C110/Final_CDR3_alignments/E2n_AR1_align_IGKL.csv") %>%
  mutate(Group = "E2neg") %>%
  select(Seq, Group)

E2neg_AR1_cdr3_align_IGKL <- rbind.data.frame(AR1_cdr3_align_IGKL, E2n_AR1_align_IGKL, stringsAsFactors = FALSE)


## Determine amino acid usage
# unable to look for substitutions since no consensus germline CDR3 exists; consequently will look at amino acids used by AR1 bnabs at each position and determine if usage is statistically significantly higher than for E2 neg clonotypes 
AR1_cdr3_usage_IGKL <- AA_usage(E2neg_AR1_cdr3_align_IGKL[E2neg_AR1_cdr3_align_IGKL$Group == "bnab",]$Seq, nchar(E2neg_AR1_cdr3_align_IGKL$Seq[1]))
E2n_AR1filt_cdr3_usage_IGKL <- AA_usage(E2neg_AR1_cdr3_align_IGKL[E2neg_AR1_cdr3_align_IGKL$Group == "E2neg",]$Seq, nchar(E2neg_AR1_cdr3_align_IGKL$Seq[1]))

## filter to only include analysis of amino acids present in AR1 bnabs; then compare usage and obtain P values
E2neg_AR1_cdr3_align_stat_pre_IGKL <- data.frame()
for(i in 1:nchar(E2neg_AR1_cdr3_align_IGKL$Seq[1])){
  AR1_filt <- AR1_cdr3_usage_IGKL[,i]
  AR1_filt2 <- AR1_filt[AR1_filt!=0]
  E2n_filt <- E2n_AR1filt_cdr3_usage_IGKL[,i] 
  E2n_filt2 <- E2n_filt[names(E2n_filt) %in% names(AR1_filt2)]

  amino_acids <- names(AR1_filt2)
  
  AR1_E2n_fish <- data.frame()
  for(j in 1:length(amino_acids)){
    AA <- amino_acids[j]
    AA_count_AR1 <- AR1_filt[names(AR1_filt) == AA]
    AA_count_E2n <- E2n_filt[names(E2n_filt) == AA]
    AA_total_AR1 <- nrow(AR1_cdr3_align_IGKL)
    AA_total_E2n <- nrow(E2n_AR1_align_IGKL)
    fish_mat <-matrix(c(AA_count_AR1, AA_count_E2n, AA_total_AR1-AA_count_AR1, AA_total_E2n-AA_count_E2n), ncol = 2)
    fish_test <- cbind.data.frame(AA = AA, 
                                  AR1_prop = AA_count_AR1/AA_total_AR1, 
                                  E2n_prop = AA_count_E2n/AA_total_E2n,
                                  P_value = fisher.test(fish_mat)$p.value)  
    AR1_E2n_fish <- rbind.data.frame(AR1_E2n_fish, fish_test, stringsAsFactors = FALSE)
  }
  
  AR1_E2n_fish2 <- AR1_E2n_fish %>%
    mutate(Position = i)
  
  E2neg_AR1_cdr3_align_stat_pre_IGKL <- rbind.data.frame(E2neg_AR1_cdr3_align_stat_pre_IGKL, AR1_E2n_fish2, stringsAsFactors = FALSE)
}

## fdr correct p-values, convert to Kabat numbering, find fold-change

# kabat numbering
kabat_AR1_cdr3_IGKL <- data.frame(Position = c(1:14), Kabat = c("89", "90", "91", "92", "93", "94", "95", "95a", "95b","95c", "95d", "95e","96", "97"))

# for fold-change calculation, have to replace 0s (replace with  1/2 of the smallest non-zero number)
cdr3_AR1_zero_corr_IGKL <- min(
  min(E2neg_AR1_cdr3_align_stat_pre_IGKL$AR1_prop[E2neg_AR1_cdr3_align_stat_pre_IGKL$AR1_prop != 0])/2,
  min(E2neg_AR1_cdr3_align_stat_pre_IGKL$E2n_prop[E2neg_AR1_cdr3_align_stat_pre_IGKL$E2n_prop != 0])/2)

# fdr correction
E2neg_AR1_cdr3_align_stat_IGKL <- E2neg_AR1_cdr3_align_stat_pre_IGKL %>%
  left_join(kabat_AR1_cdr3_IGKL, by= "Position") %>%
  mutate(AR1_prop_adj = ifelse(AR1_prop == 0, cdr3_AR1_zero_corr_IGKL, AR1_prop)) %>%
  mutate(E2n_prop_adj = ifelse(E2n_prop == 0, cdr3_AR1_zero_corr_IGKL, E2n_prop)) %>%
  mutate(FC = AR1_prop_adj/E2n_prop_adj) %>%
  mutate(P_value_adj = p.adjust(P_value, "fdr"))

# filter for significantly differences
E2neg_AR1_cdr3_align_sig_IGKL <- E2neg_AR1_cdr3_align_stat_IGKL  %>%
  filter(P_value_adj < 0.05) 

## AR2-4
## filter E2neg clonotypes to be representative of AR2-4 bnabs (same VJ genes, same SHM range, same CDR3 lengths) 
E2n_AR24_filt_IGKL <- C110_IMGT_IGKL %>%
  filter(Group == "E2neg") %>%
  filter(V_gene_IMGT %in% unique(AR24_bnab_IGKL$V_gene)) %>%
  filter(J_gene_IMGT %in% unique(AR24_bnab_IGKL$J_gene)) %>%
  mutate(CDR3_len = nchar(CDR3_aa_IMGT)) %>%
  filter(CDR3_len >= unique(min(AR24_bnab_IGKL$cdr3_len)) & 
           CDR3_len <= unique(max(AR24_bnab_IGKL$cdr3_len))) %>%
  filter(V_identity >= unique(min(AR24_bnab_IGKL$V_identity)) &
           V_identity <= unique(max(AR24_bnab_IGKL$V_identity))) %>%
  select(Group, CDR3_aa = CDR3_aa_IMGT, V_gene_IMGT, J_gene_IMGT, V_identity)

## export table of E2neg clonotypes and AR2-4 bnabs and run alignment for each on the cluster
  # write.table(E2n_AR24_filt_IGKL , "~/Desktop/E2n_AR24_filt_IGKL.txt", quote=FALSE , sep = "\t", row.names = FALSE)
  # write.table(AR24_bnab_IGKL, "~/Desktop/AR24_bnab_IGKL.txt", quote=FALSE , sep = "\t", row.names = FALSE)

## import alignments and combine to one dataframe
AR24_cdr3_align_IGKL <- read.csv("~/Desktop/Serial_C110/Final_CDR3_alignments/AR24_bnab_align_IGKL.csv") %>%
  mutate(Group = "bnab") %>%
  select(Seq, Group)

E2n_AR24_align_IGKL <- read.csv("~/Desktop/Serial_C110/Final_CDR3_alignments/E2n_AR24_align_IGKL.csv") %>%
  mutate(Group = "E2neg") %>%
  select(Seq, Group)

E2neg_AR24_cdr3_align_IGKL <- rbind.data.frame(AR24_cdr3_align_IGKL, E2n_AR24_align_IGKL, stringsAsFactors = FALSE)


## Determine amino acid usage
# unable to look for substitutions since no consensus germline CDR3 exists; consequently will look at amino acids used by AR2-4 bnabs at each position and determine if usage is statistically significantly higher than for E2 neg clonotypes 
AR24_cdr3_usage_IGKL <- AA_usage(E2neg_AR24_cdr3_align_IGKL[E2neg_AR24_cdr3_align_IGKL$Group == "bnab",]$Seq, nchar(E2neg_AR24_cdr3_align_IGKL$Seq[1]))
E2n_AR24filt_cdr3_usage_IGKL <- AA_usage(E2neg_AR24_cdr3_align_IGKL[E2neg_AR24_cdr3_align_IGKL$Group == "E2neg",]$Seq, nchar(E2neg_AR24_cdr3_align_IGKL$Seq[1]))

## filter to only include analysis of amino acids present in AR2-4; then compare usage and obtain P values
E2neg_AR24_cdr3_align_stat_pre_IGKL <- data.frame()
for(i in 1:nchar(E2neg_AR24_cdr3_align_IGKL$Seq[1])){
  AR24_filt <- AR24_cdr3_usage_IGKL[,i]
  AR24_filt2 <- AR24_filt[AR24_filt!=0]
  E2n_filt <- E2n_AR24filt_cdr3_usage_IGKL[,i] 
  E2n_filt2 <- E2n_filt[names(E2n_filt) %in% names(AR24_filt2)]

  amino_acids <- names(AR24_filt2)
  
  AR24_E2n_fish <- data.frame()
  for(j in 1:length(amino_acids)){
    AA <- amino_acids[j]
    AA_count_AR24 <- AR24_filt[names(AR24_filt) == AA]
    AA_count_E2n <- E2n_filt[names(E2n_filt) == AA]
    AA_total_AR24 <- nrow(AR24_cdr3_align_IGKL)
    AA_total_E2n <- nrow(E2n_AR24_align_IGKL)
    fish_mat <-matrix(c(AA_count_AR24, AA_count_E2n, AA_total_AR24-AA_count_AR24, AA_total_E2n-AA_count_E2n), ncol = 2)
    fish_test <- cbind.data.frame(AA = AA, 
                                  AR24_prop = AA_count_AR24/AA_total_AR24, 
                                  E2n_prop = AA_count_E2n/AA_total_E2n,
                                  P_value = fisher.test(fish_mat)$p.value)  
    AR24_E2n_fish <- rbind.data.frame(AR24_E2n_fish, fish_test, stringsAsFactors = FALSE)
  }
  
  AR24_E2n_fish2 <- AR24_E2n_fish %>%
    mutate(Position = i)
  
  E2neg_AR24_cdr3_align_stat_pre_IGKL <- rbind.data.frame(E2neg_AR24_cdr3_align_stat_pre_IGKL, AR24_E2n_fish2, stringsAsFactors = FALSE)
}

## fdr correct p-values, convert to Kabat numbering, find fold-change

# kabat numbering
kabat_AR24_cdr3_IGKL <- data.frame(Position = c(1:12), Kabat = c("89", "90", "91", "92", "93", "94", "95", "95a", "95b","95c", "96", "97"))

# for fold-change calculation, have to replace 0s (replace with  1/2 of the smallest non-zero number)
cdr3_AR24_zero_corr_IGKL <- min(
  min(E2neg_AR24_cdr3_align_stat_pre_IGKL$AR24_prop[E2neg_AR24_cdr3_align_stat_pre_IGKL$AR24_prop != 0])/2,
  min(E2neg_AR24_cdr3_align_stat_pre_IGKL$E2n_prop[E2neg_AR24_cdr3_align_stat_pre_IGKL$E2n_prop != 0])/2)

# fdr correction
E2neg_AR24_cdr3_align_stat_IGKL <- E2neg_AR24_cdr3_align_stat_pre_IGKL %>%
  left_join(kabat_AR24_cdr3_IGKL, by= "Position") %>%
  mutate(AR24_prop_adj = ifelse(AR24_prop == 0, cdr3_AR24_zero_corr_IGKL, AR24_prop)) %>%
  mutate(E2n_prop_adj = ifelse(E2n_prop == 0, cdr3_AR24_zero_corr_IGKL, E2n_prop)) %>%
  mutate(FC = AR24_prop_adj/E2n_prop_adj) %>%
  mutate(P_value_adj = p.adjust(P_value, "fdr"))

# filter for significantly differences
E2neg_AR24_cdr3_align_sig_IGKL <- E2neg_AR24_cdr3_align_stat_IGKL  %>%
  filter(P_value_adj < 0.05) 
  
  
  
  
  

    