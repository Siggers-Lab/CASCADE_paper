
### DESCRIPTION -----------------------------------------------------------

# script that generates figures used in CASCADE paper based on supplementary data

### F1 - CASCADE SCHEMATIC --------------------------------------------------------------------

# NO DATA PROCESSING OR FIGURES GENERATED FOR THIS SCHEMATIC

### F2 - CXCL10 CASCADE -------------------------------------------------------------------------

# set up environment
rm(list=ls())

# set up environment
rm(list=ls())
library(ggplot2)
library(ggseqlogo)
library(RColorBrewer)
library(cowplot)

# function to obtain reverse complement PWM
rev_compl_pwm <- function(m) {
  # reverse the input pwm
  rev_pwm <- m[,ncol(m):1]
  
  # reverse complement the pwm (A to T, C to G)
  rc_pwm <- rev_pwm
  rc_pwm[1,] <- rev_pwm[4,]
  rc_pwm[2,] <- rev_pwm[3,]
  rc_pwm[3,] <- rev_pwm[2,]
  rc_pwm[4,] <- rev_pwm[1,]
  rownames(rc_pwm) <- rownames(m)
  
  # return reversed pwm
  return(rc_pwm)
}

# function to obtain reverse complement of a DNA sequence
rev_compl <- function(DNA_seq){
  # returns the input DNA sequence in reverse
  temp_seq <- unlist(strsplit(DNA_seq, split=''))
  temp_seq <- rev(temp_seq)
  temp_seq <- paste(temp_seq, collapse='')
  temp_seq <- chartr("ACGTN", "TGCAN", temp_seq)
  return(temp_seq)
}

# function to obtain SNV value matrix (using seed-weighted mean - gain-of-function variants are also flagged and reverted to seed z-score value)
locus_SNV_matrix <- function(PBM, PBM_col, genomic_window) {
  
  # initialize a matrix to hold the SNV probe z-score values
  nucs <- c("A", "C", "G", "T")
  SNV <- matrix(nrow=4, ncol=length(genomic_window))
  
  # collect fluo (treatment/UT) values for each position in the matrix
  for (i in 1:nrow(SNV)) {
    for (j in 1:ncol(SNV)) {
      
      # determine the current genomic position
      genomic_pos <- genomic_window[j]
      
      # collect data relevant to the SNV probes at the current genomic position
      SNV_probes <- PBM[which(PBM$SNV_pos==genomic_pos),]
      
      # determine the current genomic nucleotide
      curr_seed_nuc <- unique(as.character(SNV_probes$seed_nuc))
      
      # if the current SNV nucleotide is the genomic nucleotide
      if (curr_seed_nuc == nucs[i]) {
        
        # collect tile numbers that contain the given SNV position
        tiles <- unique(SNV_probes$tile_order)
        
        # collect info about the tile probes that contain the current nucleotide at the SNV position
        curr_nuc_TILE_probes <- PBM[which(PBM$tile_order %in% tiles & PBM$SNV_pos_offset==0),]
        
        # determine weights of genomic tile probes to use for the weighted average
        weights <- curr_nuc_TILE_probes[, PBM_col]
        
        # save the z-scores of the current TILE probes for use in the weighted mean formula (in a way that preserves tile order used in prev step)
        probe_z_scores <- curr_nuc_TILE_probes[, PBM_col]
        
        # replace non-finite values with 0
        weights[!is.finite(weights)] <- 0
        probe_z_scores[!is.finite(probe_z_scores)] <- 0
        
        # replace non-positive values with a pseudo-count
        weights[weights<=0] <- 0.000001
        probe_z_scores[probe_z_scores<=0] <- 0.000001
        
        # assign the SNV value as a function of the PBM values for the genomic tiles containing the given SNV position
        SNV[i,j] <- weighted.mean(probe_z_scores, weights)
        
      } else { # if the current nucleotide corresponds to a SNV of the genomic seed nucleotide
        
        # collect tile numbers that contain the given SNV position
        tiles <- unique(SNV_probes$tile_order)
        
        # collect info about the tile probes that contain the current nucleotide at the SNV position
        curr_nuc_TILE_probes <- PBM[which(PBM$tile_order %in% tiles & PBM$SNV_pos_offset==0),]
        
        # collect info about the SNV probes that contain the current nucleotide at the SNV position
        curr_nuc_SNV_probes <- SNV_probes[which(SNV_probes$SNV_nuc==nucs[i]),]
        
        # determine weights of genomic tile probes to use
        weights <- vector(mode="numeric")
        for (k in 1:nrow(curr_nuc_SNV_probes)) {
          weight_k <- curr_nuc_TILE_probes[which(curr_nuc_TILE_probes$tile_order==curr_nuc_SNV_probes[k, "tile_order"]), PBM_col]
          weights <- c(weights, weight_k)
        }
        
        # save the z-scores of the current SNV probes for use in the weighted mean formula (in a way that preserves tile order used in prev step)
        probe_z_scores <- curr_nuc_SNV_probes[, PBM_col]
        
        # replace non-finite values with 0
        weights[!is.finite(weights)] <- 0
        probe_z_scores[!is.finite(probe_z_scores)] <- 0
        
        # replace non-positive values with a pseudo-count
        weights[weights<=0] <- 0.000001
        probe_z_scores[probe_z_scores<=0] <- 0.000001
        
        # modify the z-scores if the seed value is below threshold
        for (k in 1:length(probe_z_scores)) {
          
          # if the variant percentile is greater than 0.95 and the seed percentile is less than 0.95
          if (probe_z_scores[k] > 1.645 & weights[k] <= 1.645) {
            
            # nullify effects of the variant by setting it as the seed value
            probe_z_scores[k] <- weights[k]
            
          }
        }
        
        # assign the SNV value as a function of the PBM values for the genomic tiles containing the given SNV position
        SNV[i,j] <- weighted.mean(probe_z_scores, weights)
        
      }
    }
  }
  
  # transform fluo value matrix prior to returning to parent function to plot
  SNV_trans <- matrix(nrow=4, ncol=ncol(SNV))
  rownames(SNV_trans) <- c("A", "C", "G", "T")
  for (i in 1:nrow(SNV_trans)) {
    for (j in 1:ncol(SNV_trans)) {
      # subtract the positional column-wise median if looking at condition-specific z-scores relative to the background set of probes
      SNV_trans[i, j] <- SNV[i, j] - median(SNV[,j])
    }
  }
  
  # reverse complement the pwm if necessary
  if (dir == "REVERSE") {
    SNV_trans <- rev_compl_pwm(SNV_trans)
  }
  
  # return the matrix to be plotted as a sequence logo
  return(SNV_trans)
}

# function to transfrom a SNV energy matrix to PWM
SNV_to_PWM <- function(SNV, beta) {
  
  # transform using beta parameter
  SNV_pwm <- exp(beta * SNV)
  
  # compute probabilities for PWM using the column-wise totals
  PWM <- matrix(nrow=4, ncol=ncol(SNV_pwm))
  for (i in 1:nrow(PWM)) {
    for (j in 1:ncol(PWM)) {
      PWM[i, j] <- SNV_pwm[i, j]/sum(SNV_pwm[,j])
    }
  }
  rownames(PWM) <- c("A", "C", "G", "T")
  
  # return the PWM to parent function
  return(PWM)
  
}

# scale transformation function for rounding axes
scaleFUN <- function(x) sprintf("%.1f", x)

# create directory to compartmentalize the outputs from this block
dir.create(paste(getwd(), "Fig2_CXCL10_CASCADE", sep="/"))

# create directory to house the CASCADE logo plots
dir.create(paste(getwd(), "Fig2_CXCL10_CASCADE", "CXCL10_CASCADE_logos", sep="/"))

# create directory to house the motif matrices for PWM comparison (via TOMTOM)
dir.create(paste(getwd(), "Fig2_CXCL10_CASCADE", "CXCL10_CASCADE_PWMs", sep="/"))

# read in the CXCL10 CASCADE data
df <- read.table("suppl_data_1_CXCL10_CASCADE.dat", header=T, sep='\t', stringsAsFactors = F)

# read in the known TF binding site annotation
TF_sites <- read.table("CXCL10_TF_sites.txt", header=T, sep='\t')
TF_sites$hex_color <- paste("#", TF_sites$hex_color, sep="")

# declare seeds and antibodies used
seeds <- c("CXCL10")
antibodies <- c("IRF2",
                "IRF3",
                "IRF8",
                "p300",
                "p65",
                "RBBP5")

# for each seed sequence
for (seed in seeds) {
  
  # keep only data that matches the seed
  x <- droplevels(df[which(df$seed_names==seed),])
  
  # grab the actual genomic window (nucleotide positions within a given chromosome)
  genomic_window <- min(x$tile_starts):max(x$tile_ends)
  
  # construct the genomic sequence of the given locus
  genomic_seq <- unique(x[,c("SNV_pos", "seed_nuc")]) # collect the genomic nucleotides for each genomic position
  genomic_seq <- genomic_seq[order(genomic_seq$SNV_pos),] # ensure SNV positions are sorted
  genomic_seq <- as.character(genomic_seq[which(!is.na(genomic_seq$seed_nuc)),"seed_nuc"]) # get rid of the NA from the tile probe
  genomic_seq <- paste(genomic_seq, collapse="") # collapse the character vector into a single character string
  
  # restrict TF binding site annotation to given seed
  seed_TF_sites <- droplevels(TF_sites[which(TF_sites$seed_names==seed),])
  
  # declare the possible strand directions (FORWARD is 5' to 3' on the + strand, REVERSE is 5' to 3' on the - strand)
  directions <- c(#"FORWARD",
                  "REVERSE")
  
  # for both possible directions (FORWARD is 5' to 3' on the + strand, REVERSE is 5' to 3' on the - strand)
  for (dir in directions) {
    
    # reverse complement the current genomic sequence and reverse the window if REVERSE
    curr_genomic_seq <- genomic_seq
    curr_genomic_window <- genomic_window
    if (dir == "REVERSE") {
      curr_genomic_seq <- rev_compl(curr_genomic_seq)
      curr_genomic_window <- rev(genomic_window)
    }
    
    ### GENOMIC LADDER
    
    # plot genomic ladder
    seq_info <- data.frame(x = 1:length(genomic_window),
                           y = 1,
                           nuc = strsplit(curr_genomic_seq, "")[[1]])
    seq_info$nuc <- as.character(seq_info$nuc)
    ggp_ladder <- ggplot(seq_info, aes(x, y, label=nuc)) +
      geom_text() +
      labs(title = paste("hg38: ", as.character(unique(x$tile_chroms, sep="")))) +
      theme(plot.title = element_text(lineheight=0.8, face="bold", hjust=0)) +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
      theme(axis.line=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      # extra space in limits ensures perfect alignment with motif plotted in ggp_2 below
      scale_x_continuous(limits=c(0.5, length(genomic_window)+0.5), expand = c(0.01,0.01), position = "top", breaks = which(curr_genomic_window %% 10 == 0), labels = curr_genomic_window[which(curr_genomic_window %% 10 == 0)]) +
      theme(axis.text.x=element_text(hjust=0.15)) +
      scale_y_continuous(limits=c(-2.5, 2.5), expand = c(0.01,0.01)) +
      # remove grey background which is present for some reason
      theme(panel.background = element_blank())
    
    # plot known TF sites onto the genomic ladder
    for (i in 1:nrow(seed_TF_sites)) {
      min_x <- which(curr_genomic_window==(seed_TF_sites[i,"start"]-1))-0.5 # determines the start of the current TF binding site window
      max_x <- which(curr_genomic_window==(seed_TF_sites[i,"end"]+1))+0.5 # determines the end of the current TF binding site window
      
      # change the min_x and max_x if looking at the FORWARD logo
      if (dir == "FORWARD") {
        min_x <- which(curr_genomic_window==(seed_TF_sites[i,"start"]))-0.5 # determines the start of the current TF binding site window
        max_x <- which(curr_genomic_window==(seed_TF_sites[i,"end"]))+0.5 # determines the end of the current TF binding site window
      }
      
      # annotate the known TF site by highlighting it with a semi-transparent box
      ggp_ladder <- ggp_ladder + annotate('rect', xmin = min_x, xmax = max_x, ymin=0, ymax=2, fill=seed_TF_sites[i,"hex_color"], alpha=0.35)
      # label the current TF site using the name listed in the file (centered beneath the locus)
      ggp_ladder <- ggp_ladder + annotate('text', x = (min_x + max_x)/2, y = -1.5, label = seed_TF_sites[i,"TF_names"])
    }
    
    
    
    # for each unique antibody profiled
    for (ab in antibodies) {
      
      # determine dataframe conditions that contain the given antibody
      cond <- grep(ab, names(df)[14:length(names(df))])
      cond <- names(df)[14:length(names(df))][cond]
      
      # initialize vectors needed for future pwm_df
      z_stack <- vector(mode="numeric")
      treatment <- vector(mode="character")
      seed_pos <- vector(mode="numeric")
      
      # initialize variables needed to keep track of MAX and MIN delta z-scores
      z_min <- 999
      z_max <- -999
      
      # collect the MAX and MIN delta z-scores across all experiments for the given antibody
      for (z in 1:length(cond)) {
        
        # build the given SNV energy matrix
        pwm <- locus_SNV_matrix(x, cond[z], genomic_window)
        
        # generate copies to modify for positive and negative colsums
        pos_pwm <- pwm
        neg_pwm <- pwm
        
        # replace either the positive or negative values with 0 in the opposite pwm
        pos_pwm[pos_pwm<0] <- 0
        neg_pwm[neg_pwm>0] <- 0
        
        # get the colSums to take a cumulative stack height (for just the positive values)
        pos_pwm <- colSums(pos_pwm)
        neg_pwm <- colSums(neg_pwm)
        
        # collect the maximum and minimum stack heights across the matrix
        curr_z_min <- min(neg_pwm)
        curr_z_max <- max(pos_pwm)
        
        # remove pwm just in case
        rm(pwm, pos_pwm, neg_pwm)
        
        # update the current maximum and minimum if necessary
        if (curr_z_min < z_min) {
          z_min <- curr_z_min
        }
        if (curr_z_max > z_max) {
          z_max <- curr_z_max
        }
        
      }
      
      # use the number of conditions to determine which colors to use for the plot
      if (length(cond)<=8) {

        # there are enough colors in the default palette to plot
        PBM_colors <- brewer.pal(8, "Set2")
        
        
      } else {
        
        # construct a palette with the correct number of colors (equal to the number of seeds to be plotted) using the Set2 as a reference
        PBM_colors <- colorRampPalette(brewer.pal(8, "Set2"))
        PBM_colors <- PBM_colors(length(cond))
        
      }
      
      # obtain the plots for each condition
      for (z in 1:length(cond)) {
        
        ### SNV ENERGY LOGO 
        
        # build the given SNV energy matrix
        pwm <- locus_SNV_matrix(x, cond[z], genomic_window)
        
        # replace 0s and NaNs with pseudocount to prevent plotting errors
        pwm[pwm==0] <- 0.000001
        pwm[is.nan(pwm)] <- 0.0000001
        pwm[is.na(pwm)] <- 0.0000001
        
        # determine beta parameter value for the given experiment
        beta <- 30/max(x[, cond[z]])
        
        # save a temp version of the snv matrix
        pwm_temp <- pwm
        snv <- pwm
        pwm <- SNV_to_PWM(snv, beta)
        
        # separate into 3 defined regions
        pwm_ISRE <- pwm[ ,1:65]
        pwm_NFKB2 <- pwm[ ,66:110]
        pwm_NFKB1 <- pwm[ ,111:ncol(pwm)]
        
        # save each as a separate motif (including whole segment)
        mat_name <- paste(cond[z], "PWM.txt", sep="_")
        write.table(pwm, paste(getwd(), "Fig2_CXCL10_CASCADE", "CXCL10_CASCADE_PWMs", mat_name, sep="/"), quote=F, row.names=T, col.names=F, sep='\t')
        
        mat_name <- paste(cond[z], "ISRE", "PWM.txt", sep="_")
        write.table(pwm_ISRE, paste(getwd(), "Fig2_CXCL10_CASCADE", "CXCL10_CASCADE_PWMs", mat_name, sep="/"), quote=F, row.names=T, col.names=F, sep='\t')
        
        mat_name <- paste(cond[z], "NFKB2", "PWM.txt", sep="_")
        write.table(pwm_NFKB2, paste(getwd(), "Fig2_CXCL10_CASCADE", "CXCL10_CASCADE_PWMs", mat_name, sep="/"), quote=F, row.names=T, col.names=F, sep='\t')
        
        mat_name <- paste(cond[z], "NFKB1", "PWM.txt", sep="_")
        write.table(pwm_NFKB1, paste(getwd(), "Fig2_CXCL10_CASCADE", "CXCL10_CASCADE_PWMs", mat_name, sep="/"), quote=F, row.names=T, col.names=F, sep='\t')
        
        # clean up and restore
        pwm <- pwm_temp
        rm(pwm_temp, snv)
        
        # assign the y axis label based on the unit being plotted
        y_axis_label <- expression(paste(Delta, "z-score", sep=""))
        
        # plot the SNV energy logo
        ggp_temp <- ggseqlogo(pwm, method="custom")
        my.ggp.yrange <- ggplot_build(ggp_temp)$layout$panel_params[[1]]$y.range
        ggp_temp <- ggp_temp +
          scale_y_continuous(expand = c(0.01,0.01), labels=scaleFUN, limits=c(z_min-0.1, z_max+0.1)) +
          scale_x_continuous(expand = c(0.01,0.01)) +
          annotate('rect', xmin = 0.5, xmax = ncol(pwm)+0.5, ymin = z_min, ymax = 0.0, alpha = 0.75, col='white', fill='white') +
          annotate('rect', xmin = 0.5, xmax = ncol(pwm)+0.5, ymin = 0, ymax = z_max, col="lightgray", fill=NA, size=1.1) +
          annotate('rect', xmin = 0.5, xmax = ncol(pwm)+0.5, ymin = z_min, ymax = 0, col="lightgray", fill=NA, size=1.1) +
          ylab(y_axis_label) +
          theme(axis.text.x=element_blank()) +
          theme(axis.text.y=element_text(size=10, face="plain"), axis.title.y=element_text(size=10, face="plain")) +
          ggtitle(cond[z])
        
        # save the plot under a unique name based on the current iteration index
        assign(paste("ggp_", z, sep=""), ggp_temp)
        rm(ggp_temp)
        
        
        
        ### UPDATE VALUES NEEDED FOR LOGO MAGNITUDE COMPARISON
        
        # TEST: replace the negative values in the SNV pwm with a positive pseudocount
        pwm[pwm<=0] <- 0.000001
        
        # get the colSums to take a cumulative stack height (for just the positive values)
        colsums_pwm <- colSums(pwm)
        
        # update values in the growing vector
        z_stack <- c(z_stack, colsums_pwm)
        treatment <- c(treatment, rep(cond[z], length(colsums_pwm)))
        seed_pos <- c(seed_pos, 1:length(colsums_pwm))
        
      }
      
      
      
      ### LOGO MAGNITUDE COMPARISON
      
      # construct a data frame using the vectors
      pwm_df <- data.frame(seed_pos, treatment, z_stack)
      
      # assign the y axis label based on the unit being plotted
      y_axis_label <- expression(paste("stack height (", Delta, "z-score)", sep=""))
      
      # use ggplot to plot a scatter of the z-score vs. the position
      ggp_pwm <- ggplot(pwm_df, aes(x=seed_pos, y=z_stack, fill=treatment, color=treatment)) + geom_area(position="identity", alpha=0.35) + geom_line(size=1)
      my.ggp.yrange <- ggplot_build(ggp_pwm)$layout$panel_params[[1]]$y.range
      ggp_pwm <- ggp_pwm +
        scale_x_continuous(limits=c(0.5, length(genomic_window)+0.5), expand = c(0.01,0.01), name="") +
        scale_y_continuous(expand = c(0.01, 0.01), labels=scaleFUN) +
        annotate('rect', xmin = 0.5, xmax = length(genomic_window)+0.5, ymin = 0, ymax = my.ggp.yrange[2], col="lightgray", fill=NA, size=1.1) +
        theme(axis.line=element_blank()) +
        ylab(y_axis_label) +
        theme(legend.position="bottom") +
        theme(axis.text.x=element_blank(), axis.ticks=element_blank()) +
        theme(axis.text.y=element_text(size=10, face="plain"), axis.title.y=element_text(size=10, face="plain")) +
        scale_color_manual(values=PBM_colors) +
        scale_fill_manual(values=PBM_colors) +
        guides(fill=guide_legend(nrow=2, byrow=TRUE)) + 
        theme(legend.title = element_blank()) +
        # remove grey background which is supposed to be gone because of cowplot
        theme(panel.background = element_blank())
      
      # compile a plot vector in the order that they should be plotted
      plot_order <- paste0("ggp_", 1:length(cond))
      
      # append the genomic ladder to the top and the logo profile comparison to the bottom
      plot_order <- c("ggp_ladder", plot_order, "ggp_pwm")
      
      ### USE COWPLOT TO SAVE THE PLOTS
      ggp <- plot_grid(plotlist=mget(plot_order),
                       ncol = 1, nrow = length(plot_order), rel_heights = c(1, rep(2.5, length(plot_order)-2), 3), align = "v")
      
      # save PDF
      pdf_name <- paste(seed, ab, dir, "CASCADE.pdf", sep="_")
      save_plot(paste("Fig2_CXCL10_CASCADE", "CXCL10_CASCADE_logos", pdf_name, sep="/"), ggp, ncol = 1, nrow = length(plot_order), base_height = 2.5, base_width = length(genomic_window)/7, base_aspect_ratio = 1)
      
      # save SVG
      svg_name <- paste(seed, ab, dir, "CASCADE.svg", sep="_")
      save_plot(paste("Fig2_CXCL10_CASCADE", "CXCL10_CASCADE_logos", svg_name, sep="/"), ggp, ncol = 1, nrow = length(plot_order), base_height = 2.5, base_width = length(genomic_window)/7, base_aspect_ratio = 1)
      
    }
  }
}

### PREP A MEME-FORMATTED FILE AND AN .SH SCRIPT TO EXECUTE MOTIF ANALYSIS

# collect vector of PWM.txt files
motif_files <- list.files(path=paste(getwd(), "Fig2_CXCL10_CASCADE", "CXCL10_CASCADE_PWMs", sep="/"), pattern="_PWM.txt", full.names = T)

# if the vector isn't empty
if (length(motif_files)>0) {
  
  # intialize a MEME-formatted motif file to write the reformatted motifs to
  meme_file_name <- paste(getwd(), "Fig2_CXCL10_CASCADE", "CXCL10_CASCADE_PWMs", "CXCL10_CASCADE_motifs.meme", sep="/")
  file.create(meme_file_name)
  meme_file <- file(meme_file_name)
  writeLines(c("MEME version 5.0.3\n",
               "\n",
               "ALPHABET= ACGT\n",
               "\n",
               "strands: + -\n",
               "Background letter frequencies\n",
               "A 0.250 C 0.250 G 0.250 T 0.250 \n",
               "\n"),
             con=meme_file, sep="")
  
  # for each motif in the motif_files vector
  for (m_file in motif_files) {
    
    # read the given file
    m <- as.matrix(read.table(m_file, row.names=1))
    
    # transpose the matrix
    m <- t(m)
    
    # create a header for the given motif
    m_header <- unlist(strsplit(m_file, split="/"))[length(unlist(strsplit(m_file, split="/")))]
    m_header <- gsub("_PWM.txt", "", m_header)
    
    # create a subheader for the given motif
    m_subheader <- paste("letter-probability",
                         "matrix:",
                         "alength= 4",
                         "w=", as.character(nrow(m)),
                         "nsites= 0",
                         "E= 0")
    
    # write the header and matrix to file
    write(paste("MOTIF", m_header, sep=" "), meme_file_name, append=T)
    write(m_subheader, meme_file_name, append=T)
    write(formatC(t(m), format="f", 6), meme_file_name, ncolumns=4, sep=" ", append=T)
    write("", meme_file_name, append=T)
  }
  close(meme_file)
  
}



### F3 - SNP SCREEN - SUMMARY VOLCANO PLOTS ----------------------------------------------------------------------

# set up environment
rm(list=ls())
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# read in the SNP screen table for plotting
SNP_screen <- read.table(file="suppl_data_2_SNP_screen_COF_stats.dat", header=T, sep='\t', stringsAsFactors = F)

# read in the PBM conditions used for the SNP screen
PBM_cond <- read.table(file="SNP_screen_COF_exp.txt", stringsAsFactors = F)
PBM_cond <- PBM_cond$V1

# create directory to compartmentalize the outputs from this block
dir.create(paste(getwd(), "Fig3_SNP_screen_COF", sep="/"))

# omit random background category from plots
SNP_screen <- SNP_screen[which(SNP_screen$SNP_cat != "bg_SNP"), ]

# determine the PBM_colors based on the SNP categories
PBM_colors <- brewer.pal(6, "Set2")
PBM_colors <- PBM_colors[c(1, 3:6)] # 2 was previously used to inspect background category

# modify the levels of SNP_cat
SNP_screen$SNP_cat <- factor(SNP_screen$SNP_cat, levels = c("basal_eQTL",
                                                            "caQTL_eQTL",
                                                            "GWAS_caQTL",
                                                            "GWAS_eQTL",
                                                            "resp_eQTL"))

# scale transformation function for rounding axes
scaleFUN <- function(x) sprintf("%.1f", x)

# for each of the PBM experimental conditions tested in this set of experiments
for (i in 1:length(PBM_cond)) {
  
  # subset the SNP screen to only include data from the relevant condition
  x <- SNP_screen[which(SNP_screen$PBM_exp==PBM_cond[i]),]
  
  # determine an approximate pval cut-off equivalent to q < 0.05 for the given experiment
  thresh_pass <- min(x[which(x$qval<0.05), "neg_log10_pval"])
  thresh_fail <- max(x[which(x$qval>=0.05), "neg_log10_pval"])
  approx_cutoff <- mean(c(thresh_pass, thresh_fail))
  
  # declare a plot variable for the current scatter. Needs to be done this way because of cowplot
  ggp <- ggplot(x, aes(y=neg_log10_pval, x=delta, color=SNP_cat)) + 
    geom_point(size=1.5) + scale_color_manual(values=PBM_colors) +
    scale_x_continuous(labels=scaleFUN) + scale_y_continuous(labels=scaleFUN) +
    ggtitle(PBM_cond[i]) +
    ylab("-log10(Fisher's combined p-value)") + xlab(expression(paste(Delta, "z-score", sep=""))) +
    geom_hline(yintercept=approx_cutoff, linetype="dashed", color="grey") +
    theme_cowplot(font_size = 14)
  
  # assign as a separate variable
  assign(paste("ggp_", i, sep=""), ggp)
  
  # remove temp intermediates
  rm(x, thresh_pass, thresh_fail, approx_cutoff)
  
}

# use cowplot to grid the volcano plots
ggp_grid <- plot_grid(plotlist=mget(paste0("ggp_", 1:length(PBM_cond))), ncol = 2, nrow = length(PBM_cond)/2, label_size = 8)

# save a pdf version
pdf_name <- paste("SNP_screen_COF_volcano_plots.pdf", sep="_")
save_plot(paste(getwd(), "Fig3_SNP_screen_COF", pdf_name, sep="/"), ggp_grid, ncol = 2, nrow = length(PBM_cond)/2, base_height=4, base_aspect_ratio=1.4, limitsize=FALSE)

# save an svg version
svg_name <- paste("SNP_screen_COF_volcano_plots.svg", sep="_")
save_plot(paste(getwd(), "Fig3_SNP_screen_COF", svg_name, sep="/"), ggp_grid, ncol = 2, nrow = length(PBM_cond)/2, base_height=4, base_aspect_ratio=1.4, limitsize=FALSE)

# remove plots before continuing
rm(list=paste0("ggp_", 1:length(PBM_cond)))



### F3 - SNP SCREEN - PAIRWISE QVAL SCATTERS (labels) ----------------------------------------------------

# set up environment
rm(list=ls())
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# create a directory to house the orientation scatterplots
dir.create(paste(getwd(), "Fig3_SNP_screen_COF", "pairwise_qval_scatters_LABELS", sep="/"))

# read in the SNP screen table for plotting
SNP_screen <- read.table(file="suppl_data_2_SNP_screen_COF_stats.dat", header=T, sep='\t', stringsAsFactors = F)

# read in the PBM conditions used for the SNP screen
PBM_cond <- read.table(file="SNP_screen_COF_exp.txt", stringsAsFactors = F)
PBM_cond <- PBM_cond$V1

# omit random background category from plots
SNP_screen <- SNP_screen[which(SNP_screen$SNP_cat != "bg_SNP"), ]

# determine the PBM_colors based on the SNP categories
PBM_colors <- brewer.pal(6, "Set2")
PBM_colors <- PBM_colors[c(1, 3:6)]

# scale transformation function for rounding axes
scaleFUN <- function(x) sprintf("%.1f", x)

# plot the pairwise scatters across all PBM conditions (compared against all PBM conditions in a grid)
for (i in 1:length(PBM_cond)) {
  for (j in 1:length(PBM_cond)) {
    
    # create a temp dataframe for the current comparison
    g_df <- data.frame(SNP_screen[which(SNP_screen$PBM_exp==PBM_cond[i]), "SNP_id"],
                       SNP_screen[which(SNP_screen$PBM_exp==PBM_cond[i]), "SNP_cat"],
                       SNP_screen[which(SNP_screen$PBM_exp==PBM_cond[i]), "neg_log10_qval"],
                       SNP_screen[which(SNP_screen$PBM_exp==PBM_cond[j]), "neg_log10_qval"])
    names(g_df) <- c("SNP_id", "SNP_cat", "yval", "xval")
    
    # manually determine the legend order and associated colors for the SNP_cat
    g_df$SNP_cat <- factor(g_df$SNP_cat, levels = c("basal_eQTL",
                                                    "caQTL_eQTL",
                                                    "GWAS_caQTL",
                                                    "GWAS_eQTL",
                                                    "resp_eQTL"))
    
    # declare a plot variable for the current scatter. Needs to be done this way because of cowplot
    ggp <- ggplot(g_df, aes(y=yval, x=xval, color=SNP_cat)) +
      geom_point(size=1.5) + scale_color_manual(values=PBM_colors) +
      scale_x_continuous(labels=scaleFUN) + scale_y_continuous(labels=scaleFUN) +
      ggtitle(paste(PBM_cond[i], " vs. ", PBM_cond[j], sep="")) +
      ylab(paste(PBM_cond[i], " -log10(q)", sep="")) + xlab(paste(PBM_cond[j], " -log10(q)", sep="")) +
      theme(legend.position="right") +
      geom_text(aes(label=ifelse(yval>1.30103 & xval>1.30103, as.character(SNP_id), "")), hjust=1, vjust=1, size=3) +
      geom_hline(yintercept=1.30103, linetype="dashed", color = "grey") +
      geom_vline(xintercept=1.30103, linetype="dashed", color = "grey") +
      theme_cowplot(font_size = 14)
    
    # assign as a separate variable
    assign(paste("ggp_", j, sep=""), ggp)
    
    # use cowplot to grid the scatterplots
    ggp_grid <- plot_grid(plotlist=mget(paste("ggp_", as.character(j), sep="")), ncol = 1, nrow = 1, label_size = 8)
    
    # save pdf version
    pdf_name <- paste("pairwise", PBM_cond[i], PBM_cond[j], "scatter.pdf", sep="_")
    save_plot(paste(getwd(), "Fig3_SNP_screen_COF", "pairwise_qval_scatters_LABELS", pdf_name, sep="/"), ggp_grid, ncol = 1, nrow = 1, base_height=4, base_aspect_ratio=1.4, limitsize=FALSE)

    # remove g_df to prevent overwriting
    rm(g_df, ggp)
    
  }
}



### F3 - SNP SCREEN - PAIRWISE QVAL SCATTERS (NO labels) ----------------------------------------------------

# set up environment
rm(list=ls())
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# create a directory to house the orientation scatterplots
dir.create(paste(getwd(), "Fig3_SNP_screen_COF", "pairwise_qval_scatters_NO_LABELS", sep="/"))

# read in the SNP screen table for plotting
SNP_screen <- read.table(file="suppl_data_2_SNP_screen_COF_stats.dat", header=T, sep='\t', stringsAsFactors = F)

# read in the PBM conditions used for the SNP screen
PBM_cond <- read.table(file="SNP_screen_COF_exp.txt", stringsAsFactors = F)
PBM_cond <- PBM_cond$V1

# omit random background category from plots
SNP_screen <- SNP_screen[which(SNP_screen$SNP_cat != "bg_SNP"), ]

# determine the PBM_colors based on the SNP categories
PBM_colors <- brewer.pal(6, "Set2")
PBM_colors <- PBM_colors[c(1, 3:6)]

# scale transformation function for rounding axes
scaleFUN <- function(x) sprintf("%.1f", x)

# plot the pairwise scatters across all PBM conditions (compared against all PBM conditions in a grid)
for (i in 1:length(PBM_cond)) {
  for (j in 1:length(PBM_cond)) {
    
    # create a temp dataframe for the current comparison
    g_df <- data.frame(SNP_screen[which(SNP_screen$PBM_exp==PBM_cond[i]), "SNP_id"],
                       SNP_screen[which(SNP_screen$PBM_exp==PBM_cond[i]), "SNP_cat"],
                       SNP_screen[which(SNP_screen$PBM_exp==PBM_cond[i]), "neg_log10_qval"],
                       SNP_screen[which(SNP_screen$PBM_exp==PBM_cond[j]), "neg_log10_qval"])
    names(g_df) <- c("SNP_id", "SNP_cat", "yval", "xval")
    
    # manually determine the legend order and associated colors for the SNP_cat
    g_df$SNP_cat <- factor(g_df$SNP_cat, levels = c("basal_eQTL",
                                                    "caQTL_eQTL",
                                                    "GWAS_caQTL",
                                                    "GWAS_eQTL",
                                                    "resp_eQTL"))
    
    # declare a plot variable for the current scatter. Needs to be done this way because of cowplot
    ggp <- ggplot(g_df, aes(y=yval, x=xval, color=SNP_cat)) +
      geom_point(size=1.5) + scale_color_manual(values=PBM_colors) +
      scale_x_continuous(labels=scaleFUN) + scale_y_continuous(labels=scaleFUN) +
      ggtitle(paste(PBM_cond[i], " vs. ", PBM_cond[j], sep="")) +
      ylab(paste(PBM_cond[i], " -log10(q)", sep="")) + xlab(paste(PBM_cond[j], " -log10(q)", sep="")) +
      theme(legend.position="right") +
      #geom_text(aes(label=ifelse(yval>1.30103 & xval>1.30103, as.character(SNP_id), "")), hjust=1, vjust=1, size=3) +
      geom_hline(yintercept=1.30103, linetype="dashed", color = "grey") +
      geom_vline(xintercept=1.30103, linetype="dashed", color = "grey") +
      theme_cowplot(font_size = 14)
    
    # assign as a separate variable
    assign(paste("ggp_", j, sep=""), ggp)
    
    # use cowplot to save the scatterplots
    pdf_name <- paste("pairwise", PBM_cond[i], PBM_cond[j], "scatter.pdf", sep="_")
    ggp_grid <- plot_grid(plotlist=mget(paste("ggp_", as.character(j), sep="")), ncol = 1, nrow = 1, label_size = 8)
    save_plot(paste(getwd(), "Fig3_SNP_screen_COF", "pairwise_qval_scatters_NO_LABELS", pdf_name, sep="/"), ggp_grid, ncol = 1, nrow = 1, base_height=4, base_aspect_ratio=1.4, limitsize=FALSE)
    
    # remove g_df to prevent overwriting
    rm(g_df, ggp)
    
  }
}



### F4 - SNP CASCADE --------------------------------------------------------------------------

### PLOT MOTIF META-MATRIX GRID 

# import libraries
library(ggplot2)
library(cowplot)
library(ggseqlogo)
library(stringr)
library(RColorBrewer)

# function to obtain SNV value matrix (MODIFIED from the original CASCADE tile logo plotting script module)
zscore_to_SNV_matrix <- function(PBM, PBM_col, target_seq_col, target_seq) {
  
  # initialize a matrix to hold the SNV probe z-score values
  nucs <- c("A", "C", "G", "T")
  SNV <- matrix(nrow=4, ncol=nchar(as.character(target_seq)))
  
  # collect z-score values for each position in the matrix
  for (i in 1:nrow(SNV)) {
    for (j in 1:ncol(SNV)) {
      # construct the current sequence of interest using the current SNV
      SNV_seq <- paste(substr(target_seq, 1, j-1), nucs[i], substr(target_seq, j+1, ncol(SNV)), sep="")
      stopifnot(SNV_seq %in% as.character(PBM[, target_seq_col]))
      
      # collect the z-score for that given SNV probe sequence and condition
      SNV[i, j] <- PBM[which(as.character(PBM[, target_seq_col])==SNV_seq), PBM_col]
    }
  }
  
  # return the SNV energy matrix to the parent function
  rownames(SNV) <- nucs
  return(SNV)
}

# function to transform a SNV energy matrix using the column-wise median to plot as energy logo
SNV_to_SNV_trans <- function(SNV) {
  
  # transform z-score values to reflect deviation from column-wise median
  SNV_trans <- matrix(nrow=4, ncol=ncol(SNV))
  rownames(SNV_trans) <- c("A", "C", "G", "T")
  for (i in 1:nrow(SNV_trans)) {
    for (j in 1:ncol(SNV_trans)) {
      SNV_trans[i, j] <- SNV[i, j] - median(SNV[,j])
    }
  }
  
  # return the PWM matrix to be plotted as a sequence logo
  return(SNV_trans)
}

# function to transfrom a SNV energy matrix to PWM
SNV_to_PWM <- function(SNV, beta) {
  
  # transform using beta parameter
  SNV_pwm <- exp(beta * SNV)
  
  # compute probabilities for PWM using the column-wise totals
  PWM <- matrix(nrow=4, ncol=ncol(SNV_pwm))
  for (i in 1:nrow(PWM)) {
    for (j in 1:ncol(PWM)) {
      PWM[i, j] <- SNV_pwm[i, j]/sum(SNV_pwm[,j])
    }
  }
  rownames(PWM) <- c("A", "C", "G", "T")
  
  # return the PWM to parent function
  return(PWM)
  
}

# function to obtain reverse complement PWM
rev_compl_pwm <- function(m) {
  # reverse the input pwm
  rev_pwm <- m[,ncol(m):1]
  
  # reverse complement the pwm (A to T, C to G)
  rc_pwm <- rev_pwm
  rc_pwm[1,] <- rev_pwm[4,]
  rc_pwm[2,] <- rev_pwm[3,]
  rc_pwm[3,] <- rev_pwm[2,]
  rc_pwm[4,] <- rev_pwm[1,]
  rownames(rc_pwm) <- rownames(m)
  
  # return reversed pwm
  return(rc_pwm)
}

# function to obtain reverse complement of a DNA sequence
rev_compl <- function(DNA_seq){
  # returns the input DNA sequence in reverse
  temp_seq <- unlist(strsplit(DNA_seq, split=''))
  temp_seq <- rev(temp_seq)
  temp_seq <- paste(temp_seq, collapse='')
  temp_seq <- chartr("ACGTN", "TGCAN", temp_seq)
  return(temp_seq)
}

# function to convert an overlap_code to a real experimental expectation
convert_olap_code <- function(experiment_code, experiment_vector) {
  
  # enforce that the experiment code and the experiment vector have the same number of entries prior to conversion
  message(seed_name_list[[i]])
  message(experiment_code)
  message(paste(experiment_vector, collapse=" "))
  stopifnot(nchar(experiment_code)==length(experiment_vector))
  
  # initialize a converted code to return later
  converted_code <- ""
  
  # build a converted code by assigning each "1" in the overlap code its corresponding TF/COF name
  for (idx in 1:nchar(experiment_code)) {
    if (substr(experiment_code, idx, idx)=="1") { 
      converted_code <- paste(converted_code, experiment_vector[idx], sep=" ")
    }
  }
  
  # return the converted code
  return(converted_code)
  
}

# scale transformation function for rounding axes
scaleFUN <- function(x) sprintf("%.1f", x)

# create a directory to house the results folders for this analysis
dir.create(paste(getwd(), "Fig4_SNP_CASCADE", sep="/"))

# create a directory to house the motif grid
dir.create(paste(getwd(), "Fig4_SNP_CASCADE", "motif_grids", sep="/"))

# create a directory to house the SNV motif matrices to have on-hand
dir.create(paste(getwd(), "Fig4_SNP_CASCADE", "CASCADE_energy_matrices", sep="/"))
  
# create a directory to house the SNV motif matrices to have on-hand
dir.create(paste(getwd(), "Fig4_SNP_CASCADE", "CASCADE_PWMs", sep="/"))

# assign a seed z-score significance threshold based on a user-specified argument
zscore_thresh <- 1.5

# read in the input z-score file
df <- read.table("suppl_data_3_SNP_CASCADE.dat", sep='\t', header=T, stringsAsFactors = F)

# initialize probe_group variable
df$probe_group <- rep("", nrow(df))
for (i in 1:nrow(df)) {
  # assign as background if current probe is a background probe
  if (df[i, "probe_type"]=="BACKGROUND") {
    df[i, "probe_group"] <- "BACKGROUND"
    
    # assign as a seed probe if current probe is a seed
  } else if (df[i, "probe_type"]=="TILE" & is.na(df[i, "seed_nuc"])) {
    df[i, "probe_group"] <- "SEED"
    
    # assign as a SNV probe if it is neither a background nor a seed probe
  } else {
    df[i, "probe_group"] <- "SNV"
  }
}

# read in the gold-standard set annotation and convert overlap codes to 0-padded strings
gold_df <- read.table("SNP_SCREEN_REPRODUCIBLE_SNPs.dat", header=T, sep='\t', stringsAsFactors = F)
gold_df$overlap_code_screen <- str_pad(as.character(gold_df$overlap_code_screen), 8, pad = "0")
gold_df$overlap_code_cascade <- str_pad(as.character(gold_df$overlap_code_cascade), 8, pad = "0")
gold_df$overlap_code_shared <- str_pad(as.character(gold_df$overlap_code_shared), 8, pad = "0")

# hard code the experiment order that the overlap codes correspond to
exp_vec <- c("p300", "SMARCA4", "TBL1XR1", "RBBP5", "PU.1", "p65", "HDAC1", "GCN5")

# clone the PBM_cond variable to be able to modify it within the function
PBM_cond <- read.table("SNP_CASCADE_exp.txt", sep="\t")
PBM_cond <- as.character(PBM_cond$V1)
PBM_cond_temp <- PBM_cond

# isolate the SNP_cat in alphabetical order
SNP_cat <- unique(as.character(df$SNP_cat))
SNP_cat <- SNP_cat[order(SNP_cat[!is.na(SNP_cat)])]

# initialize PBM colors
PBM_colors <- brewer.pal(length(SNP_cat)-1, "Set2")
PBM_colors <- c("black", PBM_colors)
PBM_colors <- PBM_colors[c(1, 5, 6, 2, 3, 4, 7)]
PBM_colors <- PBM_colors[2:length(PBM_colors)]

# generate a different FULL motif file for each SNP_cat
for (s in SNP_cat) {
  
  # grab the seeds from the current CASCADE design for each SNP category
  y <- droplevels(df[which(df$SNP_cat==s),])
  seeds <- unique(as.character(y$seed_names))
  rm(y)
  
  # determine which PBM_color color to use to represent the given SNP category s
  switch_color <- switch(s,
                         "basal_eQTL" = PBM_colors[3],
                         "bg_SNP" = PBM_colors[4],
                         "caQTL_eQTL" = PBM_colors[5],
                         "fairfax_eQTL" = PBM_colors[6],
                         "GWAS_caQTL" = PBM_colors[1],
                         "GWAS_eQTL" = PBM_colors[2],
                         "resp_eQTL" = PBM_colors[6])

  # restrict to "br" orientation
  or <- "br"
  for (o in or) {
    
    # initialize a list for the plots
    snv_list <- list()
    seed_name_list <- list()
    seed_seq_list <- list()
    seed_score_list <- list()
    PBM_cond_list <- list()
    SNP_allele_list <- list()
    olap_screen_list <- list()
    olap_cascade_list <- list()
    reprod_list <- list()
    
    # modify the PBM_cond vector 
    PBM_cond <- paste("v1_a1_run1", o, PBM_cond_temp, sep="_")
    
    # for each list file on the command line
    for(q in 1:length(PBM_cond)){
      
      # determine beta parameter value for the given experiment
      beta <- 30/max(df[, PBM_cond[q]])
      
      # create a short condition name using the orientation
      short_cond <- gsub("v1_a1_run1", "", PBM_cond[q])
      short_cond <- gsub(paste("_", o, "_", sep=""), "", short_cond)
      
      # for each seed in the given CASCADE design
      for (r in 1:length(seeds)) {
        
        # isolate probes relevant to the current seed
        x <- df[which(df$seed_names==seeds[r]),]
        stopifnot(nrow(x)==79)
        
        # determine the reference seq and score for the given seed
        ref_seed <- as.character(x[which(x$probe_group=="SEED"), "target_seq"])
        ref_score <- x[which(x$probe_group=="SEED"), PBM_cond[q]]
        
        # determine the current SNP_allele using the reference seed
        curr_SNP_allele <- as.character(x[which(x$probe_group=="SEED"), "SNP_allele"])
        curr_SNP_id <- as.character(x[which(x$probe_group=="SEED"), "SNP_id"])
        
        # use current SNP_id to determine the overlap and reproducibility categories for the given SNP
        curr_olap_screen <- as.character(gold_df[which(gold_df$SNP_id==curr_SNP_id), "overlap_code_screen"])
        curr_olap_cascade <- as.character(gold_df[which(gold_df$SNP_id==curr_SNP_id), "overlap_code_cascade"])
        curr_reprod <- as.character(gold_df[which(gold_df$SNP_id==curr_SNP_id), "reprod_tier"])
        
        # replace any NAs if necessary
        if ( length(curr_olap_screen)==0 ) { curr_olap_screen <- "XXXXXXXX" }
        if ( length(curr_olap_cascade)==0 ) { curr_olap_cascade <- "XXXXXXXX" }
        if ( length(curr_reprod)==0 ) { curr_reprod <- "FAIRFAX" }
        
        # replace with a pseudocount on the off-chance that an NA was the reference score
        if (is.na(ref_score)) {
          ref_score <- 0.0000001
        }
        
        # build the SNV matrix using the current data
        snv <- zscore_to_SNV_matrix(x, PBM_cond[q], "target_seq", ref_seed)
        
        # in the event that NAs still exist, replace with pseudo
        snv[is.na(snv)] <- 0.0000001
        
        # output the untransformed SNV energy matrix
        mat_name <- paste(seeds[r], short_cond, o, "RAWSNV.txt", sep="_")
        write.table(snv, paste(getwd(), "Fig4_SNP_CASCADE", "CASCADE_energy_matrices", mat_name, sep="/"), quote=F, row.names=T, col.names=F, sep='\t')
        
        # transform the SNV matrix using the column-wise medians
        snv_trans <- SNV_to_SNV_trans(snv)
        
        # output the transformed SNV energy matrix
        mat_name <- paste(seeds[r], short_cond, o, "SNV.txt", sep="_")
        write.table(snv_trans, paste(getwd(), "Fig4_SNP_CASCADE", "CASCADE_energy_matrices", mat_name, sep="/"), quote=F, row.names=T, col.names=F, sep='\t')
        
        # transform the SNV matrix to a PWM for analysis (if it comes from a seed greater than the threshold)
        if (ref_score >= zscore_thresh & sum(is.na(snv))==0 & sum(is.na(SNV_to_PWM(snv, beta)))==0) {
          pwm <- SNV_to_PWM(snv, beta)
          
          # output the transformed SNV energy matrix
          mat_name <- paste(seeds[r], short_cond, o, "PWM.txt", sep="_")
          write.table(pwm, paste(getwd(), "Fig4_SNP_CASCADE", "CASCADE_PWMs", mat_name, sep="/"), quote=F, row.names=T, col.names=F, sep='\t')
          
          # clean up
          rm(pwm)
        }
        
        # push the current PWM to the pwm list
        snv_list[[r+length(seeds)*(q-1)]] <- snv_trans
        
        # push the current header label to the label list
        seed_name_list[[r+length(seeds)*(q-1)]] <- seeds[r]
        seed_seq_list[[r+length(seeds)*(q-1)]] <- ref_seed
        seed_score_list[[r+length(seeds)*(q-1)]] <- ref_score
        PBM_cond_list[[r+length(seeds)*(q-1)]] <- PBM_cond[q]
        SNP_allele_list[[r+length(seeds)*(q-1)]] <- curr_SNP_allele
        olap_screen_list[[r+length(seeds)*(q-1)]] <- curr_olap_screen
        olap_cascade_list[[r+length(seeds)*(q-1)]] <- curr_olap_cascade
        reprod_list[[r+length(seeds)*(q-1)]] <- curr_reprod
        
        # clean up
        rm(x, ref_seed, ref_score, curr_SNP_allele, curr_SNP_id, curr_olap_screen, curr_olap_cascade, curr_reprod, snv, snv_trans)
        
      }
      
      # clean up
      rm(beta, short_cond)
      
    }
    
    # define dimensions for pdf page
    n_col <- length(PBM_cond)
    n_row <- length(seeds)
    
    # for each possible sequence direction
    for (dir in c("FORWARD", "REVERSE")) {
      
      # save each logo into an individually numbered variable
      for (i in 1:length(snv_list)) {
        
        # save the seed sequence and motif as variables so they can be modified if necessary
        curr_snv <- snv_list[[i]]
        curr_seed_seq <- seed_seq_list[[i]]
        
        # reverse complement the seed sequence and motif if the current direction is REVERSE
        if (dir == "REVERSE") {
          curr_snv <- rev_compl_pwm(curr_snv)
          curr_seed_seq <- rev_compl(curr_seed_seq)
        }
        
        # convert overlap codes prior to plotting the title
        screen_code <- convert_olap_code(olap_screen_list[[i]], exp_vec)
        cascade_code <- convert_olap_code(olap_cascade_list[[i]], exp_vec)
        if (nchar(screen_code)==0) { screen_code <- " NONE" }
        if (nchar(cascade_code)==0) { cascade_code <- " NONE"}
        
        # plot the SNV energy logo
        ggp <- ggseqlogo(curr_snv, method="custom", seq_type="dna")
        my.ggp.yrange <- ggplot_build(ggp)$layout$panel_params[[1]]$y.range
        ggp <- ggp +
          ggtitle(paste(PBM_cond_list[[i]], paste(seed_name_list[[i]], SNP_allele_list[[i]], sep='\t'), paste("SCREEN:", screen_code, sep=""), paste("CASCADE:", cascade_code, sep=""), curr_seed_seq, seed_score_list[[i]], sep='\n')) +
          scale_x_continuous(expand = c(0.01,0.01)) +
          scale_y_continuous(expand = c(0.01,0.01), labels=scaleFUN) +
          annotate('rect', xmin = 0.5, xmax = ncol(curr_snv)+0.5, ymin = my.ggp.yrange[1], ymax = 0.0, alpha = 0.75, col='white', fill='white')
        
        # place a transparent window over the logo if the seed score is below threshold
        ggp <- ggp +
          if (seed_score_list[[i]] < zscore_thresh) {
            annotate('rect', xmin = 0.5, xmax = ncol(curr_snv)+0.5, ymin = my.ggp.yrange[1], ymax= my.ggp.yrange[2], alpha = 0.75, col='white', fill='white')
          }
        
        # plot the borders and axis label
        ggp <- ggp +
          annotate('rect', xmin = 0.5, xmax = ncol(curr_snv)+0.5, ymin = 0, ymax = my.ggp.yrange[2], col="lightgray", fill=NA, size=0.5) +
          annotate('rect', xmin = 0.5, xmax = ncol(curr_snv)+0.5, ymin = my.ggp.yrange[1], ymax = 0, col="lightgray", fill=NA, size=0.5) +
          ylab(expression(paste(Delta, "z-score", sep=""))) +
          theme(axis.text.x=element_blank()) +
          theme(axis.title.x=element_blank()) +
          theme(axis.text.y=element_text(size=4, face="plain"), axis.title.y=element_text(size=4, face="plain")) +
          theme(plot.title = element_text(size=4))
        
        # place a colored semi-transparent window over the SNP position
        ggp <- ggp +
          if (dir == "FORWARD") {
            annotate('rect', xmin = 13.5, xmax = 14.5, ymin = my.ggp.yrange[1], ymax= my.ggp.yrange[2], alpha = 0.25, fill=switch_color)
          } else {
            annotate('rect', xmin = 12.5, xmax = 13.5, ymin = my.ggp.yrange[1], ymax= my.ggp.yrange[2], alpha = 0.25, fill=switch_color)
          }
        
        # use a switch to determine what color to make the title based on the SNP's given reproducibility category
        title_color <- switch(reprod_list[[i]],
                              "GOLD" = "darkgoldenrod3",
                              "FAIRFAX" = "darkgoldenrod3",
                              "SILVER" = "darkgray",
                              "NONE" = "black")
        
        # color the title based on the switch value above
        ggp <- ggp +
          theme(plot.title = element_text(color=title_color))
        
        # assign as a separate variable
        assign(paste("ggp_", i, sep=""), ggp)
        
      }
      
      # generate a numerical pattern for cowplot by collapsing a matrix
      cowplot_indeces <- matrix(1:(n_row*n_col), n_row)
      cowplot_indeces <- as.vector(t(cowplot_indeces))
      
      # use cowplot to save the scatterplots
      ggp_grid <- plot_grid(plotlist=mget(paste0("ggp_", cowplot_indeces)), ncol = n_col, nrow = n_row, label_size = 2)
      
      # save PDF
      pdf_name <- paste("SNP_CASCADE", s, dir, o, "motif_grid.pdf", sep="_")
      save_plot(paste(getwd(), "Fig4_SNP_CASCADE", "motif_grids", pdf_name, sep="/"), ggp_grid, ncol = n_col, nrow = n_row, base_height=1.2, base_aspect_ratio=1.5, limitsize=FALSE)

      # remove temp ggp variables used for cowplot
      rm(list=paste0("ggp_", cowplot_indeces))
      
    }
    
    # remove list variables in between SNP categories
    rm(snv_list,
       seed_name_list,
       seed_seq_list,
       seed_score_list,
       PBM_cond_list,
       SNP_allele_list,
       olap_screen_list,
       olap_cascade_list,
       reprod_list)
    
  }
}

# read in the paper SNP lookup table
paper_SNPs <- read.table("SNP_CASCADE_Fig4.txt", header=T, sep='\t', stringsAsFactors = F)

# generate a motif file for just the motifs highlighted in the paper
or <- "br"
for (o in or) {
  
  # initialize vector of seeds used in the paper
  seeds <- paper_SNPs$seed_names
  
  # initialize a list for the plots
  snv_list <- list()
  seed_name_list <- list()
  seed_seq_list <- list()
  seed_score_list <- list()
  PBM_cond_list <- list()
  SNP_allele_list <- list()
  olap_screen_list <- list()
  olap_cascade_list <- list()
  reprod_list <- list()
  
  # modify the PBM_cond vector 
  PBM_cond <- paste("v1_a1_run1", o, PBM_cond_temp, sep="_")
  PBM_cond <- PBM_cond[1:(length(PBM_cond)-2)]
  
  # for each list file on the command line
  for(q in 1:length(PBM_cond)){
    
    # determine beta parameter value for the given experiment
    beta <- 30/max(df[, PBM_cond[q]])
    
    # create a short condition name using the orientation
    short_cond <- gsub("v1_a1_run1", "", PBM_cond[q])
    short_cond <- gsub(paste("_", o, "_", sep=""), "", short_cond)
    
    # for each seed in the given CASCADE design
    for (r in 1:length(seeds)) {
      
      # isolate probes relevant to the current seed
      x <- df[which(df$seed_names==seeds[r]),]
      stopifnot(nrow(x)==79)
      
      # determine the reference seq and score for the given seed
      ref_seed <- as.character(x[which(x$probe_group=="SEED"), "target_seq"])
      ref_score <- x[which(x$probe_group=="SEED"), PBM_cond[q]]
      
      # determine the current SNP_allele using the reference seed
      curr_SNP_allele <- as.character(x[which(x$probe_group=="SEED"), "SNP_allele"])
      curr_SNP_id <- as.character(x[which(x$probe_group=="SEED"), "SNP_id"])
      
      # use current SNP_id to determine the overlap and reproducibility categories for the given SNP
      curr_olap_screen <- as.character(gold_df[which(gold_df$SNP_id==curr_SNP_id), "overlap_code_screen"])
      curr_olap_cascade <- as.character(gold_df[which(gold_df$SNP_id==curr_SNP_id), "overlap_code_cascade"])
      curr_reprod <- as.character(gold_df[which(gold_df$SNP_id==curr_SNP_id), "reprod_tier"])
      
      # replace any NAs if necessary
      if ( length(curr_olap_screen)==0 ) { curr_olap_screen <- "XXXXXXXX" }
      if ( length(curr_olap_cascade)==0 ) { curr_olap_cascade <- "XXXXXXXX" }
      if ( length(curr_reprod)==0 ) { curr_reprod <- "FAIRFAX" }
      
      # replace with a pseudocount on the off-chance that an NA was the reference score
      if (is.na(ref_score)) {
        ref_score <- 0.0000001
      }
      
      # build the SNV matrix using the current data
      snv <- zscore_to_SNV_matrix(x, PBM_cond[q], "target_seq", ref_seed)
      
      # in the event that NAs still exist, replace with pseudo
      snv[is.na(snv)] <- 0.0000001
      
      # transform the SNV matrix using the column-wise medians
      snv_trans <- SNV_to_SNV_trans(snv)
      
      # push the current PWM to the pwm list
      snv_list[[r+length(seeds)*(q-1)]] <- snv_trans
      
      # push the current header label to the label list
      seed_name_list[[r+length(seeds)*(q-1)]] <- seeds[r]
      seed_seq_list[[r+length(seeds)*(q-1)]] <- ref_seed
      seed_score_list[[r+length(seeds)*(q-1)]] <- ref_score
      PBM_cond_list[[r+length(seeds)*(q-1)]] <- PBM_cond[q]
      SNP_allele_list[[r+length(seeds)*(q-1)]] <- curr_SNP_allele
      olap_screen_list[[r+length(seeds)*(q-1)]] <- curr_olap_screen
      olap_cascade_list[[r+length(seeds)*(q-1)]] <- curr_olap_cascade
      reprod_list[[r+length(seeds)*(q-1)]] <- curr_reprod
      
      # clean up
      rm(x, ref_seed, ref_score, curr_SNP_allele, curr_SNP_id, curr_olap_screen, curr_olap_cascade, curr_reprod, snv, snv_trans)
      
    }
    
    # clean up
    rm(beta, short_cond)
    
  }
  
  # define dimensions for pdf page
  n_col <- length(PBM_cond)
  n_row <- length(seeds)
  
  # save each logo into an individually numbered variable
  for (i in 1:length(snv_list)) {
    
    # save the seed sequence and motif as variables so they can be modified if necessary
    curr_snv <- snv_list[[i]]
    curr_seed_seq <- seed_seq_list[[i]]
    
    # determine which PBM_color color to use to represent the given SNP category s
    s <- paper_SNPs[which(paper_SNPs$seed_names==seed_name_list[[i]]), "SNP_cat"]
    switch_color <- switch(s,
                           "basal_eQTL" = PBM_colors[3],
                           "bg_SNP" = PBM_colors[4],
                           "caQTL_eQTL" = PBM_colors[5],
                           "fairfax_eQTL" = PBM_colors[6],
                           "GWAS_caQTL" = PBM_colors[1],
                           "GWAS_eQTL" = PBM_colors[2],
                           "resp_eQTL" = PBM_colors[6])
    
    # fetch the direction of the current seed
    dir <- paper_SNPs[which(paper_SNPs$seed_names==seed_name_list[[i]]), "plot_dir"]
    
    # reverse complement the seed sequence and motif if the current direction is REVERSE
    if (dir == "REVERSE") {
      curr_snv <- rev_compl_pwm(curr_snv)
      curr_seed_seq <- rev_compl(curr_seed_seq)
    }
    
    # convert overlap codes prior to plotting the title
    screen_code <- convert_olap_code(olap_screen_list[[i]], exp_vec)
    cascade_code <- convert_olap_code(olap_cascade_list[[i]], exp_vec)
    if (nchar(screen_code)==0) { screen_code <- " NONE" }
    if (nchar(cascade_code)==0) { cascade_code <- " NONE"}
    
    # plot the SNV energy logo
    ggp <- ggseqlogo(curr_snv, method="custom", seq_type="dna")
    my.ggp.yrange <- ggplot_build(ggp)$layout$panel_params[[1]]$y.range
    ggp <- ggp +
      ggtitle(paste(PBM_cond_list[[i]], paste(seed_name_list[[i]], SNP_allele_list[[i]], sep='\t'), paste("SCREEN:", screen_code, sep=""), paste("CASCADE:", cascade_code, sep=""), curr_seed_seq, seed_score_list[[i]], sep='\n')) +
      scale_x_continuous(expand = c(0.01,0.01)) +
      scale_y_continuous(expand = c(0.01,0.01), labels=scaleFUN) +
      annotate('rect', xmin = 0.5, xmax = ncol(curr_snv)+0.5, ymin = my.ggp.yrange[1], ymax = 0.0, alpha = 0.75, col='white', fill='white')
    
    # place a transparent window over the logo if the seed score is below threshold
    ggp <- ggp +
      if (seed_score_list[[i]] < zscore_thresh) {
        annotate('rect', xmin = 0.5, xmax = ncol(curr_snv)+0.5, ymin = my.ggp.yrange[1], ymax= my.ggp.yrange[2], alpha = 0.75, col='white', fill='white')
      }
    
    # plot the borders and axis label
    ggp <- ggp +
      annotate('rect', xmin = 0.5, xmax = ncol(curr_snv)+0.5, ymin = 0, ymax = my.ggp.yrange[2], col="lightgray", fill=NA, size=0.5) +
      annotate('rect', xmin = 0.5, xmax = ncol(curr_snv)+0.5, ymin = my.ggp.yrange[1], ymax = 0, col="lightgray", fill=NA, size=0.5) +
      ylab(expression(paste(Delta, "z-score", sep=""))) +
      theme(axis.text.x=element_blank()) +
      theme(axis.title.x=element_blank()) +
      theme(axis.text.y=element_text(size=4, face="plain"), axis.title.y=element_text(size=4, face="plain")) +
      theme(plot.title = element_text(size=4))
    
    # place a colored semi-transparent window over the SNP position
    ggp <- ggp +
      if (dir == "FORWARD") {
        annotate('rect', xmin = 13.5, xmax = 14.5, ymin = my.ggp.yrange[1], ymax= my.ggp.yrange[2], alpha = 0.25, fill=switch_color)
      } else {
        annotate('rect', xmin = 12.5, xmax = 13.5, ymin = my.ggp.yrange[1], ymax= my.ggp.yrange[2], alpha = 0.25, fill=switch_color)
      }
    
    # use a switch to determine what color to make the title based on the SNP's given reproducibility category
    title_color <- switch(reprod_list[[i]],
                          "GOLD" = "darkgoldenrod3",
                          "FAIRFAX" = "darkgoldenrod3",
                          "SILVER" = "darkgray",
                          "NONE" = "black")
    
    # color the title based on the switch value above
    ggp <- ggp +
      theme(plot.title = element_text(color=title_color))
    
    # assign as a separate variable
    assign(paste("ggp_", i, sep=""), ggp)
    
  }
  
  # generate a numerical pattern for cowplot by collapsing a matrix
  cowplot_indeces <- matrix(1:(n_row*n_col), n_row)
  cowplot_indeces <- as.vector(t(cowplot_indeces))
  
  # use cowplot to save the scatterplots
  ggp_grid <- plot_grid(plotlist=mget(paste0("ggp_", cowplot_indeces)), ncol = n_col, nrow = n_row, label_size = 2)
  
  # save PDF
  pdf_name <- paste("SNP_CASCADE", "figure", "motif_grid.pdf", sep="_")
  save_plot(paste(getwd(), "Fig4_SNP_CASCADE", "motif_grids", pdf_name, sep="/"), ggp_grid, ncol = n_col, nrow = n_row, base_height=1.2, base_aspect_ratio=1.5, limitsize=FALSE)
  
  # save SVG
  svg_name <- paste("SNP_CASCADE", "figure", "motif_grid.svg", sep="_")
  save_plot(paste(getwd(), "Fig4_SNP_CASCADE", "motif_grids", svg_name, sep="/"), ggp_grid, ncol = n_col, nrow = n_row, base_height=1.2, base_aspect_ratio=1.5, limitsize=FALSE)
  
  # remove temp ggp variables used for cowplot
  rm(list=paste0("ggp_", cowplot_indeces))
  
  # remove list variables in between SNP categories
  rm(snv_list,
     seed_name_list,
     seed_seq_list,
     seed_score_list,
     PBM_cond_list,
     SNP_allele_list,
     olap_screen_list,
     olap_cascade_list,
     reprod_list)
  
}



### prep an input file for TOMTOM analysis

# collect vector of PWM.txt files
setwd(paste(getwd(), "Fig4_SNP_CASCADE", "CASCADE_PWMs", sep="/"))
motif_files <- list.files(path=getwd(), pattern="_PWM.txt")

# if the vector isn't empty
if (length(motif_files)>0) {
  
  # intialize a MEME-formatted motif file to write the reformatted motifs to
  meme_file_name <- paste("SNP_CASCADE", "_motifs.meme", sep="")
  file.create(meme_file_name)
  meme_file <- file(meme_file_name)
  writeLines(c("MEME version 5.0.3\n",
               "\n",
               "ALPHABET= ACGT\n",
               "\n",
               "strands: + -\n",
               "Background letter frequencies\n",
               "A 0.250 C 0.250 G 0.250 T 0.250 \n",
               "\n"),
             con=meme_file, sep="")
  
  
  # for each motif in the motif_files vector
  for (m_file in motif_files) {
    
    # read the given file
    m <- as.matrix(read.table(m_file, row.names=1))
    
    # transpose the matrix
    m <- t(m)
    
    # create a header for the given motif
    m_header <- gsub("_PWM.txt", "", m_file)
    
    # create a subheader for the given motif
    m_subheader <- paste("letter-probability",
                         "matrix:",
                         "alength= 4",
                         "w=", as.character(nrow(m)),
                         "nsites= 0",
                         "E= 0")
    
    # write the header and matrix to file
    write(paste("MOTIF", m_header, sep=" "), meme_file_name, append=T)
    write(m_subheader, meme_file_name, append=T)
    write(formatC(t(m), format="f", 6), meme_file_name, ncolumns=4, sep=" ", append=T)
    write("", meme_file_name, append=T)
  }
  close(meme_file)
  
}
