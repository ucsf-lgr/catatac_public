#################################
## CROP-SEQ ANALYSIS FUNCTIONS ##
#################################

# Compilation of functions to ingest, process, and analyse CROP-seq data
# Data is stored as `SingleCellExperiment` objects 
library(data.table)
library(dplyr)
library(Matrix)
library(ggplot2)
library(iSEE)

# Miscellaneous -----------------------------------------------------------
th <- theme_bw() + theme(
  axis.text.x = element_text(size=12), 
  axis.title.x = element_text(size=12), 
  axis.text.y = element_text(size=12), 
  axis.title.y = element_text(size=12),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  axis.line = element_line(colour = "black"), 
  panel.border = element_blank(), 
  plot.title = element_text(face="bold", hjust = 0.5, size=8),
  plot.subtitle = element_text(hjust = 0.5))



# gRNA calling ------------------------------------------------------------
## function using binomial or negative binomial test
call_gRNAs_binom <- function(sce=NULL, 
                                      pval_threshold=0.001,
                                      max_to_test=5, 
                                      total_guides=NULL,
                                      prob=NULL,
                                      min_count=3, 
                                      min_fraction=0.05, 
                                      plot=TRUE, 
                                      res=200,
                                      ann=NULL, 
                                      columns=NULL, 
                                      add_stats=TRUE, 
                                      test="binomial",
                                      overdisp=NULL){
  if(is.null(total_guides)){
    stop("Please provide the total number of gRNAs to use in the testing. 
         This should be the total number of gRNAs present in the library.")
  }
  
  if(!(test %in% c("binomial", "negative binomial"))){
    stop("This function currently supports binomial and negative binomial tests only.  
         Please specify binomial or negative binomial test to use.")
  }
  
  cat(paste0("\n The ", test, " test will be used.\n\n"))
  
  ## create data table with detected gRNA-cell combinations
  gRNA_calls <- reshape2::melt(as.matrix(counts(altExp(sce, 'CRISPR'))))
  ## remove guide-cell pairs with no counts
  gRNA_calls <- gRNA_calls[gRNA_calls$value>0,]
  colnames(gRNA_calls) <- c("gRNA", "cell", "UMI_count")
  
  ## remove guides with too few UMIs
  gRNA_calls <- gRNA_calls %>% filter(UMI_count > min_count)
  
  ## annotate the targets associated with each gRNA
  if(!is.null(ann)){
    if(length(setdiff(columns, colnames(ann)))>0){
      stop(paste("Column(s)", setdiff(columns, colnames(ann)), "are not present in annotation data.frame"))
    }
    for(i in seq_along(columns)){
      gRNA_calls[,columns[i]] <- ann[match(gRNA_calls$gRNA, ann$ID), columns[i]]
    }
  }
  
  ## add total UMIs per cell, to restrict analysis to the largest counts per cell
  totals <- colSums(counts(altExp(sce, 'CRISPR')))
  gRNA_calls$cell_total_UMIs <- totals[gRNA_calls$cell]
  # compute % of total UMIs taken by each gRNA
  gRNA_calls$fraction_total_UMIs <- gRNA_calls$UMI_count/gRNA_calls$cell_total_UMIs
  
  ## convert to data.table using tibble
  gRNA_calls <- as.data.frame(gRNA_calls)
  
  ## order decreasingly by UMI_count
  gRNA_calls <- gRNA_calls %>% group_by(cell) %>% arrange(cell, -UMI_count)
  
  # add column with order per cell based on UMI abundance
  gRNA_calls <- gRNA_calls  %>% mutate(order = 1:n()) %>% ungroup()
  
  # test top gRNAs per cell as determined by max_to_test
  # if more than max_to_test integrations are true, the rest won't be called
  # higher number in max_to_test results in higher computation time
  
  ## if no probabilitis provided, use equal probabilities
  if(is.null(prob)){
    gRNA_calls$prob <- 1/total_guides
  }else{
    ## add proportions of guides in guide library to gRNA_calls
    ## provided as a named vector, with names being gRNA ids
    prob <- as.data.table(prob, keep.rownames = TRUE)
    gRNA_calls <- merge(gRNA_calls, prob, by.x="gRNA", by.y="rn")
  }
  gRNA_calls <- gRNA_calls %>% 
    mutate(gRNA_num = gRNA %>% as.factor() %>% as.numeric()) %>% 
    group_by(cell) %>% 
    mutate(means = mean(UMI_count)) %>% 
    ungroup
  
  ## use a binomial test to determine if observed counts are greater than background
  if(test == "binomial") {
    
    gRNA_calls <- gRNA_calls %>% 
      rowwise() %>% 
      mutate(binom = ifelse(order <= max_to_test, binom.test(x = UMI_count, n = cell_total_UMIs, p = prob, alternative = "greater")$p.value, NA))
    
    # bonferroni correction, per cell
    gRNA_calls <- gRNA_calls %>% 
      group_by(cell) %>% 
      mutate(n = n(),
             # reset pvals that are larger than 1 after correction
             p_adj_binom = ifelse(binom * n < 1, binom*n , 1),
             # count significant results
             n_detected_gRNAs_binom = length(which(p_adj_binom <= pval_threshold))
      ) %>% 
      ungroup()
    
    #### looking at binomial test results, because binom prop grows faster, nothing is significant using nbinom
    gRNA_calls <- gRNA_calls %>% mutate(n_detected_gRNAs = n_detected_gRNAs_binom)
    ## check if there are cells with at least max_to_test called hits
    # if so, retest a larger number
    while(max(gRNA_calls$n_detected_gRNAs) >= max_to_test){
      # get cells that need higher threshold
      done <- gRNA_calls %>% filter(n_detected_gRNAs < max_to_test)
      missing <- gRNA_calls %>% filter(n_detected_gRNAs >= max_to_test) %>%
        select(c(gRNA, cell, UMI_count, target, class, locus, cell_total_UMIs, fraction_total_UMIs,
                 order, prob, gRNA_num,means))
      
      # compute pvals for those cells with twice as many cells tested
      max_to_test <- max_to_test*2
      missing <- missing %>%
        rowwise() %>%
        mutate(binom = ifelse(order <= max_to_test, binom.test(x = UMI_count, n = cell_total_UMIs, p = prob, alternative = "greater")$p.value, NA))
      
      # bonferroni correction, per cell
      missing <- missing %>%
        group_by(cell) %>%
        mutate(n = n(),
               # reset pvals that are larger than 1 after correction
               p_adj_binom = ifelse(binom * n < 1, binom*n , 1),
               # count significant results
               n_detected_gRNAs_binom = length(which(p_adj_binom <= pval_threshold)),
               n_detected_gRNAs = n_detected_gRNAs_binom) %>%
        ungroup()
      
      # join with previous calls
      gRNA_calls <- rbind(done, missing)
    }
    
    #### working with binomial part, assign p_adj as p_adj_binom
    gRNA_calls <- gRNA_calls %>% mutate(p_adj = p_adj_binom) %>%
      group_by(cell) %>%
      mutate(n_detected_gRNAs = length(which(p_adj <= pval_threshold))) %>%
      ungroup()
  }
  
  ## use a binomial test to determine if observed counts are greater than background
  if(test == "negative binomial") {
    
    if(is.null(overdisp)) { #new
      cat("\n Estimating dispersion based on data \n")
      disp = gRNA_calls %>% summarise(glob_mean = round(mean(UMI_count),4), 
                                      glob_var = round(var(UMI_count),4), 
                                      glob_r = round(glob_mean^2/(glob_var-glob_mean),4), 
                                      glob_p = round(glob_r/(glob_mean+glob_r),4))
      cat("\n Estimates across cell and gRNAs of \n global mean: ", disp$glob_mean, "\n global variance: ", disp$glob_var, 
          "\n global dispersion: ", disp$glob_r, "\n probability of success: ", disp$glob_p, "\n")
      
      overdisp = disp$glob_r
    }
    
    gRNA_calls <- gRNA_calls %>% 
      rowwise() %>% 
      mutate(nbinom = ifelse(order <= max_to_test, pnbinom(UMI_count-1, mu = means, size = overdisp, lower.tail = F), NA)) #new
    
    # bonferroni correction, per cell
    gRNA_calls <- gRNA_calls %>% 
      group_by(cell) %>% 
      mutate(n = n(),
             # reset pvals that are larger than 1 after correction
             p_adj_nbinom = ifelse(nbinom * n < 1, nbinom*n, 1),
             # count significant results
             n_detected_gRNAs_nbinom = length(which(p_adj_nbinom <= pval_threshold))) %>% 
      ungroup()
    
    #### looking at negative binomial test results
    gRNA_calls <- gRNA_calls %>% mutate(n_detected_gRNAs = n_detected_gRNAs_nbinom) #redundant
    gRNA_calls <- gRNA_calls %>% mutate(p_adj = p_adj_nbinom) %>% #new
      group_by(cell) %>%
      mutate(n_detected_gRNAs = length(which(p_adj <= pval_threshold))) %>%
      ungroup()
    ## check if there are cells with at least max_to_test called hits
    # if so, retest a larger number
    while(max(gRNA_calls$n_detected_gRNAs) >= max_to_test){
      # get cells that need higher threshold
      done <- gRNA_calls %>% filter(n_detected_gRNAs < max_to_test)
      missing <- gRNA_calls %>% filter(n_detected_gRNAs >= max_to_test) %>%
        select(c(gRNA, cell, UMI_count, cell_total_UMIs, fraction_total_UMIs,
                 order, prob, gRNA_num,means))
      
      # compute pvals for those cells with twice as many cells tested
      max_to_test <- max_to_test*2
      missing <- missing %>%
        rowwise() %>%
        mutate(binom = ifelse(order <= max_to_test, pnbinom(UMI_count-1, mu = means, size = overdisp, lower.tail = F), NA)) #new
      
      
      # bonferroni correction, per cell
      missing <- missing %>%
        group_by(cell) %>%
        mutate(n = n(),
               # reset pvals that are larger than 1 after correction
               p_adj_nbinom = ifelse(nbinom * n < 1, nbinom*n, 1),
               # count significant results
               n_detected_gRNAs_nbinom = length(which(p_adj_nbinom <= pval_threshold)),
               n_detected_gRNAs = n_detected_gRNAs_nbinom) %>%
        ungroup()
      
      # join with previous calls
      gRNA_calls <- rbind(done, missing)
    }
    
    #### working with negative binomial part, assign p_adj as p_adj_nbinom
    gRNA_calls <- gRNA_calls %>% mutate(p_adj = p_adj_nbinom) %>%
      group_by(cell) %>%
      mutate(n_detected_gRNAs = length(which(p_adj <= pval_threshold))) %>%
      ungroup()
  }
  
  
  # for unique calls, also remove those where the fraction from total UMIs is too low
  ifelse(is.null(dim(gRNA_calls$n_detected_gRNAs == 1 &
                       gRNA_calls$fraction_total_UMIs <= min_fraction)), 
         print("no unique calls with total too low UMIs"), 
         gRNA_calls[gRNA_calls$n_detected_gRNAs == 1 &
                      gRNA_calls$fraction_total_UMIs <= min_fraction,]$p_adj <- 1)
  
  # count significant results again after adjustments above
  gRNA_calls <- gRNA_calls %>% 
    group_by(cell) %>% 
    mutate(n_detected_gRNAs = length(which(p_adj <= pval_threshold))) %>% 
    ungroup()
  # assign call based on number of significant gRNAs
  gRNA_calls$call <- ifelse(gRNA_calls$n_detected_gRNAs == 0, "unassigned",
                            ifelse(gRNA_calls$n_detected_gRNAs == 1, "unique", "multiple"))
  gRNA_long <- gRNA_calls
  # annotate significant results
  gRNA_calls <- gRNA_calls %>% group_by(cell) %>% 
    mutate(detected_gRNAs = paste0(gRNA[which(p_adj <= pval_threshold)], collapse="|"),
           detected_UMIs  = paste0(UMI_count[which(p_adj <= pval_threshold)], collapse="|"),
           detected_fraction = paste0(fraction_total_UMIs[which(p_adj <= pval_threshold)], collapse="|"),
           detected_padj = paste0(p_adj[which(p_adj <= pval_threshold)], collapse="|")
    ) %>% 
    ungroup()
  
  ## convert to datatable
  gRNA_calls <- data.table(gRNA_calls)
  
  ## add annotation if required
  if(!is.null(ann)){
    for(i in seq_along(columns)){
      gRNA_calls[, tmp := .(paste0(get(columns[i])[which(p_adj <= pval_threshold)], collapse="|")), by = "cell"]
      colnames(gRNA_calls)[ncol(gRNA_calls)] <- paste0('detected_', columns[i])
    }
  }
  
  ## QC plot
  if(plot){
    ## remove overlapping points
    idx <- subsetPointsByGrid(X = gRNA_calls$fraction_total_UMIs,
                              Y = log10(gRNA_calls$UMI_count),
                              grouping = gRNA_calls$cell,
                              resolution = res)
    df <- gRNA_calls[idx,]
    df$call <- factor(df$call, levels=c("unassigned", "multiple", "unique"))
    df <- df[order(df$call),]
    
    ## assign colors to unique/multiple/unassigned
    cols <- c(multiple="indianred3", unassigned="grey", unique="blue") #"grey10")
    df$col <- cols[match(df$call, names(cols))]
    
    ## plot 
    p <- ggplot(df, aes(fraction_total_UMIs, log10(UMI_count), colour=call)) +
      scale_color_manual(values = cols) +
      geom_point(size=0.25) +
      xlab("% of total UMIs in cell") +
      ylab(expression('log'[10]*' UMIs in gRNA')) +
      geom_hline(yintercept = log10(min_count), lty=2, col="grey10") +
      geom_vline(xintercept = min_fraction, lty=2, col="grey10") 
  }
  
  ## produce sparse matrix of cell-gRNA pairs that are significant
  # subset to significant
  gRNA_calls_sig <- gRNA_calls[gRNA_calls$p_adj < pval_threshold,]
  # gRNA indices
  rows <- match(gRNA_calls_sig$gRNA, rownames(altExp(sce, 'CRISPR')))
  # cell indices
  cols <- match(gRNA_calls_sig$cell, colnames(altExp(sce, 'CRISPR')))
  # create sparse matrix
  gRNA_cell_matrix <- sparseMatrix(i = rows, j = cols, dims = dim(altExp(sce, 'CRISPR')))
  row.names(gRNA_cell_matrix) <- row.names(altExp(sce, 'CRISPR'))
  colnames(gRNA_cell_matrix) <- colnames(altExp(sce, 'CRISPR'))
  # add to sce obect, copying over col/rowData from CRISPR counts
  altExp(sce, 'gRNA_calls') <- altExp(sce, 'CRISPR')
  assay(altExp(sce, 'gRNA_calls')) <- gRNA_cell_matrix
  
  ## produce matrix of assignments
  gRNA_calls <- gRNA_calls %>% filter(order ==1) # keep only one entry per cell
  gRNA_calls <- gRNA_calls %>% mutate("gRNA" = NULL, 
                                      "UMI_count" = NULL, 
                                      "cell_total_UMIs" = NULL, 
                                      "fraction_total_UMIs" = NULL, 
                                      "order" = NULL, 
                                      "binom" = NULL, 
                                      "p_adj" = NULL) 
  if(!is.null(columns)){
    gRNA_calls <- gRNA_calls[, (columns):=NULL]
  }
  ## add to sce colData
  if(add_stats == FALSE){
    ## only add number of gRNAs called, and call annotation
    sel <- c("n_detected_gRNAs", "call", "detected_gRNAs")
    if(!is.null(columns)){
      sel <- c(sel, paste0('detected_', columns))
    }
    colData(sce) <- cbind(colData(sce), gRNA_calls[match(colnames(sce), gRNA_calls$cell), ..sel])
  }else{
    ## add everything
    colData(sce) <- cbind(colData(sce), gRNA_calls[match(colnames(sce), gRNA_calls$cell),-1])
  }
  # remove NAs produced by cells with _no_ gRNA counts at all (missing from analysis above)
  if(sum(is.na(sce$n_detected_gRNAs))>0){
    colData(sce[,is.na(sce$n_detected_gRNAs)])[,'call'] <- "unassigned"
    colData(sce[,is.na(sce$n_detected_gRNAs)])[,'n_detected_gRNAs'] <- 0
  }
  
  ## return results
  if(plot == TRUE){
    return(list(sce = sce,
                calls = gRNA_calls,
                plot = p, 
                dat_long = gRNA_long))
  }else{
    return(list(sce = sce,
                calls = gRNA_calls, 
                dat_long = gRNA_long))
  }
  
  # end of function
}

