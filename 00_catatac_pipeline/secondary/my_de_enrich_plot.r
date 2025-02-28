

MyDEenrichRPlot <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  balanced = TRUE,
  logfc.threshold = 0.25,
  assay = NULL,
  max.genes,
  test.use = 'wilcox',
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = NULL,
  num.pathway = 10,
  return.gene.list = FALSE,
  ...
) {
  enrichr.installed <- PackageCheck("enrichR", error = FALSE)
  if (!enrichr.installed[1]) {
    stop(
      "Please install the enrichR package to use DEenrichRPlot",
      "\nThis can be accomplished with the following command: ",
      "\n----------------------------------------",
      "\ninstall.packages('enrichR')",
      "\n----------------------------------------",
      call. = FALSE
    )
  }
  if (is.null(x = enrich.database)) {
    stop("Please specify the name of enrichR database to use")
  }
  if (!is.numeric(x = max.genes)) {
    stop("please set max.genes")
  }
  assay <- assay %||% DefaultAssay(object = object)

  DefaultAssay(object = object) <- assay

  all.markers <- FindMarkers(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    only.pos = FALSE,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    assay = assay
  )

  pos.markers <- all.markers[all.markers[, 2] > logfc.threshold & all.markers[, 1] < p.val.cutoff, , drop = FALSE]

  if(nrow(pos.markers) == 0){
    message("No positive markers pass the logfc.thershold")
    pos.er <- c()
  }

  else{
  pos.markers.list <- rownames(x = pos.markers)[1:min(max.genes, nrow(x = pos.markers))]
  pos.er <- enrichR::enrichr(genes = pos.markers.list, databases = enrich.database)
  pos.er <- do.call(what = cbind, args = pos.er)
  pos.er$log10pval <- -log10(x = pos.er[, paste(enrich.database, sep = ".", "P.value")])
  pos.er$term <- pos.er[, paste(enrich.database, sep = ".", "Term")]
  pos.er <- pos.er[1:num.pathway, ]
  pos.er$term <- factor(x = pos.er$term, levels = pos.er$term[order(pos.er$log10pval)])
  gene.list <- list(pos = pos.er)
  }

  if (isTRUE(x = balanced)) {
    neg.markers <- all.markers[all.markers[, 2] < logfc.threshold & all.markers[, 1] < p.val.cutoff, , drop = FALSE]
    neg.markers.list <- rownames(x = neg.markers)[1:min(max.genes, nrow(x = neg.markers))]
    neg.er <- enrichR::enrichr(genes = neg.markers.list, databases = enrich.database)
    neg.er <- do.call(what = cbind, args = neg.er)
    neg.er$log10pval <- -log10(x = neg.er[, paste(enrich.database, sep = ".", "P.value")])
    neg.er$term <- neg.er[, paste(enrich.database, sep = ".", "Term")]
    neg.er <- neg.er[1:num.pathway, ]
    neg.er$term <- factor(x = neg.er$term, levels = neg.er$term[order(neg.er$log10pval)])

      if(isTRUE(length(neg.er$term) == 0) & isTRUE(length(pos.er == 0))){
        stop("No positive or negative marker genes identified")
      }

      else{
        if(isTRUE(length(neg.er$term) == 0)){

        gene.list <- list(pos = pos.er)

        }
        else{
          gene.list <- list(pos = pos.er, neg = neg.er)
        }
      }

  }
  if (return.gene.list) {
    return(gene.list)
  }

  if(nrow(pos.markers) == 0){
    message("No positive markers to plot")

    if (isTRUE(x = balanced)) {

      p2 <- ggplot(data = neg.er, aes_string(x = "term", y = "log10pval")) +
        geom_bar(stat = "identity", fill = "indianred2") +
        coord_flip() + xlab("Pathway") +
        scale_fill_manual(values = cols, drop = FALSE) +
        ylab("-log10(pval)") +
        ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) +
        theme_classic() +
        geom_text(aes_string(label = "term", y = 0),
                  size = 5,
                  color = "black",
                  position = position_dodge(1),
                  hjust = 0)+
        theme(axis.title.y= element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
      p <- p2

    }
    else{
      stop("Nothing to plot")
    }
  }

  else {
  p <- ggplot(data = pos.er, aes_string(x = "term", y = "log10pval")) +
    geom_bar(stat = "identity", fill = "dodgerblue") +
    coord_flip() + xlab("Pathway") +
    scale_fill_manual(values = cols, drop = FALSE) +
    ylab("-log10(pval)") +
    ggtitle(paste(enrich.database, ident.1, sep = "_", "positive markers")) +
    theme_classic() +
    geom_text(aes_string(label = "term", y = 0),
              size = 5,
              color = "black",
              position = position_dodge(1),
              hjust = 0)+
    theme(axis.title.y= element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  if (isTRUE(x = balanced)) {

    p2 <- ggplot(data = neg.er, aes_string(x = "term", y = "log10pval")) +
      geom_bar(stat = "identity", fill = "indianred2") +
      coord_flip() + xlab("Pathway") +
      scale_fill_manual(values = cols, drop = FALSE) +
      ylab("-log10(pval)") +
      ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) +
      theme_classic() +
      geom_text(aes_string(label = "term", y = 0),
                size = 5,
                color = "black",
                position = position_dodge(1),
                hjust = 0)+
      theme(axis.title.y= element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    p <- p+p2

  }
  }

  return(p)
}