library(enrichR)
library(crayon)

# Function: plot_enrichment
# Finds upregulated differentially expressed genes in specified cluster compared to
# another cluster, by default, the iPSC.
# Returns a list. First element is the patchworked plots. The second is the df of
# the list of hits.
# Access them like this:
#   enrich_ret_vals[[1]]
#   enrich_ret_vals[[2]][1]

plot_enrichment <- function(
        seurat_obj,
        auroc_cutoff = 0.70,
        cluster_id,
        ipsc_cluster_id = 'iPSC',
        enrich.databases = c(
            'GO_Biological_Process_2021',
            'GO_Molecular_Function_2021', 
            'GO_Cellular_Component_2021',
            'PanglaoDB_Augmented_2021', 
            'Descartes_Cell_Types_and_Tissue_2021',
            'Human_Gene_Atlas',
            'Tabula_Sapiens',
            'Azimuth_Cell_Types_2021'
        ),
        n_max_hits = 20,
        de_padj_cutoff = 1e-5,
        enrichment_pval_cutoff = 1e-2
    ) {

    cols = 4
    rows = 2
    cat(blue("Displaying top ",n_max_hits, " hits with adj_pval <= ", enrichment_pval_cutoff, "\n"))

    markers_rna <- presto:::wilcoxauc.Seurat(
        X = seurat_obj, 
        group_by = 'celltype', 
        assay = 'data', 
        groups_use = c(cluster_id, ipsc_cluster_id),
        seurat_assay = 'SCT'
    )

    df_pos_markers = dplyr::filter(
        markers_rna, 
        auc >= auroc_cutoff,
        group == cluster_id, 
        padj <= de_padj_cutoff, 
        logFC > 0
        ) %>% arrange(-auc)

    pos.markers.list = df_pos_markers$feature
    if(length(pos.markers.list) == 0){
        message("No positive markers pass the logfc.thershold")
        df_pos_enriched <- c()
    }
    plots <- c()
    enrichr_dbs <- c()
    i = 1
    for(enrich.database in enrich.databases) {
        df_pos_enriched <- enrichR::enrichr(
            genes = pos.markers.list, 
            databases = enrich.database
        )
        df_pos_enriched <- do.call(what = cbind, args = df_pos_enriched)
        df_pos_enriched$log10pval <- -log10(
            x = df_pos_enriched[, paste(enrich.database, sep = ".", "Adjusted.P.value")]
        )
        n_markers <- nrow(df_pos_markers)
        df_pos_enriched$term    <- df_pos_enriched[, paste(enrich.database, sep = ".", "Term")]
        df_pos_enriched$overlap <- df_pos_enriched[, paste(enrich.database, sep = ".", "Overlap")]
        df_pos_enriched$bar_label <- paste0(df_pos_enriched$term, "   ", df_pos_enriched$overlap)

        select_high_conf_hits <- df_pos_enriched$log10pval >= -log10(enrichment_pval_cutoff)
        #cat(red("1", nrow(df_pos_enriched), " rows in ", enrich.database, "\n"))
        df_pos_enriched <- df_pos_enriched[select_high_conf_hits, ]
        n_max_rows      <- min(n_max_hits, nrow(df_pos_enriched))
        #cat(red("2", nrow(df_pos_enriched), " rows in ", enrich.database, "n_max_rows=",n_max_rows,"\n"))
        if(n_max_rows > 0) {
            df_pos_enriched <- df_pos_enriched[1:n_max_rows, ]
        }

        if(nrow(df_pos_enriched) == 0) {            
            cat(yellow("\nSkipping", enrich.database, ". No hits with adj pval < ", enrichment_pval_cutoff, "\n")) 
            title_text <- paste0(enrich.database, ": No hits with adj pval < ", enrichment_pval_cutoff)
            p <- ggplot(data = df_pos_enriched, aes_string(x = "term", y = "log10pval")) +
                geom_blank() + 
                coord_flip() +
                labs(title = title_text,
                        size = 6,
                        color = "black",
                        position = position_dodge(1),
                        hjust = 0)

            plots[[i]] = p
            enrichr_dbs[[i]] = df_pos_enriched
            i = i + 1        
        } else {
            df_pos_enriched$term <- factor(x = df_pos_enriched$term, levels = df_pos_enriched$term[order(df_pos_enriched$log10pval)])
            gene.list <- list(pos = df_pos_enriched)

            p <- ggplot(data = df_pos_enriched, aes_string(x = "term", y = "log10pval")) +
                geom_bar(stat = "identity", fill = "dodgerblue") +
                coord_flip() + xlab("Pathway") +
                scale_fill_manual(values = cols, drop = FALSE) +
                ylab("-log10(adj.pval)") +
                ggtitle(paste(enrich.database, cluster_id, ": ", n_markers, "positive markers", sep = " ")) +
                theme_classic() +
                geom_text(aes_string(label = "bar_label", y = 0),
                        size = 6,
                        color = "black",
                        position = position_dodge(1),
                        hjust = 0)+
                theme(axis.title.y= element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),                 
                    axis.text=element_text(size=12),
                    axis.title=element_text(size=12,face="bold")
                )

            plots[[i]] = p
            enrichr_dbs[[i]] = df_pos_enriched
            i = i + 1        
        }
    }
    cat("\n", n_markers, "markers found for", cluster_id, "\n")
    p <- patchwork::wrap_plots(plots, nrow = rows, ncol = cols)
    list(p, enrichr_dbs)
}