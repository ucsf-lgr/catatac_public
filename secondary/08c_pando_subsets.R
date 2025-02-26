#!/home/kfeng/anaconda3/envs/r42/bin/Rscript
# Load Packages
library(Pando)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

# load obj
seurat_object <- readRDS("path_to_integrated_sobj")

# list of treatment_and_guidecapture of interest
samples <- c("DASA_HIC2", "DASA_NT", "DMSO_NT")

# set idents
Idents(seurat_object) <- seurat_object$treatment_and_guidecapture

data(motifs)

for(sample in samples){
    cat("Running:", sample, "\n")
    sobj <- subset(seurat_object, idents = sample)
    sobj <- initiate_grn(sobj, rna_assay = 'SCT', peak_assay = 'ATAC')

    print("GRN init'd")

    sobj <- find_motifs(
        sobj,
        pfm = motifs,
        genome = BSgenome.Hsapiens.UCSC.hg38
    )

    print("Found Motifs")

    sobj <- infer_grn(
        sobj,
        peak_to_gene_method = 'Signac',
        method = 'glm'
    )

    print("Inferred GRNs")

    coef(sobj)

    sobj <- find_modules(sobj)

    print("Found Modules")

    modules <- NetworkModules(sobj)

    print("Network Modules Done")

    saveRDS(sobj, paste0("path_to_save_sobj"))
    saveRDS(modules, paste0("path_to_a_folder", sample, "_pando_module.rds" ))
}


