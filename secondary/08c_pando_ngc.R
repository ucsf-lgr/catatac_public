#!/home/kfeng/anaconda3/envs/r42/bin/Rscript
# Load Packages
library(Pando)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

seurat_object <- readRDS"(path_to_integrated_sobj")
# seurat_object[["peaks"]] <- seurat_object[["ATAC"]]

# samples <- c("DASA_HIC2", "DASA_NT", "DMSO_NT")
# use all cell types 
Idents(seurat_object) <- seurat_object$treatment_and_guidecapture

data(motifs)

for(sample in samples){

}
seurat_object <- initiate_grn(seurat_object, rna_assay = 'SCT', peak_assay = 'ATAC')

print("GRN init'd")

seurat_object <- find_motifs(
    seurat_object,
    pfm = motifs,
    genome = BSgenome.Hsapiens.UCSC.hg38
)

print("Found Motifs")

seurat_object <- infer_grn(
    seurat_object,
    peak_to_gene_method = 'Signac',
    method = 'glm'
)

print("Inferred GRNs")

coef(seurat_object)

seurat_object <- find_modules(seurat_object)

print("Found Modules")

modules <- NetworkModules(seurat_object)

print("Network Modules Done")

saveRDS(seurat_object, "path_to_save_sobj")
saveRDS(modules, "path_to_save_modules")
