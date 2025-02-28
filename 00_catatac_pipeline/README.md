# catatac
CAT-ATAC primary and secondary workflows


# Primary
- 00_bcl2fq.sh
= 01_fastqc.sh
= 02_reads2tsv.sh
= 03_select_reads_by_guide.sh or 03_select_reads_by_guide_with_mismatch.sh
- 04_run_count.sh

## To run primary:
1. Git clone dir.
2. Install bbtools/bbmap and add the path.
3. Modifiy the PRIMARY_ANALYSIS_TEMPLATE.sh to contain the appropriate paths.
4. Exectute PRIMARY_ANALYSIS_TEMPLATE.sh.

# Secondary
assign_guides_to_cells_00.ipynb
seurat.ipynb
match_guide2barcode_04.ipynb
demux_guides_00.ipynb

