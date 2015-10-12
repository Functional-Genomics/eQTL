#!/bin/bash
eqtl_pipeline vcf_toplevel_dir=test/ex1 dna_rna_mapfile=test/ex1/map_file_test.tsv conf=t.conf name=aaaa expr_matrix=test/ex1/pheno.tsv vcfs="test.HG00096.vcf.gz test.HG00097.vcf.gz test.HG00099.vcf.gz test.HG00100.vcf.gz test.HG00101.vcf.gz" gtf_file=test/ex1/gencode.v19.annotation.gtf.gz chromosomes="1 2" n_folds=10 mac=0 minRD=0 minGQ=1 min_expr=0 min_perc_samples=0 cov_matrix=test/ex1/cov_test.tsv $*

