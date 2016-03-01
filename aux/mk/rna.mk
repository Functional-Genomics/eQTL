
# TODO: add checks for the variables

step2: $(step2_dir)/complete

# Filtering the pheno matrix based on fpkm threshold and number of samples
$(step2_dir)/$(expr_matrix_filename).filtered.tsv: $(matched_expr_matrix)
	filtering_pheno.py $< $(min_expr) $(min_perc_samples) $@.tmp && mv $@.tmp $@

# Quantile normalization per classes/studies/groups
$(step2_dir)/$(expr_matrix_filename).filtered.qn.tsv: $(step2_dir)/$(expr_matrix_filename).filtered.tsv
	irap_qn -i $< -m $(sample2class_file) -o $@.tmp && mv $@.tmp $@

# Transform expression values of each gene across all the samples.
$(step2_dir)/$(expr_matrix_filename).filtered.qn.trans.tsv: $(step2_dir)/$(expr_matrix_filename).filtered.qn.tsv
	normalise_pheno.py $< $(expr_transform) $@.tmp && mv $@.tmp $@


$(step2_dir)/$(expr_matrix_filename).filtered.hdf5: $(step2_dir)/$(expr_matrix_filename).filtered.qn.trans.tsv $(gtf_eqtl_tsv)  $(samples_hdf5) $(cov_hdf5)
	$(LIMIX_BINARY)/limix_converter --outfile=$@.tmp --csv=$< && \
	hdf_annotation.py $(gtf_eqtl_tsv) $@.tmp && \
	cp $(cov_hdf5) $(cov_sorted_hdf5).tmp && \
	sort_ids.py $@.tmp $(cov_sorted_hdf5).tmp $(samples_hdf5) &&\
	mv $(cov_sorted_hdf5).tmp $(cov_sorted_hdf5) && \
	mv $@.tmp $@ 

# $(step2_dir)/$(expr_matrix_filename).filtered.hdf5: $(step2_dir)/$(expr_matrix_filename).limix.hdf5  $(samples_hdf5) $(cov_hdf5)
# 	cp $(cov_hdf5) $(cov_sorted_hdf5).tmp && \
# 	filtering_pheno.py $< $(min_expr) $(min_perc_samples) $(expr_transform) $@.tmp && \
# 	sort_ids.py $@.tmp $(cov_sorted_hdf5).tmp $(samples_hdf5) &&\
# 	mv $(cov_sorted_hdf5).tmp $(cov_sorted_hdf5) && \
# 	mv $@.tmp $@ 
# $(step2_dir)/$(expr_matrix_filename).limix.hdf5: $(matched_expr_matrix) $(gtf_eqtl_tsv)
# 	$(LIMIX_BINARY)/limix_converter --outfile=$@.tmp --csv=$< && \
# 	hdf_annotation.py $(gtf_eqtl_tsv) $@.tmp && \
# 	mv $@.tmp $@ 


$(cov_sorted_hdf5): $(step2_dir)/$(expr_matrix_filename).filtered.hdf5
	if [ -e $@ ] ; then sleep 1; touch $@; fi

# $(step2_dir)/$(expr_matrix_filename).filtered.hdf5
$(step2_dir)/complete:  $(cov_sorted_hdf5) $(step2_dir)/$(expr_matrix_filename).filtered.clus.png
	$(call p_info,"Step 2 complete") touch $@
