ifeq ($(eqtl_method),limix)
step4: $(step1b_dir)/complete $(step2_dir)/complete $(step3_dir)/complete $(eqtl_dir)/step4.complete

ifeq ($(cis_window),0)
############################################################
# Limix trans-eQTL (start)
$(info trans-EQTL mode)
# Trans specific
$(step1b_dir)/all_chr.hdf5: $(foreach chr,$(chromosomes),$(step1a_dir)/$(chr)/chr$(chr).hdf5)
	$(file >$@.lst.txt,$^) \
	sed -i -E "s/^ //;s/ +/\n/g" $@.lst.txt && \
	aggregate_chromosomes.py $@.lst.txt $@.tmp.hdf5 && \
	mv $@.tmp.hdf5 $@ && rm -f $@.lst.txt

# across the whole genome
All_QTL_JOBS=$(foreach j,$(shell seq $(n_folds)), $(eqtl_dir)/all_chr/$(n_folds)_$(j).hdf5)


$(eqtl_dir)/all_chr/$(n_folds)_%.hdf5: $(step1b_dir)/all_chr.hdf5 $(step2_dir)/$(expr_matrix_filename).filtered.hdf5  $(kpop_file) $(step3_dir)/$(corr_method)/$(corr_method).hdf5 $(cov_sorted_hdf5)
	mkdir -p $(@D) && eqtl_trans.py $<   $(step2_dir)/$(expr_matrix_filename).filtered.hdf5  $(corr_method)  $(step3_dir)/$(corr_method)/$(corr_method).hdf5  $(kpop_file) $(cov_sorted_hdf5) $(limix_use_peer_covariates) $(cis_window) $(n_permutations) $(n_folds) $* $@.tmp && mv $@.tmp $@


$(eqtl_dir)/all_chr/summary.hdf5: $(All_QTL_JOBS)  $(cov_sorted_hdf5) $(step2_dir)/$(expr_matrix_filename).filtered.hdf5 $(step3_dir)/$(corr_method)/$(corr_method).hdf5
	$(file >$@.lst.txt,$(All_QTL_JOBS)) \
	sed -i -E "s/^ //;s/ +/\n/g" $@.lst.txt && \
	eqtl_summary.py $(step1b_dir)/all_chr.hdf5  $(step2_dir)/$(expr_matrix_filename).filtered.hdf5  $(corr_method) $(step3_dir)/$(corr_method)/$(corr_method).hdf5  $(kpop_file) $(cov_sorted_hdf5) $(cis_window) $(n_permutations) $(snp_alpha) $(n_folds) $@.lst.txt  $@.tmp && \
	rm -f $@.lst.txt && \
	mv $@.tmp $@

$(eqtl_dir)/summary.hdf5: $(eqtl_dir)/all_chr/summary.hdf5
	$(file >$@.lst.txt,$^) \
	sed -i -E "s/^ //;s/ +/\n/g" $@.lst.txt && \
	eqtl_aggregate.py $@.tmp $(cis_window) $(n_permutations) $@.lst.txt && \
	mv $@.tmp $@ && rm -f $@.lst.txt

# Limix trans-eQTL (end)
############################################################
else
############################################################
# Limix cis-eQTL (start)
$(info cis-EQTL mode)
eqtl_cmd1=eqtl_cis.py


# foreach chr
All_QTL_JOBS=$(foreach j,$(shell seq $(n_folds)),$(foreach chr,$(chromosomes), $(eqtl_dir)/$(chr)/$(n_folds)_$(j).hdf5))

define QTL_JOBS_chr=
$(foreach j,$(shell seq $(n_folds)), $(eqtl_dir)/$(1)/$(n_folds)_$(j).hdf5)
endef


# $(1) = chr
define make-qtl-rule-chr=
$(eqtl_dir)/$(1)/$(n_folds)_%.hdf5: $(step1b_dir)/$(1)/chr$(1).hdf5 $(step2_dir)/$(expr_matrix_filename).filtered.hdf5  $(kpop_file) $(step3_dir)/$(corr_method)/$(corr_method).hdf5 $(cov_sorted_hdf5)
	mkdir -p $$(@D) && $(eqtl_cmd1) $(step1b_dir)/$(1)/chr$(1).hdf5   $(step2_dir)/$(expr_matrix_filename).filtered.hdf5  $(corr_method)  $(step3_dir)/$(corr_method)/$(corr_method).hdf5  $(kpop_file) $(cov_sorted_hdf5) $(limix_use_peer_covariates) $(cis_window) $(n_permutations) $(n_folds) $$* $$@.tmp && mv $$@.tmp $$@

# $(step3_dir)/$(1)/summary.hdf5
$(eqtl_dir)/$(1).hdf5: $(call QTL_JOBS_chr,$(1))  $(cov_sorted_hdf5) $(step2_dir)/$(expr_matrix_filename).filtered.hdf5 $(step3_dir)/$(corr_method)/$(corr_method).hdf5
	$$(file >$$@.lst.txt,$(call QTL_JOBS_chr,$(1))) \
	sed -i -E "s/^ //;s/ +/\n/g" $$@.lst.txt && \
	eqtl_summary.py $(step1b_dir)/$(1)/chr$(1).hdf5  $(step2_dir)/$(expr_matrix_filename).filtered.hdf5  $(corr_method) $(step3_dir)/$(corr_method)/$(corr_method).hdf5  $(kpop_file) $(cov_sorted_hdf5) $(cis_window) $(n_permutations) $(snp_alpha) $(n_folds) $$@.lst.txt  $$@.tmp && \
	rm -f $$@.lst.txt && \
	mv $$@.tmp $$@

endef


$(foreach chr,$(chromosomes),$(eval $(call make-qtl-rule-chr,$(chr))))

$(eqtl_dir)/summary.hdf5: $(foreach chr,$(chromosomes),$(eqtl_dir)/$(chr).hdf5)
	$(file >$@.lst.txt,$^) \
	sed -i -E "s/^ //;s/ +/\n/g" $@.lst.txt && \
	eqtl_aggregate.py $@.tmp $(cis_window) $(n_permutations) $@.lst.txt && \
	mv $@.tmp $@ && rm -f $@.lst.txt

TARGETS8+=$(foreach chr,$(chromosomes),$(eqtl_dir)/$(chr).hdf5)
# Limix cis-eQTL (end)
############################################################
endif

#$(info $(All_QTL_JOBS))
TARGETS7+=$(All_QTL_JOBS)


$(eqtl_dir)/summary.tsv: $(eqtl_dir)/summary.hdf5
	get_results.py $(fdr_threshold) $@.tmp $^ && mv $@.tmp $@

else
#######################################################
# FastQTL
ifeq ($(eqtl_method),fastqtl)

ifeq ($(corr_method),none)
define make-fastqtl-rule-chr=
$(eqtl_dir)/$(1).tsv: $(step1b_dir)/$(1)/chr$(1)_merged.filt.vcf.gz $(step2_dir)/$(expr_matrix_filename).filtered.bed.gz.tbi $(step1b_dir)/$(1)/chr$(1)_merged.filt.vcf.gz.tbi
	mkdir -p $$(@D) && fastqtl --vcf $$< --bed $(step2_dir)/$(expr_matrix_filename).filtered.bed.gz --window $(cis_window) --threshold  $(fdr_threshold) --permute $(n_permutations) --region "$(1):1-100000000000"   --out $$@.tmp && mv $$@.tmp $$@
endef

else
# correction method used
define make-fastqtl-rule-chr=
$(eqtl_dir)/$(1).tsv: $(step1b_dir)/$(1)/chr$(1)_merged.filt.vcf.gz $(step3_dir)/$(corr_method)/$(corr_method).bed.gz.tbi $(step1b_dir)/$(1)/chr$(1)_merged.filt.vcf.gz.tbi
	mkdir -p $$(@D) && fastqtl --vcf $$< --bed $(step3_dir)/$(corr_method)/$(corr_method).bed.gz --window $(cis_window) --threshold $(fdr_threshold) --permute $(n_permutations) --region "$(1):1-100000000000"   --out $$@.tmp && mv $$@.tmp $$@
endef
endif

$(foreach chr,$(chromosomes),$(eval $(call make-fastqtl-rule-chr,$(chr))))


#TODO nf: complete
$(eqtl_dir)/summary.tsv: $(foreach chr,$(chromosomes),$(eqtl_dir)/$(chr).tsv)
	head -n 1 $< > $@.tmp && \
	tail -q -n +2 $^  >> $@.tmp && mv $@.tmp $@


TARGETS7+=$(foreach chr,$(chromosomes),$(eqtl_dir)/$(chr).tsv)
TARGETS8+=$(eqtl_dir)/summary.tsv

else
#######################################################
# matrixEQTL

# 1 - chr
# TODO: generate a pos file
ifeq ($(corr_method),none)
define make-meqtl-rule-chr=
$(eqtl_dir)/$(1).tsv: $(step1b_dir)/$(1)/chr$(1).genotype.tsv $(step2_dir)/$(expr_matrix_filename).filtered.tsv $(cov_sorted_hdf5).tsv $(gtf_eqtl_tsv)
	mkdir -p $$(@D) && run_matrix_eqtl $$< $(step2_dir)/$(expr_matrix_filename).filtered.tsv  $(cov_sorted_hdf5).tsv $(gtf_eqtl_tsv) $(cis_window) $(fdr_threshold) $$@.tmp && rename ".tmp" "" $$@.tmp*
endef

else

# meqQTL

define make-meqtl-rule-chr=
$(eqtl_dir)/$(1).tsv: $(step1b_dir)/$(1)/chr$(1).genotype.tsv $(step3_dir)/$(corr_method)/$(corr_method).tsv $(cov_sorted_hdf5).tsv $(gtf_eqtl_tsv)
	mkdir -p $$(@D) && run_matrix_eqtl $$< $(step3_dir)/$(corr_method)/$(corr_method).tsv  $(cov_sorted_hdf5).tsv $(gtf_eqtl_tsv) $(cis_window) $(fdr_threshold) $$@.tmp && rename ".tmp" "" $$@.tmp*
endef
endif

$(foreach chr,$(chromosomes),$(eval $(call make-meqtl-rule-chr,$(chr))))

# merge all files into one
$(eqtl_dir)/summary.tsv: $(foreach chr,$(chromosomes),$(eqtl_dir)/$(chr).tsv)
	head -n 1 $< > $@.tmp && \
	tail -q -n +2 $^ | grep -v "No significant" >> $@.tmp && mv $@.tmp $@


TARGETS7+=$(foreach chr,$(chromosomes),$(eqtl_dir)/$(chr).tsv)
TARGETS8+=$(eqtl_dir)/summary.tsv

endif
endif

step4: $(step1b_dir)/complete $(step2_dir)/complete $(step3_dir)/complete $(eqtl_dir)/step4.complete report

$(eqtl_dir)/step4.complete:  $(eqtl_dir)/summary.tsv
	$(call p_info,"Step 4 complete") touch $@

TARGETS9+=$(step1b_dir)/complete $(step2_dir)/complete $(step3_dir)/complete $(step3_dir)/step4.complete report

phony_targets+= setup setup_files

###################################################
