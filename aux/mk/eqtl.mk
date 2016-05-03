# =========================================================
# Copyright 2015-2016
#
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================

# complain if a file does not exist and exit
file_exists=$(if  $(realpath $(1)),,$(error $(1) not found))

############################################################
# Limix 
ifeq ($(eqtl_method),limix)

step4: $(step1b_dir)/complete $(step2_dir)/complete $(step3_dir)/complete $(eqtl_dir)/step4.complete

ifeq ($(cis_window),0)
############################################################
# Limix trans-eQTL (start)
$(info trans-EQTL mode)
# Trans specific
$(shell mkdir -p $(step1b_dir))
$(shell mkdir -p $(eqtl_dir)/all_chr/)
$(shell mkdir -p $(eqtl_dir)/all_chr/step1)
$(shell mkdir -p $(eqtl_dir)/all_chr/step2)
$(shell mkdir -p $(eqtl_dir)/all_chr/step3)
$(shell mkdir -p $(eqtl_dir)/all_chr/step4)

#
$(step1b_dir)/all_chr.hdf5: $(foreach chr,$(chromosomes),$(step1b_dir)/$(chr)/chr$(chr).hdf5)
	$(file >$@.lst.txt,$^) \
	sed -i -E "s/^ //;s/ +/\n/g" $@.lst.txt && \
	aggregate_chromosomes.py $@.lst.txt $@.tmp.hdf5 && \
	mv $@.tmp.hdf5 $@ && rm -f $@.lst.txt

# across the whole genome
All_QTL_JOBS=$(foreach j,$(shell seq $(n_folds)), $(eqtl_dir)/all_chr/step4/$(n_folds)_$(j).step4.tsv.gz)

$(eqtl_dir)/all_chr/step1/$(n_folds)_%.hdf5: $(step1b_dir)/all_chr.hdf5 $(step2_dir)/$(expr_matrix_filename).filtered.hdf5  $(kpop_file) $(step3_dir)/$(corr_method)/$(corr_method).hdf5 $(cov_sorted_hdf5)
	mkdir -p $(@D) && eqtl_trans.py $<   $(step2_dir)/$(expr_matrix_filename).filtered.hdf5  $(corr_method)  $(step3_dir)/$(corr_method)/$(corr_method).hdf5  $(kpop_file) $(cov_sorted_hdf5) $(limix_use_peer_covariates) $(cis_window) $(n_permutations) $(CHANGE_BETA_SIGN)  $(n_folds) $* $@.tmp && mv $@.tmp $@

$(eqtl_dir)/all_chr/step2/%.step2.tsv.gz $(eqtl_dir)/all_chr/step2/%.step2.tsv.gz.meta.tsv: $(eqtl_dir)/all_chr/step1/%.hdf5 $(step1b_dir)/all_chr.hdf5  $(cov_sorted_hdf5) $(step2_dir)/$(expr_matrix_filename).filtered.hdf5 $(step3_dir)/$(corr_method)/$(corr_method).hdf5
	eqtl_step2.py $(step1b_dir)/all_chr.hdf5  $(step2_dir)/$(expr_matrix_filename).filtered.hdf5  $(corr_method) $(step3_dir)/$(corr_method)/$(corr_method).hdf5  $(kpop_file) $(cov_sorted_hdf5) $(cis_window) $(n_permutations) $(snp_alpha)  $<  $@.tmp $@.meta && mv $@.meta $@.meta.tsv && mv $@.tmp $@

#$(eqtl_dir)/all_chr/step2/%.step2.tsv.gz.meta.tsv
$(eqtl_dir)/all_chr/step3/%.step3.tsv.gz: $(eqtl_dir)/all_chr/step2/%.step2.tsv.gz  
	eqtl_step3.py  $< $<.meta.tsv $@.tmp && mv $@.tmp $@

#$(eqtl_dir)/all_chr/step2/%.step2.tsv.gz.meta.tsv
$(eqtl_dir)/all_chr/step4/%.step4.tsv.gz: $(eqtl_dir)/all_chr/step3/%.step3.tsv.gz  
	eqtl_step4.py  $< $(eqtl_dir)/all_chr/step2/$*.step2.tsv.gz.meta.tsv $(fdr_threshold) $@.tmp && mv $@.tmp $@

$(eqtl_dir)/summary.tsv: $(All_QTL_JOBS)
	$(file >$@.lst.txt,$^) \
	sed -i -E "s/^ //;s/ +/\n/g" $@.lst.txt && \
	zcat $< | head -n 1 > $@.tmp &&\
	cat $@.lst.txt | while read n; do zcat $n | tail -n +2 >> $@.tmp; done &&\
	mv $@.tmp $@

TARGETS8+=$(eqtl_dir)/summary.tsv

# Limix trans-eQTL (end)
############################################################
else
############################################################
# Limix cis-eQTL (start)
$(info cis-EQTL mode)
eqtl_cmd1=eqtl_cis.py

# foreach chr
All_QTL_JOBS=$(foreach j,$(shell seq $(n_folds)),$(foreach chr,$(chromosomes), $(eqtl_dir)/$(chr)/$(n_folds)_$(j).step4.tsv.gz))

define QTL_JOBS_chr=
$(foreach j,$(shell seq $(n_folds)), $(eqtl_dir)/$(1)/$(n_folds)_$(j).step4.tsv.gz)
endef

# $(1) = chr
define make-qtl-rule-chr=
$(eqtl_dir)/$(1)/$(n_folds)_%.hdf5: $(step1b_dir)/$(1)/chr$(1).hdf5 $(step2_dir)/$(expr_matrix_filename).filtered.hdf5  $(kpop_file) $(step3_dir)/$(corr_method)/$(corr_method).hdf5 $(cov_sorted_hdf5)
	mkdir -p $$(@D) && $(eqtl_cmd1) $(step1b_dir)/$(1)/chr$(1).hdf5   $(step2_dir)/$(expr_matrix_filename).filtered.hdf5  $(corr_method)  $(step3_dir)/$(corr_method)/$(corr_method).hdf5  $(kpop_file) $(cov_sorted_hdf5) $(limix_use_peer_covariates) $(cis_window) $(n_permutations) $(CHANGE_BETA_SIGN)  $(n_folds) $$* $$@.tmp && mv $$@.tmp $$@

$(eqtl_dir)/$(1)/%.step2.tsv.gz $(eqtl_dir)/$(1)/%.step2.tsv.gz.meta.tsv: $(eqtl_dir)/$(1)/%.hdf5 $(step1b_dir)/$(1)/chr$(1).hdf5  $(cov_sorted_hdf5) $(step2_dir)/$(expr_matrix_filename).filtered.hdf5 $(step3_dir)/$(corr_method)/$(corr_method).hdf5
	eqtl_step2.py $(step1b_dir)/$(1)/chr$(1).hdf5  $(step2_dir)/$(expr_matrix_filename).filtered.hdf5  $(corr_method) $(step3_dir)/$(corr_method)/$(corr_method).hdf5  $(kpop_file) $(cov_sorted_hdf5) $(cis_window) $(n_permutations) $(snp_alpha)  $$<  $$@.tmp $$@.meta && mv $$@.meta $$@.meta.tsv && mv $$@.tmp $$@
endef

# %.step2.tsv.gz.meta.tsv
# $(eval $(call file_exists,%.step2.tsv.gz.meta.tsv))
%.step3.tsv.gz: %.step2.tsv.gz  
	eqtl_step3.py  $< $<.meta.tsv $@.tmp && mv $@.tmp $@

# %.step2.tsv.gz.meta.tsv
%.step4.tsv.gz: %.step3.tsv.gz  
	eqtl_step4.py  $< $*.step2.tsv.gz.meta.tsv $(fdr_threshold) $@.tmp && mv $@.tmp $@

$(foreach chr,$(chromosomes),$(eval $(call make-qtl-rule-chr,$(chr))))

# 
#
$(eqtl_dir)/summary.tsv: $(All_QTL_JOBS)
	$(file >$@.lst.txt,$^) \
	sed -i -E "s/^ //;s/ +/\n/g" $@.lst.txt && \
	zcat $< | head -n 1 > $@.tmp &&\
	cat $@.lst.txt | while read n; do zcat $n | tail -n +2 >> $@.tmp; done &&\
	mv $@.tmp $@

TARGETS8+=$(eqtl_dir)/summary.tsv
# Limix cis-eQTL (end)
############################################################
endif

#$(info $(All_QTL_JOBS))
TARGETS7+=$(All_QTL_JOBS)

# Limix 
############################################################	
else
############################################################
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
