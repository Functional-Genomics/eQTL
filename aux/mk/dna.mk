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

# DNA level (variants) provided as VCF files or as a matrix

step1: step0 $(step1b_dir)/complete

STEP1_TARGETS=
##############################
# start the analysis with VCFs
ifndef var_matrix

step1_a: tbi_files
step1_b: fix_vcf_headers
step1_c: filter_vcfs
step1_d: $(kpop_file) $(matched_expr_matrix) $(step1a_dir)/data_consistent

# include the samples in the expr matrix that are also in the genotype
$(matched_expr_matrix): $(dna_rna_mapfile) $(expr_matrix) $(gtf_eqtl_tsv)
	filter_phenotype.R $(expr_matrix) $(gtf_eqtl_tsv) $@.tmp && \
	pheno_preprocessing.py $(dna_rna_mapfile) $@.tmp $@.tmp2 &&  \
	mv $@.tmp2 $@ && rm -f $@.tmp 

# Fix headers if necessary
FIXEDHEADER_VCFS=$(foreach l,$(vcfs),$(name)/vcf/$(subst .vcf.gz,.fixedheader.vcf.gz,$(l)))


# required: .tbi for each vcf
TBI_FILES=$(subst .vcf.gz,.vcf.gz.tbi,$(FIXEDHEADER_VCFS))
tbi_files: $(TBI_FILES)

TARGETS1+=$(TBI_FILES)

fix_vcf_headers: $(FIXEDHEADER_VCFS) 


#bgzip
$(name)/vcf/%.fixedheader.vcf.gz: $(vcf_toplevel_dir)/%.vcf.gz
	mkdir -p $(@D) && $(fix_vcf_cmd) $< $@.tmp && mv $@.tmp $@

# Rules to process each chr
# chr=$(1) 
define make-rules-for-chr=
# TODO: remove from here
$(shell mkdir -p $(step1_dir)/$(1))
$(shell mkdir -p $(step1a_dir)/$(1))
$(shell mkdir -p $(step1b_dir)/$(1))

$(step1a_dir)/$(1)/chr$(1)_merged.vcf.gz: $$(foreach l,$(vcfs),$(step1_dir)/$(1)/$$(subst .vcf.gz,.filter.vcf.gz,$$(l))) $$(foreach l,$(vcfs),$(step1_dir)/$(1)/$$(subst .vcf.gz,.filter.vcf.gz.tbi,$$(l)))
	$$(file > $$@.lst,$$(foreach l,$(vcfs),$(step1_dir)/$(1)/$$(subst .vcf.gz,.filter.vcf.gz,$$(l))))  \
	sed -i -E "s/^ //;s/ +/\n/g" $$@.lst && \
	bcftools merge -l $$@.lst --output-type v | substitute_vcf.py | bgzip -c > $$@.tmp && \
	rm -f $$@.lst && \
	mv $$@.tmp $$@

# 
$(step1_dir)/$(1)/%.filter.vcf.gz: $(name)/vcf/%.fixedheader.vcf.gz
	mkdir -p $$(@D) && $(vcf_filter_cmd) -i $$<  -c $(1) -g $$(minGQ) -d $$(minRD) $(vcf_filter_option) | bgzip -c > $$@.tmp  && mv $$@.tmp $$@


##select MAC >=5, allow max missing sites in 20 percent of samples
$(step1a_dir)/$(1)/chr$(1)_merged.filt.vcf.gz: $(step1a_dir)/$(1)/chr$(1)_merged.vcf.gz
	vcftools --gzvcf $$< --max-missing $$(max_missing) --mac $$(mac)  --temp $$(@D) --recode --recode-INFO-all --out $$@.tmp && bgzip -c $$@.tmp.recode.vcf > $$@.tmp && mv $$@.tmp  $$@

$(step1a_dir)/$(1)/chr$(1).genotype.tsv: $(step1a_dir)/$(1)/chr$(1)_merged.filt.vcf.gz
	vcftools --012 --gzvcf $$< --out $$@.tmp &&
	generate_genotype_file $$@.tmp.012.indv $$@.tmp.012 $$@.tmp.012.pos  $$@.tmp &&	\
	mv $$@.tmp $$@


$(step1a_dir)/$(1)/chr$(1)_merged.filt.FILTER.summary: $(step1a_dir)/$(1)/chr$(1)_merged.vcf.gz
	mkdir -p $$(@D) && \
	vcftools --gzvcf $$< --max-missing $$(max_missing)  --temp $$(@D) --FILTER-summary --out $$@.tmp1 &&  \
	vcftools --gzvcf $$< --mac $$(mac)  --temp $$(@D) --FILTER-summary --out $$@.tmp2 && \
	vcftools --gzvcf $$< --mac $$(mac)  --max-missing $$(max_missing)  --temp $$(@D) --FILTER-summary --out $$@.tmp3 && \
	head -n 1 $$@.tmp1.FILTER.summary   | cut -f 1,2> $$@.tmp && \
	sum_tsv_col.sh $$@.tmp1.FILTER.summary 2 Missing  >> $$@.tmp && \
	sum_tsv_col.sh $$@.tmp2.FILTER.summary 2 Mac  >> $$@.tmp && \
	sum_tsv_col.sh $$@.tmp3.FILTER.summary 2 "Mac+Missing"  >> $$@.tmp && \
	mv $$@.tmp $$@

$(step1a_dir)/$(1)/chr$(1)_merged.filt.vcf.gz.plink: $(step1a_dir)/$(1)/chr$(1)_merged.filt.vcf.gz 
	mkdir -p $$(@D) && vcftools --gzvcf $$< --plink --out $$@.tmp && mv $$@.tmp.ped $$@.ped &&  mv $$@.tmp.map $$@.map && touch $$@


$(step1a_dir)/$(1)/plink_chr$(1).done: $(step1a_dir)/$(1)/chr$(1)_merged.filt.vcf.gz.plink
	mkdir -p $$(@D) && plink --file $$< --make-bed --fam --recode --noweb --out $$(subst .done,,$$@) && touch $$@

##convert genotypes into hdf5 file using limix binaries
##echo 'convert plink files into hdf5 files using limix binaries'
$(step1a_dir)/$(1)/chr$(1).hdf5: $(step1a_dir)/$(1)/plink_chr$(1).done 
	rm -f $$@.tmp $$@ && python $(LIMIX_BINARY)/limix_converter --outfile=$$@.tmp --plink=$$(subst .done,,$$<) && \
	geno_preprocessing.py $$@.tmp $(use_kinship_option) && \
	mv $$@.tmp $$@

ifneq ($(aggr_bed_file),)
$(step1b_dir)/$(1)/chr$(1)_filt.tsv: $$(foreach l,$(vcfs),$(step1_dir)/$(1)/$$(subst .vcf.gz,.filter.vcf.gz,$$(l))) $(aggr_bed_file)
	$$(file > $$@.lst,$$(foreach l,$(vcfs),$(step1_dir)/$(1)/$$(subst .vcf.gz,.filter.vcf.gz,$$(l))))  \
	cat $$@.lst | vcfs2matrix.sh $(aggr_bed_file) > $$@.tmp && mv $$@.tmp $$@
#&& rm -f $$@.lst

$(step1b_dir)/$(1)/chr$(1).genotype.tsv: $(step1b_dir)/$(1)/chr$(1)_filt.tsv
	cut -f 4,6- $$< |  geno_filtering.py $(geno_threshold) > $$@.tmp && mv $$@.tmp $$@

$(step1b_dir)/$(1)/chr$(1).hdf5: $(step1b_dir)/$(1)/chr$(1).genotype.tsv $(aggr_bed_file)
	cat $$< | generate_hdf5.py $(aggr_bed_file) $(1) $$@.tmp && mv $$@.tmp $$@
endif
endef

# Generate the rules per chr
$(foreach chr,$(geno_chr),$(eval $(call make-rules-for-chr,$(chr))))

STEP1_TARGETS=$(step1a_dir)/data_consistent $(samples_hdf5) $(matched_expr_matrix) filter2 filter1   $(foreach chr,$(geno_chr),$(step1b_dir)/$(chr)/chr$(chr).genotype.tsv)

$(step1a_dir)/data_consistent: $(dna_rna_mapfile) $(samples_hdf5)
#	check_consistency.py $(dna_rna_mapfile) $(kpop_file) $(expr_matrix) && touch $(step1_dir)/data_consistent
	check_consistency.py $(dna_rna_mapfile) $(samples_hdf5) && touch $(step1a_dir)/data_consistent

else
#############################################
# Variant matrix based analysis

step1_d: $(kpop_file) $(matched_expr_matrix) $(step1a_dir)/data_consistent


# filter the columns in the matrix based on their names
$(matched_expr_matrix): $(expr_matrix) $(var_matrix) 
	filter_columns.R $^ $@.tmp && mv $@.tmp $@ 

$(matched_var_matrix): $(var_matrix) $(expr_matrix)
	filter_columns.R $^ $@.tmp && mv $@.tmp $@ 


$(var_matrix).consistent: $(var_matrix) $(var_pos)
	geno_check_consistency.py $^ && touch $@

$(step1a_dir)/data_consistent: $(var_matrix).consistent $(matched_var_matrix) $(matched_expr_matrix)
	touch $(step1a_dir)/data_consistent

$(var_pos).bed4: $(var_pos)
	tail -n +2 $< | awk  'BEGIN {OFS="\t";} {print $$2,$$3,$$3,$$1;}'  >> $@.tmp && \
	mv $@.tmp $@

# var_min_freq [0,1]
$(matched_var_matrix).filt.tsv: $(matched_var_matrix) $(var_matrix).consistent
	cat $< | geno_filtering.py $(geno_threshold) > $@.tmp && mv $@.tmp $@

define make-rules-for-chr=
$(shell mkdir -p $(step1b_dir)/$(1))
$(step1b_dir)/$(1)/chr$(1).hdf5: $(matched_var_matrix).filt.tsv $(var_pos).bed4 
	cat $$< | generate_hdf5_mqtl.py  $(var_pos).bed4 $(1) $$@.tmp && mv $$@.tmp $$@

$(step1b_dir)/$(1)/chr$(1).genotype.tsv: $(matched_var_matrix).filt.tsv $(var_pos).bed4
	filter_geno_matrix.py $$^ $$@.tmp &&\
	mv $$@.tmp $$@

# deprecated
# $(step1b_dir)/$(1)/chr$(1).genotype.tsv: $(matched_var_matrix).filt.tsv $(var_pos).bed4
# 	grep "^$(1)\s" $(var_pos).bed4|cut -f 4 | sed -E 's/^/^/;s/$$$$/\\\s/' > $$@.tmp1 &&\
# 	head -n 1 $$< > $$@.tmp
# 	grep $$< -f $$@.tmp1 >> $$@.tmp 
# 	mv $$@.tmp $$@ && rm -f $$@.tmp1
endef


# Generate the rules per chr
$(foreach chr,$(geno_chr),$(eval $(call make-rules-for-chr,$(chr))))


#
STEP1_TARGETS+=$(var_matrix).consistent $(samples_hdf5) $(matched_expr_matrix) 

endif

# 
# build_Kpop.py kop.hdf5.tmp samples.hdf5 ...chr1.hdf5 ...chr2.hdf5 ...chr3.hdf5
$(kpop_file): $(foreach chr,$(geno_chr),$(step1b_dir)/$(chr)/chr$(chr).hdf5)
	build_Kpop.py $(kpop_file).tmp $(samples_hdf5).tmp $^ && \
	mv $(samples_hdf5).tmp $(samples_hdf5) &&\
	mv $(kpop_file).tmp $(kpop_file) && \
	sleep 1 && touch $(samples_hdf5)

$(samples_hdf5): $(kpop_file)
	if [ -e $@ ] ; then touch $(samples_hdf5); fi

$(step1_dir)/filter1/%.filter1.done: $(foreach chr,$(geno_chr),$(step1_dir)/$(chr)/$(subst .vcf.gz,.filter.vcf.gz.tbi,%.vcf.gz))	
	mkdir -p $(@D) && touch $@

FILTERED_VCF_FILES=$(foreach l,$(vcfs),$(step1_dir)/filter1/$(subst .vcf.gz,.filter1.done,$(l)))
FILTERED_VCF_FILES2=$(foreach chr,$(geno_chr),$(step1b_dir)/$(chr)/chr$(chr).hdf5)
#$(info $(FILTERED_VCF_FILES))
filter_vcfs:  $(FILTERED_VCF_FILES)
filter1: filter_vcfs
filter2: $(FILTERED_VCF_FILES2)

# one job per chr
TARGETS2+=$(FILTERED_VCF_FILES)
TARGETS3+=$(FILTERED_VCF_FILES2)

STEP1_TARGETS=$(step1a_dir)/data_consistent $(samples_hdf5) $(matched_expr_matrix) filter2 filter1   $(foreach chr,$(geno_chr),$(step1b_dir)/$(chr)/chr$(chr).genotype.tsv)
TARGETS4+=$(STEP1_TARGETS)




$(step1b_dir)/complete:  $(STEP1_TARGETS) vcf_stats
	$(call p_info,"Step 1 complete") touch $@

