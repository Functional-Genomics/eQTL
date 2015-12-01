# get a matrix from a specific .hdf5 file
# matrix/dataset kpop,ktot,phenotype
# row/col are switched?
$(cov_sorted_hdf5).tsv:  $(cov_sorted_hdf5)
	hdf52tsv $< "/covariates" "/row_header/sample_ID" "-"   $@.tmp y && mv $@.tmp $@

$(step2_dir)/$(expr_matrix_filename).filtered.tsv: $(step2_dir)/$(expr_matrix_filename).filtered.hdf5
	hdf52tsv $< "phenotype/Ytransformed" "phenotype/row_header/sample_ID" "/phenotype/col_header/phenotype_ID"   $@.tmp y && mv $@.tmp $@

$(step3_dir)/none/none.tsv: $(step3_dir)/none/none.hdf5
	hdf52tsv $< "Kpop" "/row_header/sample_ID" "/col_header/sample_ID"  $@.tmp  n && mv $@.tmp $@

$(step3_dir)/peer/peer.tsv: $(step3_dir)/peer/peer.hdf5
	hdf52tsv $< "/phenotype" "/row_header/sample_ID" "/col_header/phenotype_ID" $@.tmp y && mv $@.tmp $@

$(step3_dir)/panama/panama.tsv: $(step3_dir)/panama/panama.hdf5
	hdf52tsv $< "Ktot" "/row_header/sample_ID" "/col_header/sample_ID" $@.tmp y && mv $@.tmp $@


%.clus.png: %.tsv $(sample2class_file)
	generate_clustering  $< $(sample2class_file) $@.tmp && mv $@.tmp $@

%.pca.png %.pca_13.png: %.tsv $(sample2class_file)
	generate_pca  $< $(sample2class_file) $@ > $@.txt 


#########################################################################
report: plots vcf_stats $(report_dir)/settings.tsv


plots: $(report_dir)/plots

$(report_dir)/settings.tsv: $(conf) 
	mkdir -p $(@D) &&\
	( $(foreach v,$(def_vars), echo $v  $($v);)) | tr " " "\t" > $@.tmp && mv $@.tmp $@

# Copy the plots and tsv file to the report folder
$(report_dir)/plots:  $(report_dir)/expr_filtered_clus.png $(report_dir)/expr_filtered_corrected_clus.png $(report_dir)/expr_filtered_pca.png $(report_dir)/expr_filtered_corrected_pca.png $(report_dir)/vcf_filtering.png

$(report_dir)/expr_filtered_clus.png:  $(step2_dir)/$(expr_matrix_filename).filtered.clus.png $(step2_dir)/$(expr_matrix_filename).filtered.tsv
	mkdir -p $(@D) && cp $^ $(@D) && cp $< $@

$(report_dir)/expr_filtered_pca.png:  $(step2_dir)/$(expr_matrix_filename).filtered.pca.png $(step2_dir)/$(expr_matrix_filename).filtered.pca_13.png $(step2_dir)/$(expr_matrix_filename).filtered.tsv
	mkdir -p $(@D) && cp $^ $(@D) && cp $< $@

$(report_dir)/expr_filtered_corrected_clus.png: $(step3_dir)/$(corr_method)/$(corr_method).clus.png $(step3_dir)/$(corr_method)/$(corr_method).tsv
	mkdir -p $(@D) && cp $^ $(@D) && cp $< $@

$(report_dir)/expr_filtered_corrected_pca.png: $(step3_dir)/$(corr_method)/$(corr_method).pca.png $(step3_dir)/$(corr_method)/$(corr_method).pca_13.png $(step3_dir)/$(corr_method)/$(corr_method).tsv
	mkdir -p $(@D) && cp $^ $(@D) && cp $< $@


###################################################
# filtering summary stats

########
# $(1) = chr  e.g, 1, 2, 3, ...
define make-vcf-stats-for-chr=

%.vcf.gz.chr$(1).summary: %.vcf.gz  %.vcf.gz.tbi
	bcftools stats -r "$(1)" $$< > $$@.tmp && mv $$@.tmp $$@

%.vcf.gz.chr$(1).snps: %.vcf.gz.chr$(1).summary
	echo -n "$(1) " | tr " " "\t" > $$@.tmp  && grep "records:" $$< | head -n1 | cut -f 4 >> $$@.tmp && mv $$@.tmp $$@
endef

$(foreach chr,$(chromosomes),$(eval $(call make-vcf-stats-for-chr,$(chr))))

%.vcf.gz.snps: $(foreach chr,$(chromosomes),%.vcf.gz.chr$(chr).snps)
	mkdir -p $(@D) && \
	echo "Chr $(notdir $*)" | sed -E "s/\s+/\t/g" > $@.tmp.col1 &&\
	cat $@.tmp.col1  $^ >$@ && rm -f $@.tmp $@.tmp.col1
# vcf file was already split by chr 
$(step1_dir)/%.vcf.gz.chr.snps: $(foreach chr,$(chromosomes),$(step1_dir)/$(chr)/%.vcf.gz.chr$(chr).snps)
	mkdir -p $(@D) && \
	echo "Chr $(notdir $*)" | sed -E "s/\s+/\t/g" > $@.tmp.col1 &&\
	cat $@.tmp.col1  $^ >$@ && rm -f $@.tmp $@.tmp.col1


VCF_STATS_0=$(foreach vcf,$(vcfs),$(name)/vcf/$(subst .vcf.gz,.fixedheader.vcf.gz,$(vcf)).snps)

$(report_dir)/vcf_snps_0.tsv: $(VCF_STATS_0)
	mkdir -p $(@D) && \
	$(file > $@.tmp,Chr $(vcfs)) $(file > $@.tmp.lst,$^ ) \
	sed -i "s/ /\t/g" $@.tmp && \
	cat $@.tmp.lst | mjoin -stdin | tail -n +2 | tr " " "\t">> $@.tmp  && mv $@.tmp $@

#
VCF_STATS_1=$(foreach vcf,$(vcfs),$(step1_dir)/$(subst .vcf.gz,.filter.vcf.gz,$(vcf)).chr.snps)

$(report_dir)/vcf_snps_1.tsv: $(VCF_STATS_1)
	mkdir -p $(@D) && \
	$(file > $@.tmp,Chr $(vcfs))  	$(file > $@.tmp.lst,$^ ) \
	sed -i "s/ /\t/g" $@.tmp && \
	cat $@.tmp.lst | mjoin -stdin  | tail -n +2 |tr " " "\t">> $@.tmp  && mv $@.tmp $@

VCF_STATS_2=$(foreach c,$(chromosomes),$(step1a_dir)/$(c)/chr$(c)_merged.filt.FILTER.summary)	
$(report_dir)/vcf_snps_2.tsv: $(VCF_STATS_2)
	mkdir -p $(@D) && \
	echo Chr $(chromosomes) |tr " " "\t" > $@.tmp &&\
	mjoin $^ | tail -n +2 |sed -E "s/\s+/\t/g;s/\s$$//">> $@.tmp &&\
	mv  $@.tmp $@

TARGETS4+=$(VCF_STATS_0) $(VCF_STATS_1) $(VCF_STATS_2)

vcf_stats: $(report_dir)/vcf_snps_2.tsv $(report_dir)/vcf_snps_1.tsv $(report_dir)/vcf_snps_0.tsv 


$(report_dir)/vcf_filtering.png: $(report_dir)/vcf_snps_0.tsv $(report_dir)/vcf_snps_1.tsv $(report_dir)/vcf_snps_2.tsv  
	get_barplot.py $^ $@.tmp && mv $@.tmp $@


target_reports_vcfs1=$(VCF_STATS_0) $(VCF_STATS_1) $(VCF_STATS_2)

report_targets1:
	echo $(target_reports_vcfs1)
