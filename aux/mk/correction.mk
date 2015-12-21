step3: $(step3_dir)/complete

TARGETS6+= $(step3_dir)/$(corr_method)/$(corr_method).hdf5
# kpop - ktot
$(step3_dir)/panama/panama.hdf5: $(kpop_file) $(step2_dir)/$(expr_matrix_filename).filtered.hdf5
	mkdir -p $(@D) && \
	runpanama.py  $(kpop_file) $(step2_dir)/$(expr_matrix_filename).filtered.hdf5 $(hidden_k) $(snr_threshold) $@.tmp &&\
	mv $@.tmp  $@

$(step3_dir)/peer/peer.hdf5 $(step3_dir)/peer/peer_factors.hdf5: $(step2_dir)/$(expr_matrix_filename).filtered.hdf5 $(if $(subst y,,$(limix_use_peer_covariates)),,$(cov_sorted_hdf5))
	mkdir -p $(@D) && \
	runpeer.py $<  $(hidden_k) $(peer_iterations)  $@.tmp $(@D)/peer_factors.hdf5 $(if $(subst y,,$(limix_use_peer_covariates)),,$(cov_sorted_hdf5)) && \
	mv $@.tmp  $@

$(step3_dir)/none/none.hdf5: $(kpop_file)
	mkdir -p $(@D) && \
	cp $< $@.tmp && \
	mv $@.tmp  $@

#$(step3_dir)/$(corr_method)/$(corr_method).hdf5
$(step3_dir)/complete: $(step1a_dir)/complete $(step2_dir)/complete   $(step3_dir)/$(corr_method)/$(corr_method).clus.png
	$(call p_info,"Step 3 complete") touch $@
