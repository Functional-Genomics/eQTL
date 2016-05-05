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

step3: $(step3_dir)/complete

TARGETS6+= $(step3_dir)/$(corr_method)/$(corr_method).hdf5
# kpop - ktot
$(step3_dir)/panama/panama.hdf5: $(kpop_file) $(step2_dir)/$(expr_matrix_filename).filtered.hdf5
	mkdir -p $(@D) && \
	runpanama.py  $(kpop_file) $(step2_dir)/$(expr_matrix_filename).filtered.hdf5 $(hidden_k) $(snr_threshold) $@.tmp &&\
	mv $@.tmp  $@

$(step3_dir)/peer/peer.hdf5 $(step3_dir)/peer/peer_factors.hdf5: $(step2_dir)/$(expr_matrix_filename).filtered.hdf5 $(if $(subst y,,$(limix_use_peer_covariates)),,$(cov_sorted_hdf5))
	mkdir -p $(@D) && \
	run_peer.sh $<  $(hidden_k) $(peer_iterations)  $@.tmp $(@D)/peer_factors.hdf5 $(if $(subst y,,$(limix_use_peer_covariates)),,$(cov_sorted_hdf5)) && \
	mv $@.tmp  $@

$(step3_dir)/none/none.hdf5: $(kpop_file)
	mkdir -p $(@D) && \
	cp $< $@.tmp && \
	mv $@.tmp  $@

#$(step3_dir)/$(corr_method)/$(corr_method).hdf5
$(step3_dir)/complete: $(step1b_dir)/complete $(step2_dir)/complete   $(step3_dir)/$(corr_method)/$(corr_method).clus.png
	$(call p_info,"Step 3 complete") touch $@
