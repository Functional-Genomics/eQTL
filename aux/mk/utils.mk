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

#Version and license info
pname=eqtlXXXXX
version=0.0.1
contact="Add contact"
license=This pipeline is distributed  under the terms of the GNU General Public License 3

$(info *****************************************************)
$(info * $(pname) $(version))
$(info * $(contact))
$(info * $(license))
$(info *)
$(info * Initializing...)

################################################################################
# Auxiliary functions
################################################################################

# Information messages
define p_info=
$(info $(shell date "+%H:%M:%S %d/%m/%Y * ") $(1))
endef

# Error messages
define p_error=
$(info $(shell date "+%H:%M:%S %d/%m/%Y") * ERROR: $(1)) && $(error Fatal error)
endef


define get_mac=
$(shell bash -c "echo \($(words $(vcfs)) \* $(maf) \* 2 +1\)/1 | bc")
endef	

# complain if a file does not exist and exit
file_exists=$(if  $(realpath $(1)),,$(call p_error,$(1) not found))

#  check if a variable  $(1) is defined - return the variable name if it is defined or empty otherwise
is_defined=$(if $(subst undefined,,$(origin $(1))),$(1),)

##################################################################################
################################################################################
# Generic file extension rules

# 
%.vcf.gz.tbi: %.vcf.gz
	tabix -p vcf $< || ( rm -f $@ && exit 1)

%.vcf.tbi: %.vcf
	tabix -p vcf $<  || ( rm -f $@ && exit 1)

%.gtf: %.gtf.gz
	gunzip -c $< > $@.tmp && mv $@.tmp $@
#

%.tsv: %.tsv.gz
	gunzip -c $< > $@.tmp && mv $@.tmp $@


#
###############################################
# Load configuration (mandatory)
# use a configuration file?
ifdef conf
 $(call file_exists,$(conf))
 $(info * Trying to load configuration file $(conf)...)
 include $(conf)
 $(info * Configuration loaded.)
else
 $(call p_error,Configuration file missing)
endif




#
TARGETS0=
TARGETS1=
TARGETS2=
TARGETS3=
TARGETS4=
TARGETS5=
TARGETS6=
TARGETS7=
TARGETS8=
TARGETS9=
