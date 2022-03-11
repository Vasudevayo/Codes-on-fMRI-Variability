#!/bin/bash

#Written by Jeffrey Malins and Ryan Staples, with guidance from W. Einar Mencl and Peter Molfese

input_base='PictureWordfMRI'
subject_path='PictureWordfMRI/BetaSeries/Subject_IDs'

for subj in $(cat ${subject_path}/subjectlist.txt)
do

	cd ${input_base}/${subj}/stim_times
	
	#Note: dmUBLOCK was used for print and speech match and mismatch conditions
	#GAM was used for picture trials only (as there were no responses)
	n_stims=$(ls -1q ${subj}* | wc -l)
	stim_list=$(for n in $(seq 1 ${n_stims}); do printf "\'dmUBLOCK\' "; done)
	stim_list_add=$(printf "\'GAM\'")
	stim_list_final="${stim_list}${stim_list_add}"

	cd ${input_base}/${subj}
	
	#This code dynamically builds a list of condition names
	#This was variable across participants, as some did not have incorrect trials for certain types
	#All speech and mismatch conditions are listed; the name of the picture condition is added in "regress_stim_times" below
	cond_names=$(cat stim_times/labels_parsed.txt)
	no_pic=${cond_names::${#cond_names}-8}
	cond_list=$(for cond in ${no_pic}; do echo "stim_times/*${cond}*"; done)
	
	#Single subject processing in AFNI
	afni_proc.py -subj_id ${subj} -script script.${subj}.IM_RT.tcsh \
	-out_dir ${subj}.IM_RT \
	-dsets func/sn*_?_*.nii.gz func/sn*_??_*.nii.gz \
	-blocks tshift align tlrc volreg blur mask scale regress \
	-copy_anat anat/*MPRAGE*.nii.gz \
	-anat_has_skull yes \
	-tcat_remove_first_trs 4 \
	-volreg_align_e2a -volreg_tlrc_warp \
	-tlrc_opts_at -init_xform AUTO_CENTER \
	-align_opts_aea -giant_move \
	-tshift_opts_ts -tpattern alt+z2 \
	-blur_size 8 -volreg_align_to first  \
	-regress_stim_times ${cond_list} \
		stim_times/pictures.txt \
	-regress_stim_labels ${cond_names} \
	-regress_stim_types IM \
	-regress_local_times \
    -regress_basis_multi ${stim_list_final} \
	-regress_censor_outliers 0.1 \
	-regress_censor_motion 0.3 \
	-remove_preproc_files \
	-regress_opts_3dD \
	-allzero_OK \
	-GOFORIT 666 \
	-jobs 8 \
	-bash -execute
	
	done
done