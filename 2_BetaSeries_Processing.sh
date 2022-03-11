#!/bin/bash

#Written by Jeffrey Malins, Ryan Staples, and Anish Kurian, with guidance from W. Einar Mencl and Peter Molfese

#This script will extract, mask, and remove/replace outliers with zero for the beta series of a given stimulus type.

stim_tag='VisMM_Correct'
stim_label='VisMM_Correct'

input_base='PictureWordfMRI'
base_path='PictureWordfMRI/BetaSeries'
subject_path='${base_path}/Subject_IDs'
mask_dataset='${base_path}/masks/TT_3mm_mask+tlrc'

MaskedBeta_path='${base_path}/Masked_Betas/${stim_label}'
ProcessedBeta_path='${base_path}/Processed_Betas/${stim_label}'

#Make directories
mkdir ${base_path}/Extracted_Betas
mkdir ${base_path}/Extracted_Betas/${stim_label}
  
mkdir ${base_path}/Masked_Betas
mkdir ${base_path}/Masked_Betas/${stim_label}
 
mkdir ${base_path}/Processed_Betas
mkdir ${base_path}/Processed_Betas/${stim_label}

mkdir ${base_path}/QC

#Process beta series
for subj in $(cat ${subject_path}/subjectlist.txt)
do
if [ ! -f ${base_path}/Extracted_Betas/${stim_label}/Betas_${subj}+tlrc.HEAD ]; then
(	echo "Extracting Beta Series for ${subj}"
    	for briknum in {1..40}
    	do
     		# Note: This will give you errors, but it still works.
      		3dTcat -prefix ${base_path}/Extracted_Betas/${stim_label}/Betas_${subj}+tlrc ${input_base}/${subj}/${subj}.IM_RT/stats.${subj}+tlrc"[${stim_tag}_RT#0_Coef]"
      		3dTcat -glueto ${base_path}/Extracted_Betas/${stim_label}/Betas_${subj}+tlrc ${input_base}/${subj}/${subj}.IM_RT/stats.${subj}+tlrc"[${stim_tag}_RT#${briknum}_Coef]"
    	done
    	
	#Calculate number of sub-bricks per subject after removing invalid trials
  	NumBricksCorrect=$(3dinfo -nv ${base_path}/Extracted_Betas/${stim_label}/Betas_${subj}+tlrc)
  	echo $subj $NumBricksCorrect >> ${base_path}/QC/NumBricksCorrect_${stim_label}.txt

	#Remove censored trials (volumes with all beta values equal to zero)
	briktotal=$(3dinfo -nv ${base_path}/Extracted_Betas/${stim_label}/Betas_${subj}+tlrc)
   	brikend="$(($briktotal-1))"
	for numbrik in $(seq 0 1 ${brikend})
 	do
 		MeanBeta=$(3dBrickStat ${base_path}/Extracted_Betas/${stim_label}/Betas_${subj}+tlrc[${numbrik}])
 		MeanBeta_int=$(echo "($MeanBeta+.5)/1" | bc)
 		if [ $MeanBeta_int -eq 0 ]; then
 			echo "Censoring volume" $MeanBeta >> ${base_path}/QC/CensoredVolumes/CensoredVolumes_${subj}_${stim_label}.txt
 		else 
 			echo "Including volume" $MeanBeta >> ${base_path}/QC/CensoredVolumes/CensoredVolumes_${subj}_${stim_label}.txt
 			3dbucket -aglueto ${base_path}/Processed_Betas/${stim_label}/Betas_${subj}_CensoredRemoved+tlrc ${base_path}/Extracted_Betas/${stim_label}/Betas_${subj}+tlrc[${numbrik}]
 		fi
 	done

	#Calculate number of sub-bricks per subject after removing invalid trials and censored trials for motion and outliers
 	NumBricksUncensored=$(3dinfo -nv ${base_path}/Processed_Betas/${stim_label}/Betas_${subj}_CensoredRemoved+tlrc)
 	echo $subj $NumBricksUncensored >> ${base_path}/QC/NumBricksAfterCensoring_${stim_label}.txt

	#Apply a mask to extracted beta values
 	echo "Masking ${subj}"
 	3dcalc -a ${base_path}/Processed_Betas/${stim_label}/Betas_${subj}_CensoredRemoved+tlrc -b ${mask_dataset} -expr 'a*b' -prefix ${base_path}/Masked_Betas/${stim_label}/MaskedBetas_${subj}+tlrc

	#Identify outliers and relevant metrics and create outlier mask
	#echo "Identifying outliers for ${subj}"
	3dToutcount -mask ${mask_dataset} -save ${ProcessedBeta_path}/${subj}_Outliers ${MaskedBeta_path}/MaskedBetas_${subj}+tlrc
 	3dToutcount -mask ${mask_dataset} -range -fraction ${MaskedBeta_path}/MaskedBetas_${subj}+tlrc >> ${ProcessedBeta_path}/${subj}_Fraction.txt
	3dcalc -a ${ProcessedBeta_path}/${subj}_Outliers+tlrc -expr 'ispositive(a)' -prefix ${ProcessedBeta_path}/${subj}_Outlier_mask+tlrc
 	
 	#Remove trials with too many outliers (i.e., greater than 10% of voxels)
 	censbriktotal=$(3dinfo -nv ${MaskedBeta_path}/MaskedBetas_${subj}+tlrc)
  	censbrikend="$(($censbriktotal-1))"
  	for numbrik in $(seq 0 1 ${censbrikend})
 	do
 		NumOutliers=$(3dBrickStat -mask ${mask_dataset} -mean ${ProcessedBeta_path}/${subj}_Outlier_mask+tlrc[${numbrik}])
 		NumOutliers_int=$(echo "($NumOutliers*100+.5)/1" | bc)
 		if [ $NumOutliers_int -gt 10 ]; then
 			echo "Censoring volume" $NumOutliers >> ${base_path}/QC/CensoredVolumes/OutlierVols_${subj}_${stim_label}.txt
 		else 
 			echo "Including volume" $NumOutliers >> ${base_path}/QC/CensoredVolumes/OutlierVols_${subj}_${stim_label}.txt
 			3dbucket -aglueto ${ProcessedBeta_path}/${subj}_OutlierVolsRemoved+tlrc ${MaskedBeta_path}/MaskedBetas_${subj}+tlrc[${numbrik}]
 			3dbucket -aglueto ${ProcessedBeta_path}/${subj}_UpdatedOutlier_mask+tlrc ${ProcessedBeta_path}/${subj}_Outlier_mask+tlrc[${numbrik}]
 		fi
 	done
 
 	#Calculate number of sub-bricks per subject after removing invalid trials and censored trials for motion and outliers
 	NumBricksAfterOutliers=$(3dinfo -nv ${ProcessedBeta_path}/${subj}_OutlierVolsRemoved+tlrc)
 	echo $subj $NumBricksAfterOutliers >> ${base_path}/QC/NumBricksAfterOutliers_${stim_label}.txt
  	
  	#Create mask of outlier and non-outlier voxels
	3dcalc -a ${ProcessedBeta_path}/${subj}_UpdatedOutlier_mask+tlrc -expr 'iszero(a)' -prefix ${ProcessedBeta_path}/${subj}_NonOutlier_mask+tlrc	
 	
 	#Multiply outlier/non-outlier mask by original dataset
 	3dcalc -a ${ProcessedBeta_path}/${subj}_NonOutlier_mask+tlrc -b ${ProcessedBeta_path}/${subj}_OutlierVolsRemoved+tlrc -expr 'a*b' -prefix ${ProcessedBeta_path}/${subj}_OutliersZeroed+tlrc
	
)	
fi
done
