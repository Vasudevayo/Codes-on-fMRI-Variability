#!/bin/bash

#Written by Jeffrey Malins with guidance from W. Einar Mencl and Peter Molfese

base_path='PictureWordfMRI/BetaSeries'
subject_path='PictureWordfMRI/BetaSeries/Subject_IDs'
mask_path='PictureWordfMRI/BetaSeries/masks'

mkdir ${base_path}/Variability_Analyses
mkdir ${base_path}/Variability_Analyses/ROIs

#Convert MNI to TLRC
whereami -coord_file XYZ_MNI.1D'[0,1,2]' -lpi -calc_chain MNI TLRC -xform_xyz -coord_out XYZ_TLRC.1D

#Make ROIs Using Coordinates
echo "-49.500 21.498 2.609" | 3dUndump -prefix ${mask_path}/Martin2015_LIFGop_rad6+tlrc -srad 6 -orient LPI -master `@GetAfniBin`/TT_N27+tlrc.HEAD -xyz -
echo "-51.480 24.079 15.376" | 3dUndump -prefix ${mask_path}/Martin2015_LIFGpt_rad6+tlrc -srad 6 -orient LPI -master `@GetAfniBin`/TT_N27+tlrc.HEAD -xyz -
echo "-51.480 -58.716 -8.836" | 3dUndump -prefix ${mask_path}/Martin2015_LITG_rad6+tlrc -srad 6 -orient LPI -master `@GetAfniBin`/TT_N27+tlrc.HEAD -xyz -
echo "-55.440 -30.266 16.254" | 3dUndump -prefix ${mask_path}/Martin2015_LSTG_rad6+tlrc -srad 6 -orient LPI -master `@GetAfniBin`/TT_N27+tlrc.HEAD -xyz -
echo "-57.420 -25.441 -3.773" | 3dUndump -prefix ${mask_path}/Martin2015_LMTG_rad6+tlrc -srad 6 -orient LPI -master `@GetAfniBin`/TT_N27+tlrc.HEAD -xyz -
echo "-21.780 -47.986 50.305" | 3dUndump -prefix ${mask_path}/Martin2015_LSPL_rad6+tlrc -srad 6 -orient LPI -master `@GetAfniBin`/TT_N27+tlrc.HEAD -xyz -
echo "-11.880 -24.821 8.612" | 3dUndump -prefix ${mask_path}/Martin2015_LThal_rad6+tlrc -srad 6 -orient LPI -master `@GetAfniBin`/TT_N27+tlrc.HEAD -xyz -

#Resample masks
cd ${mask_path}
for aMask in Martin2015_R*
do 
	3dresample -dxyz 3 3 3 -prefix Resamp_${aMask} -inset ${aMask}
done
cd ${base_path}
	
#Extract values for ROIs
for aMask in LIFGop LIFGpt LITG LMTG LSPL LSTG LThal
do
	echo "Calculating Betas for ${aMask}"
	
	#Extract beta values
	for aSub in $(cat ${subject_path}/subjectlist.txt)
		do
 		echo "Calculating Betas for ${aSub}"
 		
		MeanA=$(3dROIstats -quiet -nzmean -nomeanout -mask ${mask_path}/Resamp_Martin2015_${aMask}_rad6+tlrc ${base_path}/Processed_Betas/AudMM_Correct/${aSub}_OutliersZeroed+tlrc.BRIK)
	   	echo $aSub $MeanA >> ${base_path}/Variability_Analyses/ROIs/${aMask}_Means_AudWordMM_Correct.txt
   		
 		MeanV=$(3dROIstats -quiet -nzmean -nomeanout -mask ${mask_path}/Resamp_Martin2015_${aMask}_rad6+tlrc ${base_path}/Processed_Betas/VisMM_Correct/${aSub}_OutliersZeroed+tlrc.BRIK)
   		echo $aSub $MeanV >> ${base_path}/Variability_Analyses/ROIs/${aMask}_Means_VisWordMM_Correct.txt
   		
   		OutliersA=$(3dROIstats -quiet -nzvoxels -nomeanout -mask ${mask_path}/Resamp_Martin2015_${aMask}_rad6+tlrc ${base_path}/Processed_Betas/AudMM_Correct/${aSub}_UpdatedOutlier_mask+tlrc.BRIK)
	   	echo $aSub $OutliersA >> ${base_path}/Variability_Analyses/ROIs/${aMask}_Outliers_AudWordMM_Correct.txt
   		
 		OutliersV=$(3dROIstats -quiet -nzvoxels -nomeanout -mask ${mask_path}/Resamp_Martin2015_${aMask}_rad6+tlrc ${base_path}/Processed_Betas/VisMM_Correct/${aSub}_UpdatedOutlier_mask+tlrc.BRIK)
   		echo $aSub $OutliersV >> ${base_path}/Variability_Analyses/ROIs/${aMask}_Outliers_VisWordMM_Correct.txt
   				
		done
done