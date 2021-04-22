#!/bin/bash

struct=$1
inDir="/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades"
outDir="/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades/BladeSampling/T2T1Ratio_Tal_Nuc"

inDir=${inDir}/${struct}/Blades

for fullID in `ls ${inDir}`; do

	id=`echo ${fullID} | cut -d '_' -f1-2`
	


	t2File=/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades/BladeSampling/T2T1Ratio_Images_Tal_Nuc/TLE_${id}_t2wt1wratio_final.mnc.gz
	
	if [ ! -e ${outDir}/TLE_${fullID}_t2wt1wratio.txt ]; then
	
		if [ -e ${t2File} ]; then
		
			echo ${fullID}
			volume_object_evaluate \
				${t2File} \
				${inDir}/${fullID}/TLE_${fullID}_skelFinal.obj \
				${outDir}/TLE_${fullID}_t2wt1wratio.txt
		
		elif [ -e ${t2File%.gz} ]; then
		
			echo ${fullID}
			volume_object_evaluate \
				${t2File%.gz} \
				${inDir}/${fullID}/TLE_${fullID}_skelFinal.obj \
				${outDir}/TLE_${fullID}_t2wt1wratio.txt
		else
			echo "No T2 for case ${id}"
		fi
	
	fi

done
