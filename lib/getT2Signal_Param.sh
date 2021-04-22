#!/bin/bash

struct=$1
inDir="/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades"
outDir="/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades/BladeSampling/T2Signal_Nuc"

inDir=${inDir}/${struct}/Blades

for fullID in `ls ${inDir}`; do

	id=`echo ${fullID} | cut -d '_' -f1-2`
	echo ${id}

	t2File=/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades/BladeSampling/T2Vol_Nuc/TLE_${id}_t2cor-0.4_stx_norm-0.4.mnc
	
	
	if [ ! -e ${outDir}/TLE_${fullID}_nnt2.txt ]; then
		if [ -e ${t2File} ]; then
		
			echo ${fullID}
			volume_object_evaluate \
				-nearest_neighbour \
				${t2File} \
				${inDir}/${fullID}/TLE_${fullID}_skelFinal.obj \
				${outDir}/TLE_${fullID}_nnt2.txt
	
		elif [ -e ${t2File%.gz} ]; then
		
			echo ${fullID}
			volume_object_evaluate \
				-nearest_neighbour \
				${t2File%.gz} \
				${inDir}/${fullID}/TLE_${fullID}_skelFinal.obj \
				${outDir}/TLE_${fullID}_nnt2.txt
		else
			echo "No T2 for case ${id}"
		fi
	fi

done
