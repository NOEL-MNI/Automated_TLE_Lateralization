#!/bin/bash

## This script extract features for one structure and for one side. The inputs are the following:
# 1. Population Prefix
# 2. Individual ID
# 3. Structure being sampled
# 4. Sampled side
# 5. Blade directory
# 6. Image directory where the images to sample are located
# 7. Directory for txt files containing sampled data

prefix=$1
id=$2
side=$3           	# L or R
subfield=$4		# CA, SUB or DG
blade_directory=$5
image_directory=$6
feature_directory=$7

# Check if the feature directory exists
if [ ! -d ${feature_directory} ]; then
	mkdir -v ${feature_directory}
fi

out_surface=${blade_directory}/${prefix}_${id}_${side}_${subfield}_surf.obj
blade_surface=${blade_directory}/${prefix}_${id}_${side}_${subfield}_skelFinal.obj

echo ${out_surface}
echo ${blade_surface}

## Get Columnar Volume
columnar_volume_file=${feature_directory}/${prefix}_${id}_${side}_${subfield}_ColVol.txt

matlab -nodisplay << EOF
	addpath('/data/noel/noel2/local/matlab/surfstat_chicago/');
	addpath('./lib')
	vertex_colvol = computeColumnarVolume('${out_surface}','${blade_surface}');
	SurfStatWriteData('${columnar_volume_file}',vertex_colvol);
EOF

## Get T2 Signal

t2_image=${image_directory}/${prefix}_${id}_t2cor-0.4_stx_norm-0.4.mnc
t2_file=${feature_directory}/${prefix}_${id}_${side}_${subfield}_nnt2.txt

# Segment the t2 file and get the mean signal in the ventricle
matlab -nodisplay << EOF
	addpath('/data/noel/noel2/local/matlab/surfstat_chicago/');
	addpath('./lib')
	getmeanventricleT2('${t2_image}',['${image_directory}','/','${prefix}','_','${id}']);
EOF

# Intersect t2 and medial surface
volume_object_evaluate \
	-nearest_neighbour \
	${t2_image} \
	${blade_surface} \
	${t2_file}


## Get FLAIR/T1w Signal
ratio_image=${image_directory}/${prefix}_${id}_t2wt1wratio_final.mnc
ratio_file=${feature_directory}/${prefix}_${id}_${side}_${subfield}_t2wt1wratio.txt

volume_object_evaluate \
	${ratio_image} \
	${blade_surface} \
	${ratio_file}
