#!/bin/bash

objFile=$1
labelFile=$2
outDir=$3

export PATH=./lib:${PATH}

matlab -nodisplay << EOF
	addpath('/data/noel/noel2/local/matlab/surfstat_chicago/');
	addpath('./lib')
	create_intrastruct_surfaces_spharm_laplacian('$1','$3','$2')
EOF
		
obj2vtk.sh $3_skelFinal.obj $3_skelFinal.vtk
obj2vtk.sh $3_surf.obj $3.vtk
