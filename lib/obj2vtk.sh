#!/bin/bash

print_help()
{
echo "
-----------------------------------------------
`basename $0` <in_file.obj> <in_file.vtk> [options]
-----------------------------------------------

Converts an obj (surface) to vtk keeping all world coordinates.
The function needs Matlab to run, and its path should have
SurfStat and the function save_surface_vtk.m to work properly. 
The function itself will load EMMA into matlab's path, and
in the Noel lab the function save_volume_vtk.m is in the SVN,
so it shoul also be available without you having anything to do.
If the function complains of things not being in its path, 
add a line at the bottom of this function to specify where things are.

Update October 2009: We can now write data per vertex! see Options.

Options:
-data <data.txt> <data_name>
  You can use as many -data switches as you want!
  It is extremely important that each -data switch is followed by 
  two arguments: 
    data.txt  : a txt file of values per vertex. This is the kind
                of file you would use for brain-view to color your surface.
    data_name : a string for baptizing the data in the vtk file. 
                It will help you keep track of what each piece of data is.

Example:
`basename $0` in_file.obj in_file.vtk \\
  -data thickness.txt thickness \\
  -data curvature.txt curvature



See also: mnc2vtk

Luis Concha
Noel Lab
April 2009
lconcha@gmail.com

-----------------------------------------------
`basename $0` <in_file.obj> <in_file.vtk> 
-----------------------------------------------

"
}


if [ $# -lt 2 ] 
then
   echo "Must specify at least 2 arguments."
	print_help
	exit 1
fi

declare -i i
i=1
nData=0
tmpFile=/tmp/obj2vtk_tmp_$$.txt
rm -f $tmpFile
for arg in "$@"
do
   case "$arg" in
	-help)
		print_help
		exit 1
	;;
	-data)
		nData=`expr $nData + 1`
		nextarg=`expr $i + 1`
		nextarg2=`expr $i + 2`
	   eval datum=\${${nextarg}}
      eval datum_name=\${${nextarg2}}
      echo "${datum_name},${datum}" >> $tmpFile
   ;;
	esac
  i=`expr $i + 1`
done



################
# BEGIN FUNCTION
###############

obj=$1
vtk=$2

# Prepare a matlab job and run it.
job=/tmp/$$.job
echo "addpath /data/noel/noel2/colliot/local/linux/imaging/install/matlab/emma" > $job
echo "[ surf, ab ] = SurfStatReadSurf1( '$obj' );" >> $job



cat $tmpFile | while read line
do
  varName=`echo $line | awk -F, '{print $1}'`
  varFilename=`echo $line | awk -F, '{print $2}'`
  echo "valuesPerVertex.${varName} = load('${varFilename}');" >> $job
  i=`expr $i + 1`
done

if [ $nData -gt 0 ]
then
   echo "save_surface_vtk(surf,'$vtk','BINARY', valuesPerVertex);" >> $job
else
   echo "save_surface_vtk(surf,'$vtk','ASCII');" >> $job
fi

cat $job

matlab -nodisplay < ${job}

# Clean up
rm -f $job
rm -f $tmpFile
