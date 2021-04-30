# Code repository to generate results for:
> B Caldairou, N Foit, et al. "Automated TLE Lateralization". Neurology 2021 (under revision)
<hr>


### 1. Software Dependencies

Software dependencies include:
 - The minc toolkit: https://bic-mni.github.io/
 - `CMRep`: https://github.com/pyushkevich/cmrep
 - `QHull`: http://www.qhull.org/download/
 - `SurfStat` MATLAB Toolbox: https://www.math.mcgill.ca/keith/surfstat/

`CMRep` and `QHull` are used together to generate the skeleton from the subfield labels.

`SurfStat` is used to read minc volumes, surfaces and .txt feature files in `MATLAB`.

Pre-processing (non-uniformity correction, registration to MNI space, inter-modality co-registration) can be performed using `FreeSurfer`, `FSL` or any other MRI processing pipeline.

Subfield segmentation can be performed through `FreeSurfer` or `ASHS`. Best alternative to `SurfPatch` is to use `ASHS` with our publicly available manual subfield labels (http://fcon_1000.projects.nitrc.org/indi/retro/mni_hipposeg.html).

Outer surface extraction and parameterization can be performed through the SPHARM extension of `3D Slicer`. 
Binaries are available at: https://www.nitrc.org/projects/spharm-pdm

Main parameters are:
- Subdivision level: `32 for CA, 24 for SUB and DG`
- SPHARM degree: `32 for CA, 24 for SUB and DG`
- Surface templates are provided (they follow the manual segmentation protocol described in: Kulaga-Yoskovitz, J. et al., Scientific Data, 2015)
- We also provide SPHARM orientation example in order to orient spherical parameters in a similar manner than ours.
<hr>


### 2. Directory Organization

The assumed organization of the directories is specified below:

	- One directory to store the surfaces and the resulting skeletons
		(One directory per case, per side and per subfield)
		[surface_directory]
		|
		|-- TLE_0362_1_L_CA.obj
		|-- TLE_0362_1_L_CA_inter.obj
		|-- TLE_0362_1_L_CA_skelFinal.obj
		|
		|-- TLE_0362_1_L_SUB.obj
		|-- TLE_0362_1_L_SUB_inter.obj
		|-- TLE_0362_1_L_SUB_skelFinal.obj
	
	- One directory to store feature images 
		[image_directory]
		|
		|-- TLE_0362_1_t2cor-0.4_stx_norm-0.4.mnc
		|-- TLE_0362_1_t2wt1wratio_final.mnc
		|-- TLE_0362_1_Brainclass.mnc
		|-- TLE_0362_1_T2norm.mat
		|-- TLE_0362_1_Ventricle.mnc
	
	- One directory to store the txt files resulting from the
		intersection between the blades and the images
		[feature_directory]
		|
		|-- TLE_0362_1_L_CA_ColVol.txt
		|-- TLE_0362_1_R_CA_ColVol.txt
		|-- TLE_0362_1_L_CA_nnt2.txt
		|-- TLE_0362_1_R_CA_nnt2.txt
		|-- TLE_0362_1_L_CA_t2wt1wratio.txt
		|-- TLE_0362_1_R_CA_t2wt1wratio.txt
<hr>

### 3. Blade Extraction

The script is used in the following way (we recommend using the nomenclature stated above):
```
./get_medial_surface \
	${subfield_surface} \
	${subfield_volume} \
	${individual_output_prefix} # Example: ${SurfaceDirectory}/TLE_0362_1_L_CA
```
<hr>


### 4. Feature Extraction

This script intersects blade surfaces with volumetric images and computes the columnar volume:
```
./get_individual_features.sh \
	${prefix} \ # e.g. TLE
	${id} \     # e.g. 0362_1
	${subfield} # CA, SUB or DG
	${side}     # L or R
	${surface_directory} \ # The same described in directory organization
	${image_directory} \
	${feature_directory}
```
<hr>


### 5. Lateralization

The following series of script can be executed in order to accomplish various tasks:

|   						|   										|
|---------------------------|-------------------------------------------|
|A_GetTrainingData.m	| loads training data and control data set |
|B_StatisticalStudy.m	| performs a statistical study between controls and training data set |
|C_Validation_TrainingSet.m	| performs the nested repeated 5-Fold validation based on the training set |
|D_Validation_TestSet.m	| performs the repeated validation on a separated the test set for generalizability |
|E_Train_Model.m	| performs the training of a final model |
|F_Test_Individual.m	| loads saved training models and ROI and performs individual lateralization |
<hr>


### License
<a href= "https://opensource.org/licenses/BSD-3-Clause"><img src="https://img.shields.io/badge/License-BSD%203--Clause-blue.svg" /></a>

```console
Copyright 2021 Neuroimaging of Epilepsy Laboratory, McGill University
```