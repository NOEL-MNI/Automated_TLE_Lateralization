THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. 
EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. 
THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

THESE ARE THE CODES USED TO GENERATE RESULTS FOR THE FOLLOWING PAPER:
B.CALDAIROU, N.FOIT, et al. "AUTOMATED TLE LATERALIZATION"....

YOU MIGHT HAVE TO ADAPT THE CODE TO YOUR NEEDS.


I. SOFTWARE DEPENDENCIES

Softwares dependencies include:
 - The minc toolkit: https://bic-mni.github.io/
 - CMRep: https://github.com/pyushkevich/cmrep
 - QHull: http://www.qhull.org/download/
 - SurfStat Matlab Toolbox: https://www.math.mcgill.ca/keith/surfstat/

CMRep and QHull are used together to generate the skeleton from the subfield labels.

SurfStat is used to read minc volumes, surfaces and .txt feature files in matlab.

Pre-processing (non-uniformity correction, registration to MNI space, inter-modality co-registration) can be performed using freeSurfer, fsl or any other processing pipeline.

Subfield segmentation can be performed through freeSurfer or ANTS.

Outer surface extraction and parameterization can be performed through the SPHARM extension of 3D Slicer. 
Binaries are also available at: https://www.nitrc.org/projects/spharm-pdm. 
Main parameters are:
	- Subdivision level: 32 for CA, 24 for SUB and DG
	- SPHARM degree: 32 for CA, 24 for SUB and DG
	- Surface templates are provided (they follow the manual segmentation protocol described in: Kulaga-Yoskovitz, J. et al., Scientific Data, 2015)
	- We also provide SPHARM orientation example in order to correctly orient the 

II. DIRECTORY ORGANIZATION

The assumed organization of the directories is the following:
	- One directory to store the surfaces and the resulting skeletons (One directory per case, per side and per subfield)
		surface_directory
		|
		|-- TLE_0362_1_L_CA.obj
		|-- TLE_0362_1_L_CA_inter.obj
		|-- TLE_0362_1_L_CA_skelFinal.obj
		|
		|-- TLE_0362_1_L_SUB.obj
		|-- TLE_0362_1_L_SUB_inter.obj
		|-- TLE_0362_1_L_SUB_skelFinal.obj
		|
		
		
	- One directory to store feature images 
		image_directory
		|
		|-- TLE_0362_1_t2cor-0.4_stx_norm-0.4.mnc
		|-- TLE_0362_1_t2wt1wratio_final.mnc
		|-- TLE_0362_1_Brainclass.mnc
		|-- TLE_0362_1_T2norm.mat
		|-- TLE_0362_1_Ventricle.mnc
		|
		
	
	- One directory to store the txt files resulting from the intersection between the blades and the images
		feature_directory
		|
		|-- TLE_0362_1_L_CA_ColVol.txt
		|-- TLE_0362_1_R_CA_ColVol.txt
		|-- TLE_0362_1_L_CA_nnt2.txt
		|-- TLE_0362_1_R_CA_nnt2.txt
		|-- TLE_0362_1_L_CA_t2wt1wratio.txt
		|-- TLE_0362_1_R_CA_t2wt1wratio.txt
		|

III. BLADE EXTRACTION

The script is used in the following way (we recommend to use the nomenclature stated above):
./get_medial_surface \
	${subfield_surface} \
	${subfield_volume} \
	${individual_output_prefix} # Example: ${SurfaceDirectory}/TLE_0362_1_L_CA

IV. Feature Extraction

This script intersects blade surfaces with images and computes the columnar volume. It is used in the following manner:
./get_individual_features.sh \
	${prefix} \ # e.g. TLE
	${id} \     # e.g. 0362_1
	${subfield} # CA, SUB or DG
	${side}     # L or R
	${surface_directory} \ # The same described in directory organization
	${image_directory} \
	${feature_directory}

V. Lateralization

A serie of script is furnished in order to perform various tasks:
	- A_GetTrainingData.m:  	loads training data and control data set
	- B_StatisticalStudy.m: 	performs a statistical study between controls and training data set
	- C_validation_trainingset.m: 	performs the nested repeated 5-Fold validation based on the training set
	- D_validation_testset.m:	performs the repeated validation on a separated the test set for generalizability
	- E_train_model.m:		performs the training of a final model
	- F_test_individual.m:		

