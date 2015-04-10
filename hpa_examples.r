#------------------------------------------------------------------------------------------------------#
# Purpose	: This code contains somes examples of using 'HPA' modules.
#
# Author	: Jong Cheol Jeong (jjeong1@bidmc.harvard.edu)
#		  Beck Lab @ Beth Israel Deaconess Medical Center, Havard Medical School.
# Date		: 4/2/2015
#------------------------------------------------------------------------------------------------------#
source('hpa.r');

query_str <- 'CD36';			# query string
down_dir  <- './hpa_test_dir';		# downloading directory 

#// downloading images and meta data from TISSUE ATLAS
mtxTissue <- bLab.protein.atlas.tissue(query=query_str, down_home=down_dir);

#// downloading images and meta data from SUBCELL ATLAS
mtxSubcell <- bLab.protein.atlas.subcell(query=query_str, down_home=down_dir);

#// downloading images and meta data from CELL LINE ATLAS
mtxCell <- bLab.protein.atlas.cell_full(query=query_str, down_home=down_dir);

#// downloading images and meta data from CANCER ATLAS
mtxCancer <- bLab.protein.atlas.cancer_full(query=query_str, down_home=down_dir);
