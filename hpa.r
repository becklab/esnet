#------------------------------------------------------------------------------------------------------#
# Purpose	: Downloading images and corresponding meta information.
#
# Usage		: copy this file into your working directory and import this file with following command in R:
#			  > source(hpa.r);
#
# Author	: Jong Cheol Jeong (jjeong1@bidmc.harvard.edu)
#		  Beck Lab @ Beth Israel Deaconess Medical Center, Havard Medical School.
#
# Date		: 4/2/2015
#------------------------------------------------------------------------------------------------------#
######################################################
# Environment setting
######################################################
Sys.setlocale('LC_ALL','C');

######################################################
# Package requirements
######################################################
#// File downloading package
if ("RCurl" %in% rownames(installed.packages()) == FALSE) {
	install.packages('RCurl');
}
require(RCurl);

#// HTML parsing package
if ("XML" %in% rownames(installed.packages()) == FALSE) {
	install.packages('XML');
}
require(XML);

######################################################
# Common functions
######################################################

#********************************************************************
# Remove html tags starts from '<' and ends with '>'
#
# ::Usage::
# 	outDATA <- bLab.data.rm_html(htmlString);
#
# ::Parameters::
#	htmlString(string): string containing HTML tags
#
#********************************************************************
bLab.data.rm_html <- function(htmlString) {
    return(gsub("<.*?>", "", htmlString));
}



#********************************************************************
# Trim leading and tailing spaces
# ::Usage::
# 	paraValue <- bLab.file.trim(x);
#
#********************************************************************
bLab.file.trim <- function( x ) {
  trimX <- gsub("(^[[:space:]]+|[[:space:]]+$)", "", x);
  return(trimX);
}


#********************************************************************
# Create random string based on given length
#
# ::Usage::
# 	outDATA <- bLab.data.random_string(len);
#
# ::Parameters::
#	len(numeric): expected length of random string
#
#********************************************************************
bLab.data.random_string <- function(str_len=1) {
    rndStr <- paste(sample(c(letters, LETTERS, 0:9), size=str_len, replace=TRUE), collapse='');
    return(rndStr);
}


#********************************************************************
# Find index of data containing a given string 
# ::Usage::
# 	infoData <- bLab.data.grep(str, vData, wildcard=TRUE);
#	infoData[[1]]: value containing the string
#	infoData[[2]]: index of matched data
#
# ::Parameters::
#	str 	  (string)  : part of string 
#	vData 	  (string)  : data array to be searched
#	whildcard (boolean) : TRUE=search a string containing 'str', FALSE= seacrh a string exactly match 'str'
#	ignore.case(boolean): TRUE= case insenstive, FASLE= case senstive
#
# :: Outputs::
#	infoData (List): containing both values [[1]] and index [[2]] of 
#   retrieved data
#	e.g.,) bLab.data.grep(str='Sample1', vData=myData);
#
#********************************************************************
bLab.data.grep <- function(str, vData, wildcard=TRUE, ignore.case=TRUE) {
	if (wildcard) {
		ptrn <- sprintf('.*%s.*', str);	
	} else {
		ptrn <- sprintf('^%s$', str);
	}
	
	Val <- grep(ptrn, vData, value=TRUE, ignore.case=ignore.case);		# get the data containing the string
	Idx <- grep(ptrn, vData, value=FALSE, ignore.case=ignore.case);		# get the index of matched data
	return (list(Val, Idx));
}


#********************************************************************
# Download a file until it gets success (i.e., try max_itr times)
#
# ::Usage::
# 	outDATA <- bLab.file.proof_download(dwn_url, fName, max_itr=1000000, verbose=TRUE, dmethod='auto', dmode='wb');
#
# ::Parameters::
#	dwn_url   (string)	: download url
#	fName     (string)	: save downloaded file as fName
#	max_itr	  (integer)	: try 'max_itr' itmes to download file 
#	verbose   (boolean) 	: TRUE= display downloading status, FALSE= do not display downloading status
#	dmethod   (string)	: method: Method to be used for downloading files.  Currently download
#			          methods '"internal"', '"wget"', '"curl"' and '"lynx"'
#			          (deprecated) are available, and there is a value '"auto"':
#			          see 'Details' and 'Note'.
#			          The method can also be set through the option
#			          '"download.file.method"': see 'options()'.
#	dmode     (string)	: character.  The mode with which to write the file.  Useful
#			          values are '"w"', '"wb"' (binary), '"a"' (append) and '"ab"'.
#			          Only used for the '"internal"' method.
#
#********************************************************************
bLab.file.proof_download <- function(dwn_url, fName, max_itr=1000000, verbose=TRUE, dmethod='auto', dmode='wb') {
    cnt <- 0;
    while(1) {
	cnt <- cnt + 1;
	err_flag <- tryCatch(download.file(dwn_url, fName, method=dmethod, quiet=(!verbose), mode=dmode),
			     warning = function(war) {print(war); flush.console(); return(1)},
			     error = function(war) { print(war); flush.console(); return(1)},
			     finally = function(war) {print(war); flush.console(); return(0)}
			     );
	if (!err_flag) {
		break;
	} else if (cnt > max_itr) {
		break;
	}
    }
    return(0);
}



#********************************************************************
# find index of NA and non NA values
# ::Usage::
# 	naInfo <- bLab.data.findNA(vData);
#
# ::Parameters::
#	vData (string)		: array containing 'NA' values
#
# :: Outputs::
#	naInfo[[1]]	: index of 'NA' values
#	naInfo[[2]]	: index of 'non-NA' values
#
#********************************************************************
bLab.data.findNA <- function (vData) {
    NAidx <- which(is.na(vData) == TRUE);
    noNAidx <- which(!(is.na(vData)) == TRUE);
    
    naInfo <- list();
    naInfo[[1]] <- NAidx;
    naInfo[[2]] <- noNAidx;
    
    return(naInfo)
}


#********************************************************************
# Split strings in a matrix and extract targeted item
# Matrix version of strsplit()
#
# ::Usage::
# 	outDATA <- bLab.data.mtxsplit(Mtx, sep, tg_idx)
#
# ::Parameters::
#	Mtx (matrix/data.frame)	: Matrix containing strings
#	sep (character)		: seplit tag (e.g., '-')
#	tg_idx	 (integer)	: index of a targeted item after strings are splitted 
#
# :: Outputs::
#	outDATA		: list
#	outDATA[[1]]	: selected items
#	outDATA[[2]]	: original coordinates of selcted items in Mtx
#			  - 1st column: row index
#			  - 2nd column: column index
#
#********************************************************************
bLab.data.mtxsplit <- function(Mtx, sep, tg_idx) {
    # if array is given then convert it to matrix
    if (is.null(dim(Mtx))) {
	Mtx <- matrix(Mtx, nrow=length(Mtx), ncol=1);
    }
    
    # for ouput
    outDATA <- list();
    selData <- array();
    
    # split the string
    X <- apply(Mtx, 1, strsplit, split=sep);
    mtxIdx <- matrix(data=0, nrow=nrow(Mtx), ncol=2);
    # select data 
    cnt <- 0;
    ridx <- 0;
    for (i in X) {
	ridx <- ridx + 1;
	cidx <- 0;
	for (j in i) {
	    cnt <- cnt + 1;
	    cidx <- cidx + 1;
	    selData[cnt] <- j[tg_idx];
	    mtxIdx[cnt, 1] <- ridx;
	    mtxIdx[cnt, 2] <- cidx;
	}
    }
    outDATA[[1]] <- selData;
    outDATA[[2]] <- mtxIdx;
    return(outDATA);
}


#********************************************************************
# Insert matrix into a existing matrix by given column index
# ::Usage::
# 	newMtx <- bLab.data.InsertCol(tg_idx, inData, orgMtx);
#
# ::Parameters::
#	tg_idx 	(integer)	: row index where the given data (i.e, inData) to be inserted
#	inData	(array/matrix)  : data to be inserted into original matrix (i.e, orgMtx)
#	orgMtx 	(matrix)	: original matrix to be extended 
#
# :: Outputs::
#	newMtx (matrix): extended matrix
#
#
#********************************************************************
bLab.data.InsertCol <- function(tg_idx, inData, orgMtx) {
    #// convert array into matrix
    if (!is.matrix(inData)) {
	inData <- matrix(inData, nrow=length(orgMtx), ncol=1);
    }

    in_nrow <- nrow(inData);
    in_ncol <- ncol(inData);
    org_nrow <- nrow(orgMtx);
    org_ncol <- ncol(orgMtx);
    
    if (in_nrow != org_nrow) {
	cmt <- sprintf('The number of rows between new data (%d) and original data (%d) is different!', in_nrow, org_nrow);
	stop(cmt);
    }
    
    # add inData to the left side of orgMtx
    if ((tg_idx == 1) || (tg_idx == 0)) {
	newMtx <- cbind(inData, orgMtx);
    # add inData to the right side of orgMtx
    } else if ((tg_idx < 0) || (tg_idx > org_ncol)) {
	newMtx <- cbind(orgMtx, inData, deparse.level=0);
    } else {
	leftMtx <- matrix(orgMtx[, 1:tg_idx-1], nrow=org_nrow);
	rightMtx <- matrix(orgMtx[, tg_idx:org_ncol], nrow=org_nrow);
	newMtx <- cbind(leftMtx, inData, rightMtx, deparse.level=0);
    }
    return(newMtx);
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MODUEL GROUP: bLab.protein
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#********************************************************************
# Search protein atlas DB and retrieving primary information excluding internal links
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.data(query=NULL);
#
# ::Parameters::
#	query     (string): search query (e.g., insulin, PGR, CD36)
#
# ::Outputs::
#	outDATA	(matrix): containing all information about the searched results
#	[COLUMN]	[INFORMATION]
#	1		Gene
#	2		Gene description
#	3		Protein class
#	4		URL of tissue information
#	5		URL of subcell information
#	6		URL of cell line information
#	7		URL of cancer information
#	8		RNA tissue category
#
#********************************************************************
bLab.protein.atlas.data <- function(query=NULL) {
    
    #// information about website and search url
    home_url   <- 'http://www.proteinatlas.org'; 
    search_url <- sprintf('%s/search_result.php', home_url);
    
    #// catetory informaion in main data matrix (i.e., initial search results: mtxDATA)
    col_ts 	   <- 4;	# column index of tissue in mtxDATA
    col_sc	   <- 5;	# column index of subcell in mtxDATA
    col_cl	   <- 6;	# column index of cell line in mtxDATA
    col_cc	   <- 7;	# column index of cancer in mtxDATA
    
    #// path for mtxFILE that is temporary datafile to store temporary data
    mtxFILE <- sprintf('./.mtxFILE_%s%s.tmp', format(Sys.time(), '%Y%m%H%M%S'), bLab.data.random_string(6));
    
    #// query must be given
    if (is.null(query)) {
	stop('bLab.protein.atlas.data:: query is missing!');
    }
    
    #// feed parser the retrieved HTML codes
    query    <- paste(strsplit(query, ' ')[[1]], collapse='+');
    htmlMain <- sprintf('%s/search/%s', home_url, query);
    rPage    <- htmlTreeParse(htmlMain, useInternalNodes=TRUE);
    
    #// extract maximum page number
    pages <- xpathSApply(rPage, "//th[@colspan='2' and @class='right']/p", xmlValue);
    pages <- gsub('\\s', '', bLab.data.rm_html(pages));
    #pages <- as.numeric(gsub('.+f([0-9]+).+', '\\1', pages));
    pages <- as.numeric(gsub('.+of(.+)', '\\1', pages));
    
    if (is.na(pages)) {
	pages <- 1;
    }
    
    #// search all pages
    for (p_num in 1: pages) {
	
	htmlPage <- sprintf('%s&page=%d', htmlMain, p_num);
	cmt <- sprintf('[%d/%d] searching %s...\n', p_num, pages, htmlPage);
	cat(cmt); flush.console();
	
	rPage    <- htmlTreeParse(htmlPage, useInternalNodes=TRUE);
	
	#// extract value between <tr><td></td></tr> section
	selTable <- xpathSApply(rPage, "//tbody[@class='hover']/tr/td", xmlValue);
	if (is.null(selTable)) {
	    mtxDATA <- NULL;
	    return(mtxDATA);
	}
	
	#// creating primary data matrix storing overall information of search results
	mtxColNames <- c('gene', 'description', 'protein_class', 'tissue', 'subcell', 'cell_line', 'cancer', 'RNA_tissue_category');
	num_col <- 8;
	num_row <- round(length(selTable)/num_col);
	mtxDATA <- matrix(selTable, nrow=num_row, ncol=num_col, byrow=TRUE);
	colnames(mtxDATA) <- mtxColNames;
	
	#// extracting gene associated URLs (i.e., tissue category)
	cmt <- sprintf('\tbLab.protein.atlas.data:: searching category status...'); cat(cmt); flush.console();
	
	subURLs  <- xpathSApply(rPage, "//tbody[@class='hover']/tr/td[1]//a", xmlGetAttr, 'href');
	
	groups <- c('tissue', 'subcellular', 'cell', 'cancer')
	for (i in 1:length(subURLs)) {
	    # extracting base URL (i.e., remove category information 'tissue' in URL);
	    base_url <- strsplit(subURLs[i], '/')[[1]];
	    base_url <- paste(base_url[1:length(base_url)-1], collapse='/');
	    subURLs[i] <- sprintf('%s%s', home_url, base_url);
	    
	    # check if a category exists: tissue=mtxDATA[,4], subcell=mtxDATA[,5], cell_line=mtxDATA[,6], cancer=mtxDATA[,7]
	    subPage    <- htmlTreeParse(subURLs[i], useInternalNodes=TRUE);
	    selTable <- xpathSApply(subPage, "//div[@class='atlas_header']//ul//li", xmlGetAttr, 'style');
	    # identifying not applicatable category
	    for (j in 1:length(selTable)) {
		if (length(bLab.data.grep('/na.gif', selTable[j])[[1]]) > 0) {
		    mtxDATA[i, col_ts-1+j] <- NA;
		} else {
		    mtxDATA[i, col_ts-1+j] <- sprintf('%s/%s', subURLs[i], groups[j]);
		}
	    }
	}
	
	write.table(mtxDATA, file=mtxFILE, append=TRUE, sep='\t', row.names=FALSE, col.names=FALSE);
	cat('DONE!\n'); flush.console();
    }
    #// writing main information
    
    mtxDATA <- read.table(mtxFILE, sep='\t');
    colnames(mtxDATA) <- mtxColNames;
    file.remove(mtxFILE);
    return(mtxDATA);
}


#********************************************************************
# Parsing tissue category in protein ATLAS 
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.tissue_parse(ts_url);
#
# ::Parameters::
#	ts_url    (string): url for tissue webpage (e.g., http://www.proteinatlas.org/ENSG00000254647-INS/tissue)
#
# ::Outputs::
#	outDATA	(list): containing all information about the searched results
#		outDATA[[1]]: tissue information 
#		outDATA[[2]]: retrieved complete url for images
#		outDATA[[3]]: retrieved short url for images
#		outDATA[[4]]: retrieved image labels
#
#********************************************************************
bLab.protein.atlas.tissue_parse <- function(ts_url, home_url='http://www.proteinatlas.org') {
    #// store data into a matrix 
    tsDATA <- matrix(data=NA, nrow=1, ncol=10);
    
    #// searching tissue URL and get informaion 
    rPage <- htmlTreeParse(ts_url, useInternalNodes=TRUE);
    
    #// varifying missing data & get information
    tissue_cols   <- c('gene', 'description', 'RNA_tissue_category', 'summary', 'expression', 'class', 'localization', 'evidence', 'reliability', 'images');
    tissue_groups <- c("Gene description", "RNA tissue category", "Protein summary", "Protein expression", "Protein class",
		       "Predicted localization", "Protein evidence", "Protein reliability")
    
    tsDATA[1, 1]  <- bLab.file.trim(xpathSApply(rPage, "//div[@class='atlas_header']/div[@class='gene_name']", xmlValue));
    Group 	  <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark staticwidth']/tr//th", xmlValue));
    Info  	  <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark staticwidth']/tr//td", xmlValue));
    for (g in 1:length(Group)) {
	info_idx <- which(tolower(tissue_groups) == tolower(Group[g]));
	tsDATA[1, info_idx+1] <- Info[g]
    }
    
    colnames(tsDATA) <- tissue_cols;
    
    #// Extracting image URLs
    imgDATA <- NULL;
    img_url   <- xpathSApply(rPage, "//table[@style='width: 100%;']//td[@style='background: #C2C2C2; padding-left: 0px; padding-right: 0px; text-align: center;']//a", xmlGetAttr, 'href');
    img_label <- bLab.file.trim(xpathSApply(rPage, "//table[@style='width: 100%;']//td[@style='background: #C2C2C2; padding-left: 0px; padding-right: 0px; text-align: center;']", xmlValue));
    
    if (length(img_url) > 0) {
	for (j in 1:length(img_url)) {
	    imgDATA[j]  <- sprintf('%s%s', home_url, img_url[j]);
	}
    }
    
    return(list(tsDATA, imgDATA, img_url, img_label));
}



#********************************************************************
# Parsing organ section of tissue category in protein ATLAS 
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.tissue_organ_parse(ts_url);
#
# ::Parameters::
#	ts_url    (string): url for tissue webpage (e.g., http://www.proteinatlas.org/ENSG00000254647-INS/tissue)
#
# ::Outputs::
#	outDATA	(matrix):
#		1st column: corresponding gene
#		2nd column: organ name
#		3rd column: stained tissue image URL
#
#********************************************************************
bLab.protein.atlas.tissue_organ_parse <- function(ts_url, home_url='http://www.proteinatlas.org') {
    #// searching tissue URL and get informaion 
    rPage <- htmlTreeParse(ts_url, useInternalNodes=TRUE);
    
    #// varifying missing data & get information
    tissue_cols   <- c('gene', 'organs', 'image_url');
    
    gene <- bLab.file.trim(xpathSApply(rPage, "//div[@class='atlas_header']/div[@class='gene_name']", xmlValue));
    organ <- bLab.file.trim(xpathSApply(rPage, "//tbody[@class='hover']//td[@class='organ']", xmlValue));
    organ_url <- xpathSApply(rPage, "//tbody[@class='hover']//td[@class='organ']//a", xmlGetAttr, 'href');
    
    tsDATA <- matrix(data=NA, nrow=length(organ), ncol=3);
    for (i in 1:nrow(tsDATA)) {
	tsDATA[i, 1] <- gene;
	tsDATA[i, 2] <- organ[i];
	tsDATA[i, 3] <- sprintf('%s%s', home_url, organ_url[i]);
    }
    
    colnames(tsDATA) <- tissue_cols;
    return(tsDATA);
}



#********************************************************************
# Parsing meta information which is hidden inside javascript codes
#
# In html (e.g., http://www.proteinatlas.org/ENSG00000103522-IL21R/tissue/urinary+bladder) 
# "onclick="return imageLoad(1, 11477143, this.href, 2824);" is passed to
# javascript function "function imageLoad(assayId, imgId, imgSrc, geneId, queryString)"
# located in "http://www.proteinatlas.org/search.js?version=13.0.3"
#
# Finally data are retrieved by "function imageMeta(assayId, imgId, queryString)"
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.tissue_meta_parse(meta_url);
#
# ::Parameters::
#	meta_url (string): "imageLoad(1, 11477143, this.href, 2824)" is converted to
#              "http://www.proteinatlas.org/image.php?meta=1&assay_id=1&image_id=11477143"
#
# ::Outputs::
#	outDATA	(matrix):
#    		Column information:
#			[1]'tissue',    [2]'antibody',      [3]'patient_id', [4]'gender',  [5]'age',
#		        [6]'tissue_id', [7]'staining_rate', [8]'intensity',  [9]'quantity',[10]'location'
#
#********************************************************************
bLab.protein.atlas.tissue_meta_parse <- function(meta_url) {
    #// extract meta info
    mtxMeta <- matrix(data=NA, nrow=length(meta_url), ncol=10);
    
    colNames <- c('tissue',    'antibody',      'patient_id', 'gender',  'age',
		  'tissue_id', 'staining_rate', 'intensity',  'quantity', 'location');
    
    for (i in 1:length(meta_url)) {
	tmpPage <- as(htmlTreeParse(meta_url[i], useInternalNodes=TRUE), 'character');
	arrPage <- strsplit(tmpPage, '<br>')[[1]];
	
	# find end of </html>
	htm_idx <- grep('</html>', arrPage, ignore.case=TRUE);
	if (length(htm_idx) > 0) {
	    arrPage[htm_idx] <- "";
	}
	
	# grep information
	mtxMeta[i, 1]  <- bLab.file.trim(gsub('.+<b>(.+)</b>.*', '\\1', arrPage[1]));
	mtxMeta[i, 2]   <- bLab.file.trim(arrPage[2]);

	# patient id
	p_id <- grep('Patient id', arrPage, ignore.case=TRUE)
	if (length(p_id) > 0) {
	    mtxMeta[i, 3] <- bLab.file.trim(bLab.data.rm_html(strsplit(arrPage[p_id], ':')[[1]][2]));
	}
	
	# gender information
	if (length(grep('Female', arrPage)) > 0) {
	    mtxMeta[i, 4] <- 'Female';
	} else if (length(grep('Male', arrPage)) > 0) {
	    mtxMeta[i, 4] <- 'Male';
	}
	
	#// age information
	age_idx <- grep('age', arrPage, ignore.case=TRUE);
	if (length(age_idx) > 0) {
	    age_idx <- age_idx[1];	# age is always located at the first line!
	    mtxMeta[i, 5] <- bLab.file.trim(gsub('.+\\s+([0-9]+).*', '\\1',  arrPage[age_idx]));
	}

	
	#// tissue information
	tss_idx <- grep('\\(.*\\)', arrPage);
	if (length(tss_idx) > 0) {
	    mtxMeta[i, 6] <- bLab.file.trim(bLab.data.rm_html(paste(arrPage[tss_idx], collapse=';')));
	}

	#// antibody staining information
	as_idx <- grep('Staining:', arrPage, ignore.case=TRUE)
	if (length(as_idx) > 0) {
	    mtxMeta[i, 7] <- bLab.file.trim(bLab.data.rm_html(strsplit(arrPage[as_idx], ':')[[1]][2]));
	}

	#// Intensity information
	it_idx <- grep('Intensity:', arrPage, ignore.case=TRUE)
	if (length(it_idx) > 0) {
	    mtxMeta[i, 8] <- bLab.file.trim(bLab.data.rm_html(strsplit(arrPage[it_idx], ':')[[1]][2]));
	}

	#// Quantity information
	qt_idx <- grep('Quantity:', arrPage, ignore.case=TRUE)
	if (length(qt_idx) > 0) {
	    qtt <- bLab.file.trim(bLab.data.rm_html(strsplit(arrPage[qt_idx], ':')[[1]][2]));
	    mtxMeta[i, 9] <- gsub('&gt;', '>', qtt);
	    
	}

	#// Location information
	lc_idx <- grep('Location:', arrPage, ignore.case=TRUE);
	if (length(lc_idx) >0) {
	    block_lc_idx <- grep('</span>', arrPage, ignore.case=TRUE);
	    if (length(block_lc_idx) > 0) {
		tmpInfo <- paste(arrPage[lc_idx[1]:max(block_lc_idx)], collapse="");	
	    } else {
		tmpInfo <- paste(arrPage[lc_idx[1]:length(arrPage)], collapse="");
	    }
	}
	
	if (length(lc_idx) > 0) {
	    tmp_lc <- bLab.file.trim(bLab.data.rm_html(strsplit(tmpInfo, ':')[[1]][2]));
	    mtxMeta[i, 10] <- bLab.file.trim(gsub("(.*)'.*", '\\1', tmp_lc));
	}
    }
    colnames(mtxMeta) <- colNames;
    return(mtxMeta);
}


#********************************************************************
# Parsing stained tissue image in tissue category at protein ATLAS 
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.tissue_stain_image_parse(ts_url);
#
# ::Parameters::
#	ts_url    (string): url for tissue webpage (e.g., http://www.proteinatlas.org/ENSG00000103522-IL21R/tissue/urinary+bladder)
#			    this can be obtained from "bLab.protein.atlas.tissue_organ_parse"
#
# ::Outputs::
#	outDATA	(matrix): containing all information about images
#	column information
#    		[1]'gene',      [2]'tissue',        [3]'antibody',  [4]'patient_id', [5]'gender',   [6]'age',
#		[7]'tissue_id', [8]'staining_rate', [9]'intensity', [10]'quantity',  [11]'location',[12]'img_url'
#
#********************************************************************
bLab.protein.atlas.tissue_stain_image_parse <- function(ts_url, home_url='http://www.proteinatlas.org') {
    #// parsing url 
    rPage <- htmlTreeParse(ts_url, useInternalNodes=TRUE);
    
    gname  <- bLab.file.trim(xpathSApply(rPage, "//div[@class='atlas_header']/div[@class='gene_name']", xmlValue));
    antiB  <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark']/tr//th[@class='head']//p", xmlValue));

    image_url <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark']//td[@class='nopadd']/a", xmlGetAttr, 'href'));
    meta_info <-bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark']//td[@class='nopadd']/a", xmlGetAttr, 'onclick'));
    
    if ((length(image_url) < 1) || (length(meta_info) < 1)) {
	return(NA);
    }
    #// get meta informaton which is located in function imageMeta(assayId, imgId, queryString) in http://www.proteinatlas.org/image.js?version=13.0.3
    meta_url  <- array();
    for (i in 1:length(meta_info)) {
        tmp_info <- strsplit(meta_info[i], ',')[[1]];
	if (length(tmp_info) >= 4) {
	    assay_id <- gsub('return.+[\\(\\)](.+)', '\\1', tmp_info[1]);
	    image_id <- bLab.file.trim(tmp_info[2]);
	    dummy1   <- bLab.file.trim(tmp_info[3]);
	    dummy2   <- bLab.file.trim(gsub('(.+)[\\(\\)].+', '\\1', tmp_info[4]));
	    meta_url[i] <- sprintf('%s/image.php?meta=1&assay_id=%s&image_id=%s', home_url, assay_id, image_id);
	} else {
	    meta_url[i] <- NA;
	}
    }
    
    #// get meta information
    mtxMeta <- bLab.protein.atlas.tissue_meta_parse(meta_url);

    
    #// gathering image information
    num_col  <- ncol(mtxMeta) + 2;
    imgInfo  <- matrix(data=NA, nrow=length(image_url), ncol=num_col);
    
    imgInfo[, 2:(num_col-1)] <- mtxMeta;				# put mtxMeta between first and last column of new matrix (i.e., imgInfo)
    
    for (i in 1:nrow(imgInfo)) {
	imgInfo[i, 1]       <- gname;
	imgInfo[i, num_col] <- sprintf('%s%s', home_url, image_url[i]);
    }
    
    colNames <- c('gene',      'tissue',        'antibody',  'patient_id', 'gender',   'age',
		  'tissue_id', 'staining_rate', 'intensity', 'quantity',   'location', 'img_url');
    
    colnames(imgInfo) <- colNames;
    return (imgInfo);
}

#********************************************************************
# Search protein atlas DB and download tissue images
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.tissue(query=NULL, down_home='./');
#
# ::Parameters::
#	query     (string): search query (e.g., insulin, PGR, CD36)
#	down_home (string): home directory where files can be downloaded
#		(e.g., query='CD36', down_home='/my/download' will create directory with 'query'
#                so, file will be stored in '/my/download/CD36/...');
#	verbose  (boolean): TRUE= display downloading status, FALSE= do not show downloading status
#
# ::Outputs::
#	outDATA	(list): containing all information about the searched results
#		outDATA[[1]]: main information about search results
#		outDATA[[2]]: details about tissue and corresponding genes
#		outDATA[[3]]: downloaded tissue image information
#		outDATA[[4]]: downlaoded high resolution tissue image information
#
#	TEXT files: downloaded information in text file
#		/down_home/<query>/main_tissue_<query>.txt : tab delimited text file of outDATA[[1]]
#		/down_home/<query>/image_tissue_<query>.txt: tab delimited text file of downloaded images
#		/down_home/<query>/info_tissue_<query>.txt : tab delimited text file of retrieved tissue information
#		/down_home/<query>/stain_tissue_<query>.txt : tab delimited text file of overall tissue data of
#							      (e.g., http://www.proteinatlas.org/ENSG00000103522-IL21R/tissue/urinary+bladder)
#							      and meta data. For details please refer function 'bLab.protein.atlas.tissue_stain_image_parse'
#		/down_home/<query>/stain_info_tissue_<query>.txt : tab delimited text file of organ system information
#								   (e.g., http://www.proteinatlas.org/ENSG00000103522-IL21R/tissue)
#		/down_home/<query>/stain_img_tissue_<query>.txt : tab delimited text file of downloaded images 
#
#********************************************************************
bLab.protein.atlas.tissue <- function (query=NULL, down_home='./', verbose=TRUE) {
    #------------------------------------------------------------------------------
    # configuration for downloading environment
    #------------------------------------------------------------------------------
    # defining path for downloaded images
    img_type <- 'tissue';
    
    #down_path <- sprintf('%s/%s', down_home, query);
    str_query <- paste(strsplit(query, ' ')[[1]], collapse='_');
    down_path <- sprintf('%s/%s', down_home, str_query);
    
    dir.create(down_path, showWarnings = FALSE, recursive = TRUE);

    down_tissue <- sprintf('%s/%s', down_path, img_type);
    dir.create(down_tissue, showWarnings = FALSE, recursive = TRUE);
    
    # log file to record downloaded images
    img_log <- sprintf('%s/image_%s_%s.txt', down_path, img_type, str_query);
    fid_img <- file(img_log, 'w');
    cmt <- 'query\tgene\tcategory\tloci\timage_url\timage_file';	# image log file header
    writeLines(cmt, fid_img);
    
    # other log files
    main_log <- sprintf('%s/main_%s_%s.txt', down_path, img_type, str_query);		# overall information about searched results
    ts_log   <- sprintf('%s/info_%s_%s.txt', down_path, img_type, str_query);		# tissue information

    stain_log <- sprintf('%s/stain_%s_%s.txt', down_path, img_type, str_query);		# stain information
    if (file.exists(stain_log)) {
	file.remove(stain_log);
    }
    
    stain_info_log <- sprintf('%s/stain_info_%s_%s.txt', down_path, img_type, str_query);	# stain information
    if (file.exists(stain_info_log)) {
	file.remove(stain_info_log);
    }

    stain_img_log  <- sprintf('%s/stain_img_%s_%s.txt', down_path, img_type, str_query);	# stain information
    if (file.exists(stain_img_log)) {
	file.remove(stain_img_log);
    }
    
    #// get retrieved information
    mtxDATA <- bLab.protein.atlas.data(query=query);
    if (is.null(mtxDATA)) {
	return(list(NULL, NULL, NULL));
    }
    num_row <- nrow(mtxDATA);
    col_ts  <- 4;		# column containing URL of tissue information
    
    #// writing main information
    write.table(mtxDATA, main_log, row.names=FALSE);
    cat('DONE!\n'); flush.console();
    
    #------------------------------------------------------------------------------
    # Downloading TISSUE data
    #------------------------------------------------------------------------------
    #// retrieving data from each category
    tsDATA <- matrix(data=NA, nrow=num_row, ncol=10);	# tissue data
    tsDATA[, 1] <- mtxDATA[, 1]				# tsDATA and mtxDATA has same number of rows
    
    for (i in 1:num_row) {
	# if information exists for tissue then retrieving information
	ts_url <- mtxDATA[i, col_ts];
	if (!is.na(ts_url)) {
	    tmpDATA 	<- bLab.protein.atlas.tissue_parse(ts_url);
	    tsDATA[i, ] <- tmpDATA[[1]];	# update tsDATA based on individual tissue URLs
	    if (!is.null(tmpDATA[[2]])) {
		img_url   <- tmpDATA[[2]];	# image full url
		img_label <- tmpDATA[[4]];	# image label 
		for (j in 1:length(img_url)) {
		    tmp_img_url  <- img_url[j];
		    tmp_name <- strsplit(img_url[j], '/')[[1]];
		    tmp_name <- tmp_name[length(tmp_name)];
		    
		    # file image name = tissue_path + image_label() + original name 
		    tmp_img_name <- sprintf('%s/%s_%s_%s',down_tissue, tsDATA[i, 1], gsub(" ", "_", img_label[j]), tmp_name);
		    cmt <- sprintf('%s\t%s\t%s\t%s\t%s\t%s', query,  tsDATA[i, 1], img_type, img_label[j], tmp_img_url, tmp_img_name);
		    writeLines(cmt, fid_img);
		    #download.file(tmp_img_url, tmp_img_name, method='auto', quiet=(!verbose));
		    bLab.file.proof_download(tmp_img_url, tmp_img_name, max_itr=1000000, verbose=TRUE);
		}
		tsDATA[i, 10] <- img_log;
	    }
	    
	    #// writing stained image information
	    tmpDATA2    <- bLab.protein.atlas.tissue_organ_parse(ts_url);	# get url for stain images
	    write.table(tmpDATA2, file=stain_info_log, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t');
	    
	    for (i_url in 1:nrow(tmpDATA2)) {
		tmp_url <- tmpDATA2[i_url, 3];
		if (!is.na(tmp_url)) {
		    tmp_imgInfo <- bLab.protein.atlas.tissue_stain_image_parse(tmp_url);
		    if (!is.na(tmp_imgInfo[1])) {
			write.table(tmp_imgInfo, file=stain_log, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t');	
		    }
	    	    
		}
	    }
	}
    }
    tissue_cols   <- c('gene', 'description', 'RNA_tissue_category', 'summary', 'expression', 'class', 'localization', 'evidence', 'reliability', 'images');
    colnames(tsDATA) <- tissue_cols;
    write.table(tsDATA, ts_log, row.names=FALSE);
    cat('DONE!\n'); flush.console();
    close(fid_img);
    imgDATA <- read.table(img_log, header = TRUE, sep = "\t",);
    
    #// downloading stain images
    imgS      <- as.matrix(read.table(stain_log, sep = "\t",));			# read all stain images
    ncol_imgS <- ncol(imgS);
    nrow_imgS <- nrow(imgS);
    
    imgDATA2  <- matrix(data=NA, nrow=nrow(imgS), ncol=(ncol_imgS + 1));
    imgDATA2[ ,1:ncol_imgS] <- as.matrix(imgS);		# adding one more column (i.e., path to store images) into imgH
    
    for (s in 1:nrow_imgS) {
	tmp_img_url  <- imgS[s, ncol_imgS];
	tmp_img_name <- strsplit(tmp_img_url, '/')[[1]];
	img_label <- paste(imgS[s, c(1,2,3,5)], collapse='_');
	img_label <- gsub(" ", "_", img_label);
	tmp_img_name <- sprintf('%s/%s_%s', down_tissue, img_label, tmp_img_name[length(tmp_img_name)]);
	bLab.file.proof_download(tmp_img_url, tmp_img_name, max_itr=1000000, verbose=TRUE);
	imgDATA2[s, ncol_imgS+1] <- tmp_img_name;
    }
    colNames <- c('gene',     'tissue',    'antibody',      'patient_id', 'gender',
		  'age',      'tissue_id', 'staining_rate', 'intensity',  'quantity',
		  'location', 'image_url', 'image_path');
    colnames(imgDATA2) <- colNames;
    write.table(imgDATA2, file=stain_img_log, append=FALSE, row.names=FALSE, col.names=TRUE, sep='\t');
    
    #// return all values
    return(list(mtxDATA, tsDATA, imgDATA, imgDATA2));
}



#********************************************************************
# Parsing subcellular category in protein ATLAS 
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.subcell_parse(sc_url);
#
# ::Parameters::
#	sc_url    (string): url for tissue webpage (e.g., http://www.proteinatlas.org/ENSG00000171105-INSR/subcellular)
#
# ::Outputs::
#	outDATA	(list): containing all information about the searched results
#		outDATA[[1]]: tissue information 
#		outDATA[[2]]: retrieved complete url for images
#		outDATA[[3]]: retrieved short url for images
#		outDATA[[4]]: retrieved image labels
#
#********************************************************************
bLab.protein.atlas.subcell_parse <- function(sc_url, home_url='http://www.proteinatlas.org') {
    #// store data into a matrix 
    sc_num_col <- 8;
    scDATA <- matrix(data=NA, nrow=1, ncol=sc_num_col);	# subcell data
    
    #// searching subcellular URL and get informaion 
    rPage <- htmlTreeParse(sc_url, useInternalNodes=TRUE);
    
    #// varifying missing data & get information
    subcell_cols   <- c('gene', 'summary', 'main_location', 'add_location', 'reliability', 'evidence', 'assay_summary', 'images');
    subcell_groups <- c("Summary", "Main location", "Additional location", "Reliability", "Protein evidence", "Assay summary")
    
    scDATA[1, 1]  <- bLab.file.trim(xpathSApply(rPage, "//div[@class='atlas_header']/div[@class='gene_name']", xmlValue));
    Group 	  <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark staticwidth']/tr//th", xmlValue));
    Info  	  <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark staticwidth']/tr//td", xmlValue));
    for (g in 1:length(Group)) {
	info_idx <- which(tolower(subcell_groups) == tolower(Group[g]));
	scDATA[1, info_idx+1] <- Info[g]
    }
    
    #// Extracting image URLs
    imgDATA   <- NULL;
    img_url   <- xpathSApply(rPage, "//table[@class='dark']//table//a[@class='draggable']", xmlGetAttr, 'href');
    img_label <- xpathSApply(rPage, "//table[@class='dark']//table//img", xmlGetAttr, 'alt');

    if (length(img_url) > 0) {
	for (j in 1:length(img_url)) {
	    imgDATA[j]  <- sprintf('%s%s', home_url, img_url[j]);
	}
    }
    
    return(list(scDATA, imgDATA, img_url, img_label));
}


#********************************************************************
# Search protein atlas DB and download subcellular images
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.subcell(query=NULL, down_home='./');
#
# ::Parameters::
#	query     (string): search query (e.g., insulin, PGR, CD36)
#	down_home (string): home directory where files can be downloaded
#		(e.g., query='CD36', down_home='/my/download' will create directory with 'query'
#                so, file will be stored in '/my/download/CD36/...');
#	verbose  (boolean): TRUE= display downloading status, FALSE= do not show downloading status
#
# ::Outputs::
#	outDATA	(list): containing all information about the searched results
#		outDATA[[1]]: main information about search results
#		outDATA[[2]]: details about tissue and corresponding genes
#		outDATA[[3]]: downloaded tissue image information
#
#	TEXT files: downloaded information in text file
#		/down_home/<query>/main_subcellular_<query>.txt : tab delimited text file of outDATA[[1]]
#		/down_home/<query>/image_subcellular_<query>.txt: tab delimited text file of downloaded images
#		/down_home/<query>/info_subcellular_<query>.txt : tab delimited text file of retrieved tissue information
#
#********************************************************************
bLab.protein.atlas.subcell <- function (query=NULL, down_home='./', verbose=TRUE) {
    
    #------------------------------------------------------------------------------
    # configuration for downloading environment
    #------------------------------------------------------------------------------
    # defining path for downloaded images
    img_type <- 'subcellular';
    
    #down_path <- sprintf('%s/%s', down_home, query);
    str_query <- paste(strsplit(query, ' ')[[1]], collapse='_');
    down_path <- sprintf('%s/%s', down_home, str_query);
    dir.create(down_path, showWarnings = FALSE, recursive = TRUE);

    down_subcell <- sprintf('%s/%s', down_path, img_type);
    dir.create(down_subcell, showWarnings = FALSE, recursive = TRUE);
    
    # log file to record downloaded images
    img_log <- sprintf('%s/image_%s_%s.txt', down_path, img_type, str_query);
    fid_img <- file(img_log, 'w');
    cmt <- 'query\tgene\tcategory\tID\timage_url\timage_file';	# image log file header
    writeLines(cmt, fid_img);
    
    # other log files
    main_log <- sprintf('%s/main_%s_%s.txt', down_path, img_type, str_query);	# overall information about searched results
    sc_log   <- sprintf('%s/info_%s_%s.txt', down_path, img_type, str_query);	# subcellular information
    
    #// get retrieved information
    mtxDATA <- bLab.protein.atlas.data(query=query);
        if (is.null(mtxDATA)) {
	return(list(NULL, NULL, NULL));
    }

    num_row <- nrow(mtxDATA);
    col_sc  <- 5;		# column containing URL of subcellular information
    
    #// writing main information
    write.table(mtxDATA, main_log, row.names=FALSE);
    cat('DONE!\n'); flush.console();
    
    #------------------------------------------------------------------------------
    # Downloading Subcellular data
    #------------------------------------------------------------------------------
    #// retrieving data from each category
    scDATA <- matrix(data=NA, nrow=num_row, ncol=8);	# subcellular data
    scDATA[, 1] <- mtxDATA[, 1]				# scDATA and mtxDATA has same number of rows
    
    for (i in 1:num_row) {
	# if information exists for subcellular then retrieving information
	sc_url <- mtxDATA[i, col_sc];
	if (!is.na(sc_url)) {
	    tmpDATA 	<- bLab.protein.atlas.subcell_parse(sc_url);
	    scDATA[i, ] <- tmpDATA[[1]];	# update scDATA based on individual tissue URLs
	    if (!is.null(tmpDATA[[2]])) {
		img_url   <- tmpDATA[[2]];	# image full url
		img_label <- tmpDATA[[4]];	# image label 
		for (j in 1:length(img_url)) {
		    tmp_img_url  <- img_url[j];
		    tmp_name <- strsplit(img_url[j], '/')[[1]];
		    tmp_name <- tmp_name[length(tmp_name)];
		    
		    # file image name = subcell_path + image_label() + original name
		    tmp_label_name  <- gsub('[/,: ]', '-', img_label[j]);
		    tmp_label_name1 <- gsub('[/,: ]', '-', tmp_label_name);
		    tmp_label_name  <- gsub('[/,: ]', '-', tmp_label_name1);
		    
		    tmp_img_name <- sprintf('%s/%s_%s_%s',down_subcell, scDATA[i, 1], tmp_label_name, tmp_name);
		    cmt <- sprintf('%s\t%s\t%s\t%s\t%s\t%s', query,  scDATA[i, 1], img_type, tmp_label_name1[1], tmp_img_url, tmp_img_name);
		    writeLines(cmt, fid_img);
		    #download.file(tmp_img_url, tmp_img_name, method='auto', quiet=(!verbose));
		    bLab.file.proof_download(tmp_img_url, tmp_img_name, max_itr=1000000, verbose=TRUE);
		}
		scDATA[i, 8] <- img_log;
	    }
	}
    }
    subcell_cols <- c('gene', 'summary', 'main_location', 'add_location', 'reliability', 'evidence', 'assay_summary', 'images');
    colnames(scDATA) <- subcell_cols;
    write.table(scDATA, sc_log, row.names=FALSE);
    cat('DONE!\n'); flush.console();
    close(fid_img);
    imgDATA <- read.table(img_log, header = TRUE, sep = "\t",);
    
    #// return all values
    return(list(mtxDATA, scDATA, imgDATA));
}




#********************************************************************
# Parsing cell line category in protein ATLAS 
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.cell_parse(cl_url);
#
# ::Parameters::
#	cl_url    (string): url for cell line webpage (e.g., http://www.proteinatlas.org/ENSG00000146247-PHIP/cell/HPA019838)
#
# ::Outputs::
#	outDATA	(list): containing all information about the searched results
#		outDATA[[1]]: cell line information 
#		outDATA[[2]]: retrieved complete url for images
#		outDATA[[3]]: retrieved short url for images
#		outDATA[[4]]: retrieved image labels
#
#********************************************************************
bLab.protein.atlas.cell_parse <- function(cl_url, home_url='http://www.proteinatlas.org') {
    #// store data into a matrix 
    cl_num_col <- 7;
    clDATA <- matrix(data=NA, nrow=1, ncol=cl_num_col);	# cell line data
    
    #// searching cell line URL and get informaion 
    rPage <- htmlTreeParse(cl_url, useInternalNodes=TRUE);
    
    #// varifying missing data & get information
    cell_cols   <- c('gene', 'description', 'rna_expression', 'protein_expression', 'protein_class', 'reliability', 'images');
    cell_groups <- c("Gene description", "RNA expression", "Protein expression", "Protein class", "Reliability");
    
    clDATA[1, 1]  <- bLab.file.trim(xpathSApply(rPage, "//div[@class='atlas_header']/div[@class='gene_name']", xmlValue));
    Group 	  <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark staticwidth']/tr//th", xmlValue));
    Info  	  <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark staticwidth']/tr//td", xmlValue));
    for (g in 1:length(Group)) {
	info_idx <- which(tolower(cell_groups) == tolower(Group[g]));
	clDATA[1, info_idx+1] <- Info[g]
    }
    
    #// Extracting image URLs
    imgDATA   <- NULL;
    img_url   <- xpathSApply(rPage, "//table[@class='dark']//th[@style='width: 150px; border-top: 0; padding: 0px']/a", xmlGetAttr, 'href');
    img_label <- xpathSApply(rPage, "//table[@class='dark']//th[@style='width: 150px; border-top: 0; padding: 0px']//p", xmlValue);

    if (length(img_url) > 0) {
	for (j in 1:length(img_url)) {
	    imgDATA[j]  <- sprintf('%s%s', home_url, img_url[j]);
	}
    }
    
    return(list(clDATA, imgDATA, img_url, img_label));
}




#********************************************************************
# Parsing cell line category in protein ATLAS 
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.cell_summary_parse(cl_url);
#
# ::Parameters::
#	cl_url    (string): url for cell line webpage (e.g., http://www.proteinatlas.org/ENSG00000146247-PHIP/cell/HPA019838)
#
# ::Outputs::
#	outDATA	(matrix): ['gene', 'cell line', 'cell line category', 'antibody', 'image url']; 
#
#********************************************************************
bLab.protein.atlas.cell_summary_parse <- function(cl_url) {
    #// extract home_url
    home_url <- strsplit(cl_url, '/')[[1]]
    home_url <- paste(home_url[1:(length(home_url)-2)], collapse='/');
    
    #// extracting base url
    cl_home <- strsplit(cl_url, '/')[[1]]
    cl_home <- paste(cl_home[1:(length(cl_home)-1)], collapse='/');
    
    #// searching tissue URL and get informaion 
    rPage <- htmlTreeParse(cl_url, useInternalNodes=TRUE);
    
    #// varifying missing data & get information
    cell_cols   <- c('gene', 'cell_line', 'cLines', 'antibody', 'image_url');
    
    gene       <- bLab.file.trim(xpathSApply(rPage, "//div[@class='atlas_header']/div[@class='gene_name']", xmlValue));
    cl_img_url <- xpathSApply(rPage, "//tbody[@class='hover']//td[@class='organ']//a", xmlGetAttr, 'href');
    
    #// extracting antibody
    antiB  <- bLab.file.trim(xpathSApply(rPage, "//div[@class='menu']//span/a", xmlGetAttr, 'href'));
    a_idx  <- grep('/cell', antiB);
    antiB  <- antiB[a_idx];
    
    cell_line  <- array();
    cLines     <- array();
    for (i in 1: length(cl_img_url)) {
	synX <- sprintf("//tbody[@class='hover']/tr//td[@class='organ' and @style='text-align: center;']/a[@href='%s']", cl_img_url[i]);
	cell_line[i]  <- bLab.file.trim(xpathSApply(rPage, synX, xmlValue));
	synX <- sprintf("//tbody[@class='hover']/tr//td[@class='organ' and @style='text-align: left; padding-left:20px; padding-right: 2px;']/a[@href='%s']", cl_img_url[i]);
	cLines[i]     <- bLab.file.trim(xpathSApply(rPage, synX, xmlValue));
    }
    
    #// remove duplicated items
    uqIdx <- which(!duplicated(cl_img_url) == TRUE);
    cl_img_url <- cl_img_url[uqIdx];
    cell_line  <- cell_line[uqIdx];
    cLines     <- cLines[uqIdx];
    cl_search  <- bLab.data.mtxsplit(cl_img_url, '/', 2)[[1]]
    
    
    clDATA <- matrix(data=NA, nrow=length(cl_search) * length(antiB), ncol=5);
    cnt <- 0;
    for (i in 1:length(cl_search)) {
	for (j in 1:length(antiB)) {
	    cnt <- cnt + 1;
	    clDATA[cnt, 1] <- gene;
	    clDATA[cnt, 2] <- cell_line[i];
	    clDATA[cnt, 3] <- cLines[i];
	    clDATA[cnt, 4] <- strsplit(antiB[j], '/')[[1]][4];
	    clDATA[cnt, 5] <- sprintf('%s%s/%s', home_url, antiB[j], cl_search[i]);

	}
    }
    
    colnames(clDATA) <- cell_cols;
    return(clDATA);
}


#********************************************************************
# Parsing meta information which is hidden inside javascript codes
#
# In html (e.g., http://www.proteinatlas.org/ENSG00000138760-SCARB2/cell/HPA018014/HL-60) 
# "onclick="return imageLoad(3, 5381314, this.href, 7686);" is passed to
# javascript function "function imageLoad(assayId, imgId, imgSrc, geneId, queryString)"
# located in "http://www.proteinatlas.org/search.js?version=13.0.3"
#
# Finally data are retrieved by "function imageMeta(assayId, imgId, queryString)"
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.cell_meta_parse(meta_url);
#
# ::Parameters::
#	meta_url (string): "imageLoad((3, 5381314, this.href, 7686)" is converted to
#              "http://www.proteinatlas.org/image.php?meta=1&assay_id=3&image_id=5381314"
#
# ::Outputs::
#	outDATA	(matrix):
#    		Column information:
#			[1]'cell_line_name', [2]'gender', [3]'age', [4]'cell_line_id', [5]'number_of_cell', [6]'positive_cells'
#
#********************************************************************
bLab.protein.atlas.cell_meta_parse <- function(meta_url) {
    #// extract meta info
    mtxMeta <- matrix(data=NA, nrow=length(meta_url), ncol=6);
    
    colNames <- c('cell_line_name', 'gender', 'age', 'cell_line_id', 'number_of_cell', 'positive_cells');
    
    for (i in 1:length(meta_url)) {
	tmpPage <- as(htmlTreeParse(meta_url[i], useInternalNodes=TRUE), 'character');
	arrPage <- strsplit(tmpPage, '<br>')[[1]];
	
	# find end of </html>
	htm_idx <- grep('</html>', arrPage, ignore.case=TRUE);
	if (length(htm_idx) > 0) {
	    arrPage[htm_idx] <- "";
	}
	
	# grep information
	mtxMeta[i, 1]  <- bLab.file.trim(gsub('.+<b>(.+)</b>.*', '\\1', arrPage[1]));

	# gender information
	if (length(grep('Female', arrPage)) > 0) {
	    mtxMeta[i, 2] <- 'Female';
	} else if (length(grep('Male', arrPage)) > 0) {
	    mtxMeta[i, 2] <- 'Male';
	}
	
	#// age information
	age_idx <- grep('age', arrPage, ignore.case=TRUE);
	if (length(age_idx) > 0) {
	    age_idx <- age_idx[1];	# age is always located at the first line!
	    mtxMeta[i, 3] <- bLab.file.trim(gsub('.+\\s+([0-9]+).*', '\\1',  arrPage[age_idx]));
	}
	
	# cell line number
	cell_id <- grep('cell line', arrPage, ignore.case=TRUE);
	if (length(cell_id) > 0) {
	    mtxMeta[i, 4] <- bLab.file.trim(bLab.data.rm_html(strsplit(arrPage[cell_id], ' ')[[1]][3]));
	}
	
	# number of cells
	cell_num <- grep('cells', arrPage, ignore.case=TRUE);
	if (length(cell_num) > 0) {
	    mtxMeta[i, 5] <- bLab.file.trim(bLab.data.rm_html(strsplit(arrPage[cell_num], ' ')[[1]][1]));
	}

	# positive rate
	pos_id<- grep('\\%', arrPage, ignore.case=TRUE);
	if (length(pos_id) > 0) {
	    mtxMeta[i, 6] <- bLab.file.trim(bLab.data.rm_html(strsplit(arrPage[pos_id], ' ')[[1]][1]));
	}
    }
    colnames(mtxMeta) <- colNames;
    return(mtxMeta);
}



#********************************************************************
# Parsing stained cell line image in cell line category at protein ATLAS 
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.cell_stain_parse(cl_url);
#
# ::Parameters::
#	cl_url    (string): url for tissue webpage (e.g., http://www.proteinatlas.org/ENSG00000138760-SCARB2/cell/CAB015415/U-266%2F70)
#			    this can be obtained from "bLab.protein.atlas.tissue_organ_parse"
#
# ::Outputs::
#	outDATA	(matrix): containing all information about images
#	column information
#    		[1]'gene',           [2]'cell_line_name', [3]'gender', [4]'age', [5]'cell_line_id',
#		[6]'number_of_cell', [7]'positive_cells', [8]'img_url'
#
#********************************************************************
bLab.protein.atlas.cell_stain_parse <- function(cl_url, home_url='http://www.proteinatlas.org') {
    #// parsing url 
    rPage <- htmlTreeParse(cl_url, useInternalNodes=TRUE);
    
    gname  <- bLab.file.trim(xpathSApply(rPage, "//div[@class='atlas_header']/div[@class='gene_name']", xmlValue));
    antiB  <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark']/tr//th[@class='head']//p", xmlValue));
    ab_idx <- grep('Antibody ', antiB);
    antiB  <- antiB[ab_idx];

    image_url <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark']//td[@class='nopadd']/a", xmlGetAttr, 'href'));
    meta_info <-bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark']//td[@class='nopadd']/a", xmlGetAttr, 'onclick'));
	
    if ((length(image_url) < 1) || (length(meta_info) < 1)) {
	return(NA);
    }
    
    #// get meta informaton which is located in function imageMeta(assayId, imgId, queryString) in http://www.proteinatlas.org/image.js?version=13.0.3
    meta_url  <- array();
    for (i in 1:length(meta_info)) {
	tmp_info <- strsplit(meta_info[i], ',')[[1]];
	if (length(tmp_info) >= 4) {
	    assay_id <- gsub('return.+[\\(\\)](.+)', '\\1', tmp_info[1]);
	    image_id <- bLab.file.trim(tmp_info[2]);
	    dummy1   <- bLab.file.trim(tmp_info[3]);
	    dummy2   <- bLab.file.trim(gsub('(.+)[\\(\\)].+', '\\1', tmp_info[4]));
	    meta_url[i] <- sprintf('%s/image.php?meta=1&assay_id=%s&image_id=%s', home_url, assay_id, image_id);
	} else {
	    meta_url[i] <- NA;
	}
    }
    
    #// get meta information
    mtxMeta <- bLab.protein.atlas.cell_meta_parse(meta_url);
    
    #// gathering image information
    num_col  <- ncol(mtxMeta) + 2;
    imgInfo  <- matrix(data=NA, nrow=length(image_url), ncol=num_col);
    imgInfo[, 2:(num_col-1)] <- mtxMeta;				# put mtxMeta between first and last column of new matrix (i.e., imgInfo)
    for (i in 1:nrow(imgInfo)) {
	imgInfo[i, 1]       <- gname;
	imgInfo[i, num_col] <- sprintf('%s%s', home_url, image_url[i]);
    }
    
    
    colNames <- c('gene', 'cell_line_name', 'gender', 'age', 'cell_line_id', 'number_of_cell', 'positive_cells', 'img_url');
    
    colnames(imgInfo) <- colNames;
    return (imgInfo);
}


#********************************************************************
# Search protein atlas DB and download cell line images
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.cell(query=NULL, down_home='./');
#
# ::Parameters::
#	query     (string): search query (e.g., insulin, PGR, CD36, SCARB2)
#	down_home (string): home directory where files can be downloaded
#		(e.g., query='CD36', down_home='/my/download' will create directory with 'query'
#                so, file will be stored in '/my/download/CD36/...');
#	verbose  (boolean): TRUE= display downloading status, FALSE= do not show downloading status
#
# ::Outputs::
#	outDATA	(list): containing all information about the searched results
#		outDATA[[1]]: main information about search results
#		outDATA[[2]]: details about cell line and corresponding genes
#		outDATA[[3]]: downloaded cell line image information
#
#	TEXT files: downloaded information in text file
#		/down_home/<query>/main_cell_<query>.txt : tab delimited text file of outDATA[[1]]
#		/down_home/<query>/image_cell_<query>.txt: tab delimited text file of downloaded images
#		/down_home/<query>/info_cell_<query>.txt : tab delimited text file of retrieved tissue information
#
#********************************************************************
bLab.protein.atlas.cell <- function (query=NULL, down_home='./', verbose=TRUE) {
    
    #------------------------------------------------------------------------------
    # configuration for downloading environment
    #------------------------------------------------------------------------------
    # defining path for downloaded images
    img_type <- 'cell';
    
    #down_path <- sprintf('%s/%s', down_home, query);
    str_query <- paste(strsplit(query, ' ')[[1]], collapse='_');
    down_path <- sprintf('%s/%s', down_home, str_query);
    
    dir.create(down_path, showWarnings = FALSE, recursive = TRUE);

    down_cell <- sprintf('%s/%s', down_path, img_type);
    dir.create(down_cell, showWarnings = FALSE, recursive = TRUE);
    
    # log file to record downloaded images
    img_log <- sprintf('%s/image_%s_%s.txt', down_path, img_type, str_query);
    fid_img <- file(img_log, 'w');
    cmt <- 'query\tgene\tcategory\tID\timage_url\timage_file';	# image log file header
    writeLines(cmt, fid_img);
    
    # other log files
    main_log <- sprintf('%s/main_%s_%s.txt', down_path, img_type, str_query);	# overall information about searched results
    cl_log   <- sprintf('%s/info_%s_%s.txt', down_path, img_type, str_query);		# tissue information
    
    #// get retrieved information
    mtxDATA <- bLab.protein.atlas.data(query=query);
        if (is.null(mtxDATA)) {
	return(list(NULL, NULL, NULL));
    }

    num_row <- nrow(mtxDATA);
    col_cl  <- 6;		# column containing URL of cell line information
    
    #// writing main information
    write.table(mtxDATA, main_log, row.names=FALSE);
    cat('DONE!\n'); flush.console();
    
    #------------------------------------------------------------------------------
    # Downloading cell line data
    #------------------------------------------------------------------------------
    #// retrieving data from each category
    clDATA <- matrix(data=NA, nrow=num_row, ncol=7);	# cell line data
    clDATA[, 1] <- mtxDATA[, 1]				# clDATA and mtxDATA has same number of rows
    
    for (i in 1:num_row) {
	# if information exists for tissue then retrieving information
	cl_url <- mtxDATA[i, col_cl];
	if (!is.na(cl_url)) {
	    tmpDATA 	<- bLab.protein.atlas.cell_parse(cl_url);
	    clDATA[i, ] <- tmpDATA[[1]];	# update clDATA based on individual tissue URLs
	    if (!is.null(tmpDATA[[2]])) {
		img_url   <- tmpDATA[[2]];	# image full url
		img_label <- tmpDATA[[4]];	# image label 
		for (j in 1:length(img_url)) {
		    tmp_img_url  <- img_url[j];
		    tmp_name <- strsplit(img_url[j], '/')[[1]];
		    tmp_name <- tmp_name[length(tmp_name)];
		    
		    # file image name = tissue_path + image_label() + original name
		    tmp_label_name <- gsub('[/,: ]', '-', img_label[j]);
		    
		    tmp_img_name <- sprintf('%s/%s_%s_%s',down_cell, clDATA[i, 1], tmp_label_name, tmp_name);
		    cmt <- sprintf('%s\t%s\t%s\t%s\t%s\t%s', query,  clDATA[i, 1], img_type, tmp_label_name, tmp_img_url, tmp_img_name);
		    writeLines(cmt, fid_img);
		    #download.file(tmp_img_url, tmp_img_name, method='auto', quiet=(!verbose));
		    bLab.file.proof_download(tmp_img_url, tmp_img_name, max_itr=1000000, verbose=TRUE);
		}
		clDATA[i, 7] <- img_log;
	    }
	}
    }
    cell_cols   <- c('gene', 'description', 'rna_expression', 'protein_expression', 'protein_class', 'reliability', 'images');
    colnames(clDATA) <- cell_cols;
    write.table(clDATA, cl_log, row.names=FALSE);
    cat('DONE!\n'); flush.console();
    close(fid_img);
    imgDATA <- read.table(img_log, header = TRUE, sep = "\t",);
    
    #// return all values
    return(list(mtxDATA, clDATA, imgDATA));
}




#********************************************************************
# Search protein atlas DB and download cell line images including
# antibody stained images and corresponding meta images
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.cell_full(query=NULL, down_home='./');
#
# ::Parameters::
#	query     (string): search query (e.g., insulin, PGR, CD36, SCARB2)
#	down_home (string): home directory where files can be downloaded
#		(e.g., query='CD36', down_home='/my/download' will create directory with 'query'
#                so, file will be stored in '/my/download/CD36/...');
#	verbose  (boolean): TRUE= display downloading status, FALSE= do not show downloading status
#
# ::Outputs::
#	outDATA	(list): containing all information about the searched results
#		outDATA[[1]]: main information about search results
#		outDATA[[2]]: details about cell line and corresponding genes
#		outDATA[[3]]: downloaded cell line image information
#		outDATA[[4]]: dataframe containing downloaded stained-images and meta information
#
#	TEXT files: downloaded information in text file
#		/down_home/<query>/main_cell_<query>.txt : tab delimited text file of outDATA[[1]]
#		/down_home/<query>/image_cell_<query>.txt: tab delimited text file of downloaded images
#		/down_home/<query>/info_cell_<query>.txt : tab delimited text file of retrieved tissue information
#		/down_home/<query>/stain_cell_<query>.txt      : overall summaries of cell lines that belong to a gene
#		/down_home/<query>/stain_info_cell_<query>.txt : antibody stained images and meta information
#		/down_home/<query>/stain_img_cell_<query>.txt  : all information is same as 'stain_info_cell_<query>.txt'
#								 addintionally containing the path of downloaded images
#
#********************************************************************
bLab.protein.atlas.cell_full <- function (query=NULL, down_home='./', verbose=TRUE) {
   
    #------------------------------------------------------------------------------
    # configuration for downloading environment
    #------------------------------------------------------------------------------
    # defining path for downloaded images
    img_type <- 'cell';
    
    #down_path <- sprintf('%s/%s', down_home, query);
    str_query <- paste(strsplit(query, ' ')[[1]], collapse='_');
    down_path <- sprintf('%s/%s', down_home, str_query);
    
    dir.create(down_path, showWarnings = FALSE, recursive = TRUE);

    down_cell <- sprintf('%s/%s', down_path, img_type);
    dir.create(down_cell, showWarnings = FALSE, recursive = TRUE);
    
    # log file to record downloaded images
    img_log <- sprintf('%s/image_%s_%s.txt', down_path, img_type, str_query);
    fid_img <- file(img_log, 'w');
    cmt <- 'query\tgene\tcategory\tID\timage_url\timage_file';	# image log file header
    writeLines(cmt, fid_img);
    
    # other log files
    main_log <- sprintf('%s/main_%s_%s.txt', down_path, img_type, str_query);		# overall information about searched results
    cl_log   <- sprintf('%s/info_%s_%s.txt', down_path, img_type, str_query);		# tissue information
    
    stain_log      <- sprintf('%s/stain_%s_%s.txt', down_path, img_type, str_query);	# stain information
    if (file.exists(stain_log)) {
	file.remove(stain_log);
    }
    
    stain_info_log <- sprintf('%s/stain_info_%s_%s.txt', down_path, img_type, str_query);	# stain information
    if (file.exists(stain_info_log)) {
	file.remove(stain_info_log);
    }

    stain_img_log  <- sprintf('%s/stain_img_%s_%s.txt', down_path, img_type, str_query);	# stain information
    if (file.exists(stain_img_log)) {
	file.remove(stain_img_log);
    }
    

    #// get retrieved information
    mtxDATA <- as.matrix(bLab.protein.atlas.data(query=query));
        if (is.null(mtxDATA)) {
	return(list(NULL, NULL, NULL));
    }

    num_row <- nrow(mtxDATA);
    col_cl  <- 6;		# column containing URL of cell line information
    
    #// writing main information
    write.table(mtxDATA, main_log, row.names=FALSE);
    cat('DONE!\n'); flush.console();
    
    #------------------------------------------------------------------------------
    # Downloading cell line data
    #------------------------------------------------------------------------------
    #// retrieving data from each category
    clDATA <- matrix(data=NA, nrow=num_row, ncol=7);	# cell line data
    clDATA[, 1] <- mtxDATA[, 1]				# clDATA and mtxDATA has same number of rows
    
    for (i in 1:num_row) {
	# if information exists for tissue then retrieving information
	cl_url <- mtxDATA[i, col_cl];
	if (!is.na(cl_url)) {
	    tmpDATA 	<- bLab.protein.atlas.cell_parse(cl_url);
	    clDATA[i, ] <- tmpDATA[[1]];	# update clDATA based on individual tissue URLs
	    if (!is.null(tmpDATA[[2]])) {
		img_url   <- tmpDATA[[2]];	# image full url
		img_label <- tmpDATA[[4]];	# image label 
		for (j in 1:length(img_url)) {
		    tmp_img_url  <- img_url[j];
		    tmp_name <- strsplit(img_url[j], '/')[[1]];
		    tmp_name <- tmp_name[length(tmp_name)];
		    
		    # file image name = tissue_path + image_label() + original name
		    tmp_label_name <- gsub("-", "_", img_label[j]);
		    
		    tmp_img_name <- sprintf('%s/%s_%s_%s',down_cell, clDATA[i, 1], tmp_label_name, tmp_name);
		    cmt <- sprintf('%s\t%s\t%s\t%s\t%s\t%s', query,  clDATA[i, 1], img_type, tmp_label_name, tmp_img_url, tmp_img_name);
		    writeLines(cmt, fid_img);
		    #download.file(tmp_img_url, tmp_img_name, method='auto', quiet=(!verbose));
		    bLab.file.proof_download(tmp_img_url, tmp_img_name, max_itr=1000000, verbose=TRUE);
		}
		clDATA[i, 7] <- img_log;
	    }
	    
	    # extracting cell line summary ['gene', 'cell line', 'cell line category', 'antibody', 'image url']; 
	    tmpDATA2 	<- bLab.protein.atlas.cell_summary_parse(cl_url);
	    write.table(tmpDATA2, file=stain_log, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t');
	    
	    # parsing stain page
	    cmt <- sprintf('\t>>> Retrieving high resolution image and meta tags...\n'); cat(cmt); flush.console();
	    s_num_row <- nrow(tmpDATA2);
	    s_num_col <- ncol(tmpDATA2);
	    for (s in 1:s_num_row) {
		cl_url2  <- tmpDATA2[s, s_num_col];
		#// Extracing stain image information
		#	column information
		#	[1]'gene',           [2]'cell_line_name', [3]'antibody',       [4]'gender', [5]'age',
		#       [6]'cell_line_id',   [7]'number_of_cell', [8]'positive_cells', [9]'img_url'
		imgStain <- bLab.protein.atlas.cell_stain_parse(cl_url2);
		if (!is.na(imgStain[1])) {
		    AntiBody <- matrix(data=tmpDATA2[s, 4], nrow=nrow(imgStain), ncol=1);
		    imgStain <- bLab.data.InsertCol(3, AntiBody, imgStain);
			cmt <- sprintf('[%d/%d] ncol(imgStain)=%d\n', s, s_num_row, ncol(imgStain));
			cat(cmt);
			print(cl_url2);
			print(imgStain);
			cat('\n'); flush.console();

		    write.table(imgStain, file=stain_info_log, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t');
		}
	    }
	}
    }
    
    #**********************************************************
    # Downloading high resolution (i.e., antibody stain) images
    #**********************************************************
    if (file.exists(stain_info_log)) {
	imgH      <- as.matrix(read.table(stain_info_log, sep='\t'));
	nrow_imgH <- nrow(imgH);
	ncol_imgH <- ncol(imgH);
	
	imgDATA2 <- matrix(data=NA, nrow=nrow(imgH), ncol=(ncol_imgH+1));
	imgDATA2[ ,1:ncol_imgH] <- as.matrix(imgH);		# adding one more column (i.e., path to store images) into imgH 
	
	for (i in 1:nrow_imgH) {
	    #// extracting image name
	    tmp_img_url  <- imgH[i, ncol_imgH];
	    tmp_img_name <- strsplit(tmp_img_url, '/')[[1]];
	    tmp_img_name <- tmp_img_name[length(tmp_img_name)];
	    
	    g_name   <- gsub('[/, ]', '-', imgH[i, 1]);
	    cell_line<- gsub('[/, ]', '-', imgH[i, 2]);
	    cell_line<- gsub('[/, ]', '-', cell_line);
	    antibody <- gsub('[/, ]', '-', imgH[i, 3]);
	    gender   <- gsub('[/, ]', '-', imgH[i, 4]);
	    tmp_img_name <- sprintf('%s/%s_%s_%s_%s_%s', down_cell, g_name, cell_line, antibody, gender, tmp_img_name);
	    bLab.file.proof_download(tmp_img_url, tmp_img_name, max_itr=1000000, verbose=TRUE);
	    imgDATA2[i, ncol_imgH+1] <- tmp_img_name;
	}
	
	colNames <- c('gene',     'cell_line_name', 'antibody',       'gender',
		      'age',      'cell_line_id',   'number_of_cell', 'positive_cells',  
		      'img_url',  'image_file');
	colnames(imgDATA2) <- colNames;
	write.table(imgDATA2, file=stain_img_log, append=FALSE, row.names=FALSE, col.names=TRUE, sep='\t');
    } else {
	imgDATA2 <- NA;
    }
    
    cell_cols   <- c('gene', 'description', 'rna_expression', 'protein_expression', 'protein_class', 'reliability', 'images');
    colnames(clDATA) <- cell_cols;
    write.table(clDATA, cl_log, row.names=FALSE);
    cat('DONE!\n'); flush.console();
    close(fid_img);
    imgDATA <- read.table(img_log, header = TRUE, sep = "\t",);

    #// return all values
    return(list(mtxDATA, clDATA, imgDATA, imgDATA2));
}




#********************************************************************
# Parsing cancer category in protein ATLAS 
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.cancer_parse(cc_url);
#
# ::Parameters::
#	cc_url    (string): url for cancer webpage (e.g., http://www.proteinatlas.org/ENSG00000169047-IRS1/cancer)
#
# ::Outputs::
#	outDATA	(list): containing all information about the searched results
#		outDATA[[1]]: cancer information 
#		outDATA[[2]]: retrieved complete url for images
#		outDATA[[3]]: retrieved short url for images
#		outDATA[[4]]: retrieved image labels
#		outDATA[[5]]: staining summary section (list of cancers)
#
#********************************************************************
bLab.protein.atlas.cancer_parse <- function(cc_url, home_url='http://www.proteinatlas.org') {
    #// store data into a matrix 
    ccDATA <- matrix(data=NA, nrow=1, ncol=5);
    
    #// searching tissue URL and get informaion 
    rPage <- htmlTreeParse(cc_url, useInternalNodes=TRUE);
    
    #// varifying missing data & get information
    cancer_cols   <- c('gene', 'description', 'protein_class', 'protein_evidence', 'images');
    cancer_groups <- c("Gene description", "Protein class", "Protein evidence")
    
    ccDATA[1, 1]  <- bLab.file.trim(xpathSApply(rPage, "//div[@class='atlas_header']/div[@class='gene_name']", xmlValue));
    Group 	  <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark staticwidth']/tr//th", xmlValue));
    Info  	  <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark staticwidth']/tr//td", xmlValue));
    for (g in 1:length(Group)) {
	info_idx <- which(tolower(cancer_groups) == tolower(Group[g]));
	ccDATA[1, info_idx+1] <- Info[g]
    }
    
    colnames(ccDATA) <- cancer_cols;
    
    #// Extracting image URLs
    imgDATA <- NULL;
    img_url   <- xpathSApply(rPage, "//table[@style='margin: auto; margin-bottom: 5px;']//td[@style='padding-right: 10px; padding-left: 10px; background: #C2C2C2; text-align: center;']//a", xmlGetAttr, 'href');
    img_label <- bLab.file.trim(xpathSApply(rPage, "//table[@style='margin: auto; margin-bottom: 5px;']//td[@style='padding-right: 10px; padding-left: 10px; background: #C2C2C2; text-align: center;']", xmlValue));
    
    if (length(img_url) > 0) {
	for (j in 1:length(img_url)) {
	    imgDATA[j]  <- sprintf('%s%s', home_url, img_url[j]);
	}
    }
    
    #--------------------------------------------
    # for STAINING SUMMARY section
    #--------------------------------------------
    ccSTAI     <- matrix(data=NA, nrow=1, ncol=3);  # matrix for STAINING SUMMARY
    st_summary <- xpathSApply(rPage, "//table[@class='border dark']/tr//td[@style='width:100%;']", xmlValue);
    cc_tag     <- "//tbody[@onclick]/tr//td[@class='nowrap tight'][1]";
    cancers    <- xpathSApply(rPage, cc_tag, xmlValue);

    return(list(ccDATA, imgDATA, img_url, img_label, cancers));
}



#********************************************************************
# Parsing cancer category in protein ATLAS 
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.cancer_stain_parse(cc_url);
#
# ::Parameters::
#	cc_url    (string): url for cancer webpage (e.g., http://www.proteinatlas.org/ENSG00000082175-PGR/cancer/tissue)
#
# ::Outputs::
#	outDATA	(list): containing all information about the searched results
#		outDATA[[1]]: cancer information per each tissue (i.e, row), columns corresponds to 'antibody stains'
#			      last column contains URL for each tissue
#		outDATA[[2]]: cancer information per each tissue (i.e, row), columns corresponds to 'total number of antibody stains' 
#		outDATA[[3]]: summary and literature conformity 
#
#********************************************************************
bLab.protein.atlas.cancer_stain_parse <- function(cc_url) {
    #// searching tissue URL and get informaion
    tmp   <- strsplit(cc_url, '/')[[1]];
    p_url <- paste(tmp[1:length(tmp)-1], collapse='/');# parent url of current url
    rPage <- htmlTreeParse(cc_url, useInternalNodes=TRUE);
    
    gname  <- bLab.file.trim(xpathSApply(rPage, "//div[@class='atlas_header']/div[@class='gene_name']", xmlValue));
    stains <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark no_styled_a tight']/tr//p[@class='last']", xmlValue));
    cancers<- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark no_styled_a tight']/tbody[@class='hover']//tr//td[2]", xmlValue));
    
    #// define full url for each cancers
    tmp_url <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark no_styled_a tight']/tbody[@class='hover']//tr//td[2]//a", xmlGetAttr, 'href'));
    cc_url2 <- array();
    for (u in 1:length(tmp_url)) {
	cc_url2[u] <- sprintf('%s/%s', p_url, tmp_url[u]);
    }
    
    #// extracting distribution of gene expression [# of high, # of medium, # of low ]
    tmp_gexp <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark no_styled_a tight']/tbody[@class='hover']//div[@class='bar']//img", xmlGetAttr, 'onmouseover'));
    gExp   <- array();
    numExp <- array();		# total number of expressed samples
    exp_cnt <- 0;
    for (i in 1:length(tmp_gexp)) {
	exp_cnt <- exp_cnt + 1;
	
	tmp <- strsplit(tmp_gexp[i], ':')[[1]]; 			# split ':'
	tmp1 <- bLab.file.trim(strsplit(tmp[2], '\'')[[1]][1]);		# extract x out of yy (High, Medium, Low, Not detected);
	tmp2 <- as.numeric(strsplit(gsub('\\D', ' ', tmp1), ' ')[[1]])
	no_na_idx <- bLab.data.findNA(tmp2)[[2]];	# find non 'NA' index
	nums <- tmp2[no_na_idx];
	
	#// reformating expression
	exv    <- nums[1];					# expression value
	tot    <- nums[2];					# total number of staining samples
	
	#// the variable 'aff' is just for debugging purpose
	if (grepl('high', tmp[2], ignore.case = TRUE)) {
	    aff <- 'H';
	} else if (grepl('medium', tmp[2], ignore.case = TRUE)) {
	    aff <- 'M';
	} else if (grepl('low', tmp[2], ignore.case = TRUE)) {
	    aff <- 'L';
	} else if (grepl('not', tmp[2], ignore.case = TRUE)) {
	    aff <- 'NA';
	} else {
	    stop('bLab.protein.atlas.cancer_stain_parse:: staining affinity is missed..');
	}
	
	#// get toal number of samples and expressed samples
	numExp[exp_cnt] <- tot;
	gExp[exp_cnt]   <- exv;
    }
    
    mtx_gExp   <- matrix(gExp,   ncol=length(stains)*4, byrow=TRUE);	# matrix containing expressed distribution
    mtx_numExp <- matrix(numExp, ncol=length(stains)*4, byrow=TRUE);	# matrix containing totol number of samples
    
    #// extraing summary parts
    info   <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark no_styled_a tight']//th[@colspan='2']//p", xmlValue));
    desc   <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark no_styled_a tight']//th[@colspan='2']//following-sibling::td", xmlValue));
    
    
    stainNames <- array();
    for (i in 1:length(stains)) {
	stainNames[i] <- tolower(paste(strsplit(stains[i], ' ')[[1]], collapse='_'));
    }

    mtxDUMMY<- matrix(desc, ncol=length(stains), byrow=TRUE);
    mtxINFO<- matrix(data=NA, nrow=length(stains), ncol=3);
    s_cnt <- 0;
    for (i in 1:length(stains)) {
	mtxINFO[i, 1] <- stains[i];
	mtxINFO[i, 2] <- mtxDUMMY[1, i];
	mtxINFO[i, 3] <- mtxDUMMY[2, i];
    }
    colnames(mtxINFO) <- c('antibody', 'summary','literature_conformity');
    
    
    #// merging data into a matrix
    # gene, stains1~n, url
    mtxDATA <- matrix(data=NA, nrow=length(cancers), ncol=(length(stains) * 4 + 3));
    mtxTOTA <- matrix(data=NA, nrow=length(cancers), ncol=(length(stains) * 4 + 2));
    for (i in 1:nrow(mtxDATA)) {
	#// for data
	mtxDATA[i, 1] <- gname;
	mtxDATA[i, 2] <- cancers[i];
	mtxDATA[i, 3:(ncol(mtx_gExp)+2)] <- mtx_gExp[i, ]; 
	mtxDATA[i, ncol(mtxDATA)] <- cc_url2[i];
	
	#// for information about total number of samples
	mtxTOTA[i, 1] <- gname;
	mtxTOTA[i, 2] <- cancers[i];
	mtxTOTA[i, 3:(ncol(mtx_numExp)+2)] <- mtx_numExp[i, ]; 
    }
    
    #// generating column names
    affLabel <- c('high', 'medium', 'low', 'na');
    colNames <- array();
    colNames[1] <- 'gene';
    colNames[2] <- 'tissue';
    
    c_idx <- 2;
    for (i in stainNames) {
	for (j in affLabel) {
	    c_idx <- c_idx + 1;
	    tmp <- sprintf('%s_%s', i, j);
	    colNames[c_idx] <- tmp;
	}
    }
    colNames[length(colNames)+1] <- 'image_url';
    
    #// defining column names
    colnames(mtxDATA) <- colNames;
    colnames(mtxTOTA) <- colNames[1:(length(colNames)-1)];
    
    return (list(mtxDATA, mtxTOTA, mtxINFO));
}    




#********************************************************************
# Parsing cancer category in protein ATLAS 
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.cancer_stain_image_parse(cc_url);
#
# ::Parameters::
#	cc_url    (string): url for cancer webpage (e.g., http://www.proteinatlas.org/ENSG00000254647-INS/cancer/tissue/colorectal+cancer)
#
# ::Outputs::
#	outDATA	(matrix): containing all information about images
#	column information
#	[1]'gene',    [2]'antibody',          [3]'patient_id', [4] 'gender',  [5] 'age'
#	[6]'tissue',  [7]'antibody_staining', [8]'intensity',  [9]'quantity', [10]'location',
#	[11]'img_url'
#
#********************************************************************
bLab.protein.atlas.cancer_stain_image_parse <- function(cc_url, home_url='http://www.proteinatlas.org') {
    
    rPage <- htmlTreeParse(cc_url, useInternalNodes=TRUE);
    
    gname  <- bLab.file.trim(xpathSApply(rPage, "//div[@class='atlas_header']/div[@class='gene_name']", xmlValue));
    antiB  <- bLab.file.trim(xpathSApply(rPage, "//table[@class='border dark']/tr//th[@class='head']/p[@class='last']", xmlValue));

    image_url  <- bLab.file.trim(xpathSApply(rPage, "//table[@class='noborder']//div[@name='annotationPair']/a", xmlGetAttr, 'href'));
    image_info <- bLab.file.trim(xpathSApply(rPage, "//table[@class='noborder']//div[@name='annotationPair']//a/img", xmlGetAttr, 'onmouseover'));
    
    #// return NA if neither image_url nor image_info exists
    if ((length(image_url) < 1) || (length(image_info) < 1))  {
	return(NA);
    }
    
    #// get image information 
    imgInfo  <- matrix(data=NA, nrow=length(image_url), ncol=12);
    colNames <- c('gene',   'disease', 'antibody',          'patient_id', 'gender',   'age',
		  'tissue', 'antibody_staining', 'intensity',  'quantity', 'location', 'img_url'); 


    for (i in 1:length(image_info)) {
	#// gene name
	imgInfo[i, 1] <- gname;
	disease <- strsplit(cc_url, '/')[[1]];
	disease <- gsub('[+]', ' ', disease[length(disease)]);
	imgInfo[i, 2]  <- disease;
	
	#// get name of corresponding antibody
	sp_url        <- strsplit(image_url[i], '/')[[1]];
	anti_id       <- sp_url[3];
	anti_idx      <- bLab.data.grep(anti_id, antiB)[[2]];
	imgInfo[i, 3] <- antiB[anti_idx];
	
	#*** Split information based on <br> tag ****
	tmpInfo <- strsplit(image_info[i], '<br>')[[1]];
	
	#cmt <- sprintf('*** PARSING tmpINFO ***\n\t==>%s\n\n', image_info[i]); cat(cmt); flush.console();
	
	#// patient_id
	p_id <- grep('Patient id', tmpInfo, ignore.case=TRUE)
	if (length(p_id) > 0) {
	    imgInfo[i, 4] <- bLab.file.trim(bLab.data.rm_html(strsplit(tmpInfo[p_id], ':')[[1]][2]));
	}
	
	#// gender information
	if (length(grep('Female', tmpInfo)) > 0) {
	    imgInfo[i, 5] <- 'Female';
	} else if (length(grep('Male', tmpInfo)) > 0) {
	    imgInfo[i, 5] <- 'Male';
	}
	
	#// age information
	age_idx <- grep('age', tmpInfo, ignore.case=TRUE);
	#cmt <- '\t==>::age_idx:: \n'; cat(cmt); flush.console();
	#print(tmpInfo[age_idx]); flush.console();
	if (length(age_idx) > 0) {
	    age_idx <- age_idx[1];	# age is always located at the first line!
	    imgInfo[i, 6] <- bLab.file.trim(gsub('.* ([0-9]+)<.*', '\\1',  tmpInfo[age_idx]));
	}
	
	#// tissue information
	tss_idx <- grep('\\(.*\\)', tmpInfo);
	if (length(tss_idx) > 0) {
	    imgInfo[i, 7] <- bLab.file.trim(bLab.data.rm_html(paste(tmpInfo[tss_idx], collapse=';')));
	}

	#// antibody staining information
	as_idx <- grep('Antibody staining', tmpInfo, ignore.case=TRUE)
	if (length(as_idx) > 0) {
	    imgInfo[i, 8] <- bLab.file.trim(bLab.data.rm_html(strsplit(tmpInfo[as_idx], ':')[[1]][2]));
	}

	#// Intensity information
	it_idx <- grep('Intensity', tmpInfo, ignore.case=TRUE)
	if (length(it_idx) > 0) {
	    imgInfo[i, 9] <- bLab.file.trim(bLab.data.rm_html(strsplit(tmpInfo[it_idx], ':')[[1]][2]));
	}

	#// Quantity information
	qt_idx <- grep('Quantity', tmpInfo, ignore.case=TRUE)
	if (length(qt_idx) > 0) {
	    imgInfo[i, 10] <- bLab.file.trim(bLab.data.rm_html(strsplit(tmpInfo[qt_idx], ':')[[1]][2]));
	}

	#// Location information
	lc_idx <- grep('Location', tmpInfo, ignore.case=TRUE)
	if (length(lc_idx) > 0) {
	    tmp_lc <- bLab.file.trim(bLab.data.rm_html(strsplit(tmpInfo[lc_idx], ':')[[1]][2]));
	    imgInfo[i, 11] <- bLab.file.trim(gsub("(.*)'.*", '\\1', tmp_lc));
	}

	#// image URL
	imgInfo[i,12] <- sprintf('%s%s', home_url, image_url[i]);
    }
    colnames(imgInfo) <- colNames;
    return (imgInfo);
}

#********************************************************************
# Search protein atlas DB and download cancer images
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.cancer(query=NULL, down_home='./');
#
# ::Parameters::
#	query     (string): search query (e.g., insulin, PGR, CD36)
#	down_home (string): home directory where files can be downloaded
#		(e.g., query='CD36', down_home='/my/download' will create directory with 'query'
#                so, file will be stored in '/my/download/CD36/...');
#	verbose  (boolean): TRUE= display downloading status, FALSE= do not show downloading status
#
# ::Outputs::
#	outDATA	(list): containing all information about the searched results
#		outDATA[[1]]: main information about search results
#		outDATA[[2]]: details about cancer and corresponding genes
#		outDATA[[3]]: downloaded cancer image information
#
#	TEXT files: downloaded information in text file
#		/down_home/<query>/main_cancer_<query>.txt : tab delimited text file of outDATA[[1]]
#		/down_home/<query>/image_cancer_<query>.txt: tab delimited text file of downloaded images
#		/down_home/<query>/info_cancer_<query>.txt : tab delimited text file of retrieved tissue information
#
#********************************************************************
bLab.protein.atlas.cancer <- function (query=NULL, down_home='./', verbose=TRUE) {
    
    #------------------------------------------------------------------------------
    # configuration for downloading environment
    #------------------------------------------------------------------------------
    # defining path for downloaded images
    img_type <- 'cancer';
    
    #down_path <- sprintf('%s/%s', down_home, query);
    str_query <- paste(strsplit(query, ' ')[[1]], collapse='_');
    down_path <- sprintf('%s/%s', down_home, str_query);
    
    dir.create(down_path, showWarnings = FALSE, recursive = TRUE);

    down_cancer <- sprintf('%s/%s', down_path, img_type);
    dir.create(down_cancer, showWarnings = FALSE, recursive = TRUE);
    
    # log file to record downloaded images
    img_log <- sprintf('%s/image_%s_%s.txt', down_path, img_type, str_query);
    fid_img <- file(img_log, 'w');
    cmt <- 'query\tgene\tcategory\tloci\timage_url\timage_file';	# image log file header
    writeLines(cmt, fid_img);
    
    # other log files
    main_log <- sprintf('%s/main_%s_%s.txt', down_path, img_type, str_query);	# overall information about searched results
    cc_log   <- sprintf('%s/info_%s_%s.txt', down_path, img_type, str_query);	# cancer information
    
    #// get retrieved information
    mtxDATA <- bLab.protein.atlas.data(query=query);
        if (is.null(mtxDATA)) {
	return(list(NULL, NULL, NULL));
    }

    num_row <- nrow(mtxDATA);
    col_cc  <- 7;		# column containing URL of cancer information
    
    #// writing main information
    write.table(mtxDATA, main_log, row.names=FALSE);
    cat('DONE!\n'); flush.console();
    
    #------------------------------------------------------------------------------
    # Downloading cancer data
    #------------------------------------------------------------------------------
    #// retrieving data from each category
    ccDATA <- matrix(data=NA, nrow=num_row, ncol=5);	# tissue data
    ccDATA[, 1] <- mtxDATA[, 1]				# tsDATA and mtxDATA has same number of rows
    
    for (i in 1:num_row) {
	# if information exists for cancer then retrieving information
	cc_url <- mtxDATA[i, col_cc];
	if (!is.na(cc_url)) {
	    tmpDATA 	<- bLab.protein.atlas.cancer_parse(cc_url);
	    ccDATA[i, ] <- tmpDATA[[1]];	# update ccDATA based on individual cancer URLs
	    if (!is.null(tmpDATA[[2]])) {
		img_url   <- tmpDATA[[2]];	# image full url
		img_label <- tmpDATA[[4]];	# image label 
		for (j in 1:length(img_url)) {
		    tmp_img_url  <- img_url[j];
		    tmp_name <- strsplit(img_url[j], '/')[[1]];
		    tmp_name <- tmp_name[length(tmp_name)];
		    
		    # file image name = cancer_path + image_label() + original name 
		    tmp_img_name <- sprintf('%s/%s_%s_%s',down_cancer, ccDATA[i, 1], gsub(" ", "_", img_label[j]), tmp_name);
		    cmt <- sprintf('%s\t%s\t%s\t%s\t%s\t%s', query,  ccDATA[i, 1], img_type, img_label[j], tmp_img_url, tmp_img_name);
		    writeLines(cmt, fid_img);
		    #download.file(tmp_img_url, tmp_img_name, method='auto', quiet=(!verbose));
		    bLab.file.proof_download(tmp_img_url, tmp_img_name, max_itr=1000000, verbose=TRUE);
		}
		ccDATA[i, 5] <- img_log;
	    }
	}
    }
    cancer_cols   <- c('gene', 'description', 'protein_class', 'protein_evidence', 'images');
    colnames(ccDATA) <- cancer_cols;
    write.table(ccDATA, cc_log, row.names=FALSE);
    cat('DONE!\n'); flush.console();
    close(fid_img);
    imgDATA <- read.table(img_log, header = TRUE, sep = "\t",);
    
    #// return all values
    return(list(mtxDATA, ccDATA, imgDATA));
}



#********************************************************************
# Search protein atlas DB and download cancer images including stain images
# This requires more time than 'bLab.protein.atlas.cancer' which downloads
# represenative small images only.
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.cancer_full(query=NULL, down_home='./');
#
# ::Parameters::
#	query     (string): search query (e.g., insulin, PGR, CD36)
#	down_home (string): home directory where files can be downloaded
#		(e.g., query='CD36', down_home='/my/download' will create directory with 'query'
#                so, file will be stored in '/my/download/CD36/...');
#	verbose  (boolean): TRUE= display downloading status, FALSE= do not show downloading status
#
# ::Outputs::
#	outDATA	(list): containing all information about the searched results
#		outDATA[[1]]: main information about search results
#		outDATA[[2]]: details about cancer and corresponding genes
#		outDATA[[3]]: downloaded cancer image information
#		outDATA[[4]]: downloaded high resolution cancer images
#
#	TEXT files: downloaded information in text file
#		/down_home/<query>/main_cancer_<query>.txt : tab delimited text file of outDATA[[1]]
#		/down_home/<query>/image_cancer_<query>.txt: tab delimited text file of downloaded images
#		/down_home/<query>/info_cancer_<query>.txt : tab delimited text file of retrieved tissue information
#		/down_home/<query>/stain_cancer_<query>.txt : tab delimited text file of overall cancer data of
#							      (e.g., http://www.proteinatlas.org/ENSG00000170989-S1PR1/cancer/tissue)
#		/down_home/<query>/stain_info_cancer_<query>.txt : tab delimited text file of antibody stain information
#		/down_home/<query>/stain_img_cancer_<query>.txt : tab delimited text file of downloaded images 
#
#********************************************************************
bLab.protein.atlas.cancer_full <- function (query=NULL, down_home='./', verbose=TRUE) {
    
    #------------------------------------------------------------------------------
    # configuration for downloading environment
    #------------------------------------------------------------------------------
    # defining path for downloaded images
    img_type <- 'cancer';
    
    #down_path <- sprintf('%s/%s', down_home, query);
    str_query <- paste(strsplit(query, ' ')[[1]], collapse='_');
    down_path <- sprintf('%s/%s', down_home, str_query);
    
    dir.create(down_path, showWarnings = FALSE, recursive = TRUE);

    down_cancer <- sprintf('%s/%s', down_path, img_type);
    dir.create(down_cancer, showWarnings = FALSE, recursive = TRUE);
    
    # log file to record downloaded images
    img_log <- sprintf('%s/image_%s_%s.txt', down_path, img_type, str_query);
    fid_img <- file(img_log, 'w');
    cmt <- 'query\tgene\tcategory\tloci\timage_url\timage_file';	# image log file header
    writeLines(cmt, fid_img);
    
    # other log files
    main_log    <- sprintf('%s/main_%s_%s.txt', down_path, img_type, str_query);		# overall information about searched results
    cc_log      <- sprintf('%s/info_%s_%s.txt', down_path, img_type, str_query);		# cancer information
    
    stain_log      <- sprintf('%s/stain_%s_%s.txt', down_path, img_type, str_query);		# stain information
    if (file.exists(stain_log)) {
	file.remove(stain_log);
    }
    
    stain_info_log <- sprintf('%s/stain_info_%s_%s.txt', down_path, img_type, str_query);	# stain information
    if (file.exists(stain_info_log)) {
	file.remove(stain_info_log);
    }

    stain_img_log  <- sprintf('%s/stain_img_%s_%s.txt', down_path, img_type, str_query);	# stain information
    if (file.exists(stain_img_log)) {
	file.remove(stain_img_log);
    }
    
    #// get retrieved information
    mtxDATA <- bLab.protein.atlas.data(query=query);
        if (is.null(mtxDATA)) {
	return(list(NULL, NULL, NULL));
    }

    num_row <- nrow(mtxDATA);
    col_cc  <- 7;		# column containing URL of cancer information
    
    #// writing main information
    write.table(mtxDATA, main_log, row.names=FALSE);
    
    #------------------------------------------------------------------------------
    # Downloading cancer data
    #------------------------------------------------------------------------------
    #// retrieving data from each category
    ccDATA <- matrix(data=NA, nrow=num_row, ncol=5);	# tissue data
    ccDATA[, 1] <- mtxDATA[, 1]				# tsDATA and mtxDATA has same number of rows
    
    for (i in 1:num_row) {
	# if information exists for cancer then retrieving information
	cc_url <- mtxDATA[i, col_cc];
	if (!is.na(cc_url)) {
	    tmpDATA 	<- bLab.protein.atlas.cancer_parse(cc_url);
	    ccDATA[i, ] <- tmpDATA[[1]];	# update ccDATA based on individual cancer URLs
	    if (!is.null(tmpDATA[[2]])) {
		img_url   <- tmpDATA[[2]];	# image full url
		img_label <- tmpDATA[[4]];	# image label 
		for (j in 1:length(img_url)) {
		    tmp_img_url  <- img_url[j];
		    tmp_name <- strsplit(img_url[j], '/')[[1]];
		    tmp_name <- tmp_name[length(tmp_name)];
		    
		    # file image name = cancer_path + image_label() + original name 
		    tmp_img_name <- sprintf('%s/%s_%s_%s',down_cancer, ccDATA[i, 1], gsub(" ", "_", img_label[j]), tmp_name);
		    cmt <- sprintf('%s\t%s\t%s\t%s\t%s\t%s', query,  ccDATA[i, 1], img_type, img_label[j], tmp_img_url, tmp_img_name);
		    writeLines(cmt, fid_img);
		    #download.file(tmp_img_url, tmp_img_name, method='auto', quiet=(!verbose), mode='wb');
		    bLab.file.proof_download(tmp_img_url, tmp_img_name, max_itr=1000000, verbose=TRUE);
		}
		ccDATA[i, 5] <- img_log;
	    }
	    
	    #***********************************************************************************************
	    # preparing for donwloading high resolution (i.e., antibody stain) images with antibody staining
	    #***********************************************************************************************
	    cc_url2  <- sprintf('%s/tissue', cc_url);
	    cmt <- sprintf('Reading %s...', cc_url2); cat(cmt); flush.console();
	    # parsing stain page
	    myStain <- bLab.protein.atlas.cancer_stain_parse(cc_url2);
	    write.table(myStain[[1]], file=stain_log, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t');
	    write.table(myStain[[3]], file=stain_info_log, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t');
	    
	    imgURLs <- myStain[[1]][, ncol(myStain[[1]])];
	    for (st_img in imgURLs) {
		myImage <- bLab.protein.atlas.cancer_stain_image_parse(st_img);
		if (!is.na(myImage[1])) {
		   write.table(myImage, file=stain_img_log, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t');	
		}
	    }
	    cmt <- sprintf('DONE!\n'); cat(cmt); flush.console();
	}
    }
    
    #**********************************************************
    # Downloading high resolution (i.e., antibody stain) images
    #**********************************************************
    if (file.exists(stain_img_log)) {
	imgH      <- as.matrix(read.table(stain_img_log, sep='\t'));
	nrow_imgH <- nrow(imgH);
	ncol_imgH <- ncol(imgH);
	
	imgDATA2 <- matrix(data=NA, nrow=nrow(imgH), ncol=(ncol_imgH+1));
	imgDATA2[ ,1:ncol_imgH] <- as.matrix(imgH);		# adding one more column (i.e., path to store images) into imgH 
	
	for (i in 1:nrow_imgH) {
	    #// extracting image name
	    tmp_img_url  <- imgH[i, ncol_imgH];
	    tmp_img_name <- strsplit(tmp_img_url, '/')[[1]];
	    tmp_img_name <- tmp_img_name[length(tmp_img_name)];
	    
	    g_name   <- gsub('[/,: ]', '-', imgH[i, 1]);
	    disease  <- gsub('[/,: ]', '-', imgH[i, 2]);
	    antibody <- gsub('[/,: ]', '-', strsplit(imgH[i, 3], ' ')[[1]][2]);
	    gender   <- gsub('[/,: ]', '-', imgH[i, 5]);
	    tissue   <- gsub('[/,: ]', '-', strsplit(imgH[i, 7], ' ')[[1]][1]);
	    tmp_img_name <- sprintf('%s/%s_%s_%s_%s_%s_%s', down_cancer, g_name, disease, antibody, gender, tissue, tmp_img_name);
	    #download.file(tmp_img_url, tmp_img_name, method='auto', quiet=(!verbose), mode='wb');
	    bLab.file.proof_download(tmp_img_url, tmp_img_name, max_itr=1000000, verbose=TRUE);
	    imgDATA2[i, ncol_imgH+1] <- tmp_img_name;
	}
	
	colNames <- c('gene',     'disease', 'antibody',          'patient_id', 'gender',
		      'age',      'tissue',  'antibody_staining', 'intensity',  'quantity',
		      'location', 'img_url', 'image_file');
	colnames(imgDATA2) <- colNames;
	write.table(imgDATA2, file=stain_img_log, append=FALSE, row.names=FALSE, col.names=TRUE, sep='\t');

    } else {
	imgDATA2 <- NA;
    }
    
    cancer_cols   <- c('gene', 'description', 'protein_class', 'protein_evidence', 'images');
    colnames(ccDATA) <- cancer_cols;
    write.table(ccDATA, cc_log, row.names=FALSE, col.names=TRUE);
    cat('DONE!\n'); flush.console();
    close(fid_img);
    imgDATA <- read.table(img_log, header = TRUE, sep = "\t",);
    
    #// return all values
    return(list(mtxDATA, ccDATA, imgDATA, imgDATA2));
}



#********************************************************************
# Searching data downloaded from atlas.cancer_full
# 
# ::Usage::
# 	outDATA <- bLab.protein.atlas.cancer_stain_search(keyword, field, inFile);
#
# ::Parameters::
#	query     (string): search keyword (e.g., 'breast cancer')
#	fields    (string): fields where keyword to be searched (e.g., 'disease')
#			    [1] gene          	[2] disease       	[3] antibody
#			    [4] patient_id    	[5] gender        	[6] age
#			    [7] tissue        	[8] staining_rate 	[9] intensity
#			    [10] quantity      	[11] location      	[12] image_url
#			    [13] image_file
#	inFile   (string): file location where stain image information is located (i.g., /<query>/stain_img_cancer_<query>.txt)
#
# ::Outputs::
#	outDATA	(list): containing information about the matched data and its original index in the file
#		outDATA[[1]]: matched data
#		outDATA[[2]]: original row index of the matched data
#
#********************************************************************
bLab.protein.atlas.cancer_stain_search <- function(keyword, field, inFile) {
    # define a field to search query
    Fields <- c('gene',     'disease',   'antibody',      'patient_id', 'gender',
	        'age',      'tissue',    'staining_rate', 'intensity',  'quantity',
		'location', 'image_url', 'image_file');
    
    fidx <- which(tolower(field) == Fields);
    if (length(fidx) < 1)  {
	cmt <- bLab.data.printArray(Fields, 3); cat('Possible Fields:\n', cmt); flush.console();
	stop('bLab.protein.atlas.cancer_search:: field does not existed!')
    }
    
    #// reading file
    mtxStain  <- read.table(inFile, sep='\t');
    
    refSearch <- mtxStain[, fidx];
    query  <- sprintf('*%s*', keyword);		# creating query for the search
    rstIdx <- grep(query, refSearch);		# index of matched items
    rstData<- as.matrix(mtxStain[rstIdx, ]);	# matched items
    
    return(list(rstData, rstIdx));    
}

