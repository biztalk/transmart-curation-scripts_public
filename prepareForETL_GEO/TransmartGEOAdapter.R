#This R file contains convinience functions that can be used when downloading GEO files and converting them to files that can be loaded by transmart.
#These functions should be applicable to every dataset downloaded from GEO.
#To use these functions start your R script with  | source("/home/username/parentdirectory/TransmartGEOAdapter.R") |


##########
# write the expressionTable from a gse object to an output file
##########
writeExpressionTable <- function(outputFileExpression, gse){

    #get the expression data
    expressionTable =  exprs(gse)

    #add the rownames as the first ID_REF column
    expressionTable <- cbind(rownames(expressionTable),expressionTable)
    colnames(expressionTable)[1] <- "ID_REF"

    #write the expression table
    write.table(expressionTable, file= outputFileExpression, sep= "\t", row.names=FALSE, quote=FALSE )
}


##########
# write the phenotypic data of a gse object to an output file#
##########
writeClinicalTable <- function(outputFileClinical, pdSubset){

    #write the phenotypic table
    write.table(pdSubset, file=outputFileClinical , sep= "\t", row.names=FALSE, quote=FALSE, na="")
}

##########
# split the levels of certain columns in a dataframe on a certain splitString and take use an index of the splitresult as the new levels
#ie "time: day 1, time: day 2" becomes "day1, day 2" . Splitstring is ": ", index is 2 and columnNameList = (time) in this example
##########
splitLevelsOfDataFrameColums <- function(dataFrame, columnNameList, splitString, index ){

    for (columnName in columnNameList) {

        levels(dataFrame[[columnName]]) <- sapply(strsplit(levels(dataFrame[[columnName]]), splitString), function(x) x[index] )
    }
    return(dataFrame)
}

##########
# split the levels of certain columns in a dataframe on a certain splitString and take using the last element of the splitresult as the new levels
##########
splitLevelsOfDataFrameColumsTail <- function(dataFrame, columnNameList, splitString ){

    for (columnName in columnNameList) {

        levels(dataFrame[[columnName]]) <- sapply(strsplit(levels(dataFrame[[columnName]]), splitString), function(x) tail(x, n=1) )
    }
    return(dataFrame)
}



##########
# create column annotation for the clinical columns
# every column in the clinicalTable is annotated here with extra information
# every column in the clinicalTable is a record in the annotation dataFrame
# pdSubject should be the dataframe that is contains the individuals and the clinical information for every individuals
# the subjectIdColumn should be the name of the column that contains the unique identifiers for the individuals.
# subjectIdColumn is annotated with a special reserved categoryname and label
# categeory codes should be a named list that specifices a category code for each columnName. For example Subjects or Demographics#
##########
writeClinicalColumns <- function(outputFileClinicalColumns, outputFileClinical, pdSubset, subjectIdColumnName, categoryCodes){

    # set column names for mapping file. In these columns the extra annotation information is stored
    annotationColumnNames <- c("Filename", "Category_Code (tranSMART)", "Column", "DataLabel", "data_label_src", "ControlVocab_cd")

    #create the data frame. every record stand for a column in the clinicalTable
    annotation <- data.frame(matrix(nrow = ncol(pdSubset), ncol = length(annotationColumnNames)))
    colnames(annotation) <- annotationColumnNames
    rownames(annotation) <- colnames(pdSubset)

    annotation$Filename <- basename(outputFileClinical)

    # Loop over all the columns of the pdSubset and add  category code, data label and index to the annotation
    #if the columnName is the subjectIdColumnName set the subjectID values for categoryCode and dataLabel
    for (columnName in colnames(pdSubset)) {
        if(columnName == subjectIdColumnName) {
            annotation[columnName, "Category_Code (tranSMART)"] <- "Subjects"
            annotation[columnName, "DataLabel"] <- "SUBJ_ID"
        }else {
            annotation[columnName, "Category_Code (tranSMART)"] <- categoryCodes[[columnName]]
            annotation[columnName, "DataLabel"] <- columnName
        }
        #set the index of the column
        annotation[columnName, "Column"] <- match(columnName,colnames(pdSubset))
    }

    # remove NAs from last two columns
    annotation[["data_label_src"]] <- ""
    annotation[["ControlVocab_cd"]] <- ""

    # write column mapping table to file
    write.table(annotation, outputFileClinicalColumns, row.names = FALSE, sep = "\t", quote=FALSE)
}


##########
# create Column mapping file for expression data
# writes a sample to indivual mapping file. Indicates which samples in the expression data belong to which indivual
# a indivual can have one sample or multiple
# pdExpressionSubset arguments should be a dataframe with 6 columns
# 1) SUBJECT_ID
# 2) SAMPLE_ID
# 3) PLATFORM_ID of the mircroarray used
# 4) TISSUETYPE
# 5) ATTR1. For example a column with timepoints , ie day1 day2 day3 etc. Can also be NA if each indivual only has one sample
# 6) ATTR2. Same as ATTR1. For example a column with timepoints or tissue type , ie lesional skin, non lesional skin.

##########
writeExpressionMapping <- function(outputFileExpressionMapping, pdExpressionSubset, geo_id  )
{
    sampleCount = nrow(pdExpressionSubset)
    columnCount = ncol(pdExpressionSubset)

    #if ATTR1 and or ATTR2 are not completely NA use them in the category cd
    categoryCd <- "Biomarker_data+PLATFORM+TISSUETYPE"
    if(!all(is.na(pdExpressionSubset$ATTR1))){ categoryCd <- paste(categoryCd,"ATTR1", sep="+" )  }
    if(!all(is.na(pdExpressionSubset$ATTR2))){ categoryCd <- paste(categoryCd,"ATTR2", sep="+" )  }

    # set column names for mapping file
    mappingFileColNames <- c("STUDY_ID", "SITE_ID", "SUBJECT_ID", "SAMPLE_CD", "PLATFORM", "TISSUETYPE", "ATTR1",  "ATTR2",  "CATEGORY_CD")

    outColsMapping <- data.frame(matrix(nrow = sampleCount, ncol = length(mappingFileColNames)))
    colnames(outColsMapping) <- mappingFileColNames
    rownames(outColsMapping) <- pdExpressionSubset[[2]]
    outColsMapping$STUDY_ID <- geo_id
    outColsMapping$SUBJECT_ID <- pdExpressionSubset[[1]]
    outColsMapping$SAMPLE_CD <- pdExpressionSubset[[2]]
    outColsMapping$PLATFORM <- pdExpressionSubset[[3]]
    outColsMapping$TISSUETYPE <- pdExpressionSubset[[4]]
    outColsMapping$ATTR1 <- pdExpressionSubset[[5]]
    outColsMapping$ATTR2 <- pdExpressionSubset[[6]]
    outColsMapping$CATEGORY_CD <- categoryCd

    # cleanup NAs in SITE_ID column
    outColsMapping[["SITE_ID"]] <- ""

    # write column mapping table to file
    write.table(outColsMapping, outputFileExpressionMapping, row.names = FALSE, sep = "\t", quote=FALSE, na="")
}



writeClinicalParams <- function(outputFile, columnsFile, wordMapFileName="", recordExclusionFileName="" ) {

  # define clinical params
  clinicalParams <- paste(
    "COLUMN_MAP_FILE=", basename(columnsFile), "\r\n",
    "WORD_MAP_FILE=", wordMapFileName, "\r\n",
    "RECORD_EXCLUSION_FILE=", recordExclusionFileName,
    sep=""
  )
  write(clinicalParams, outputFile)
}



#the datatype argument defines in which space the expression values are. "R" means raw, "L" means base2 logspace and "T" means load the data as Z score
writeExpressionParams <- function(outputFile, geo_id, expressionMappingFileName, dataType="R" ) {

 # define expression params
 expressionParams <- paste(
   "DATA_FILE_PREFIX=",geo_id, "_expression", "\r\n",
   "DATA_TYPE=", dataType, "\r\n",
   "MAP_FILENAME=", basename(expressionMappingFileName),
   sep=""
 )
 write(expressionParams, outputFile)
}


#get the type of the expressionValues in the expressionTable
#test if the values are log transformed. If so returns "L"
#code to test for log transformatoin is copied from NCBI geo2r : http://www.ncbi.nlm.nih.gov/geo/geo2r/

#otherwise returns "R"
getExpressoinValueType <- function(gse)
{

    #get the expression data
    expressionTable =  exprs(gse)

    # log2 transform
    qx <- as.numeric(quantile(expressionTable, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
              (qx[6]-qx[1] > 50 && qx[2] > 0) ||
              (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

    #by default we set expressionValueType to log transformed
    expressionValueType <- "L"

    if(LogC){expressionValueType <- "R"}

    return(expressionValueType)
}





##########
# create output subdirectories and return the filePaths to outpufFiles in those subdirectories
##########
createOutputDirsAndFiles <- function(basedir, geo_id){


    outputDirGSE = file.path(outputDir, geo_id)
    outputDirGSEClinical = file.path(outputDirGSE, "clinical")
    outputDirGSEExpression = file.path(outputDirGSE, "expression")

    dir.create(outputDirGSE)
    dir.create(outputDirGSEClinical)
    dir.create(outputDirGSEExpression)

    outputFileClinicalPath = file.path(outputDirGSEClinical, paste(geo_id,  "_clinical.txt", sep=""))
    outputFileClinicalColumnsPath = file.path(outputDirGSEClinical, paste(geo_id,  "_columns.txt", sep=""))
    outputFileExpressionPath = file.path(outputDirGSEExpression, paste(geo_id,  "_expression.txt", sep=""))
    outputFileExpressionMappingPath = file.path(outputDirGSEExpression, paste(geo_id, "_mapping.txt", sep=""))
    outputFileClinicalParamsPath = file.path(outputDirGSE, "clinical.params")
    outputFileExpressionParamsPath = file.path(outputDirGSE, "expression.params")

    outputInfo = list(
                        outputFileClinical=outputFileClinicalPath,
                        outputFileClinicalColumns = outputFileClinicalColumnsPath,
                        outputFileExpression = outputFileExpressionPath,
                        outputFileExpressionMapping = outputFileExpressionMappingPath,
                        outputFileClinicalParams = outputFileClinicalParamsPath,
                        outputFileExpressionParams = outputFileExpressionParamsPath
     )
     return(outputInfo)
}
