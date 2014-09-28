library(GEOquery)
library("R.utils")
source("/home/wim/Projects/GEOQuery/TransmartGEOAdapter.R")


#set the GEO study to load and the base output directory
geo_id <- "GSE58558"

outputDir = "/home/wim//Projects/compare_GEO_R_Scripts/bioConductor/"

#load the GSE
gseList <- getGEO(geo_id)

#create the output subdirectories and the outputFilePaths
geo_id <-paste(geo_id, "-noIndivExtraction",sep="")

#create the output subdirectories and the outputFilePaths
outputInfo <- createOutputDirsAndFiles(outputDir,geo_id )

#get the gse from the list by index or name
gse <- gseList$GSE58558_series_matrix.txt.gz
gse <- gseList[[1]]

# now get the phenotypic data (covariates etc.) using pData()
pd <- pData(gse)

#show the available phenotypic data
names(pd)

#take a subset of the phenotypic data we would like to use
pdSubset <- pd[c( "title", "geo_accession", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1","characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.4", "characteristics_ch1.5" )]
colnames(pdSubset) <- c( "individual", "sample_id", "organism", "time","lesional","responder", "age", "gender", "scorad" )


#remove the keys from the phenotypic columns by splitting the levels for the column on ": " and selecting index 2
#ie "time: day 1, time: day 2" becomes "day1, day 2"
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c("time","lesional",  "responder", "age", "gender", "scorad" ), ": ", 2)

levels(pdSubset$time) <- c("week_00","week_12", "week_02")

#select a further subset of the phenotypic data that describe the individuals. One column should uniquely identify the individuals
pdClinicalSubset = pdSubset[c( "individual", "organism", "time","lesional","responder", "age", "gender", "scorad"  )]

#set category codes. Under these category codes the phenotypic data will be listed in transmart
categoryCodes <- list()
categoryCodes$organism       <- "Subjects"
categoryCodes$time           <- "Subjects"
categoryCodes$lesional       <- "Subjects"
categoryCodes$responder      <- "Subjects"
categoryCodes$age            <- "Subjects"
categoryCodes$gender         <- "Subjects"
categoryCodes$scorad         <- "Subjects"



#select the individual sample mapping from the first subset
#add the platform_id and tissueType
#add time and lesional columns as the  attributes attributes
pdExpressionSubset = pdSubset[c( "individual", "sample_id" )]
pdExpressionSubset["PLATFORM"] <- attr(gse, "annotation")
pdExpressionSubset["TISSUETYPE"] <- as.character(pd$source_name_ch1)
pdExpressionSubset["ATTR1"] <- NA
pdExpressionSubset["ATTR2"] <- NA


#write the expressionTable
writeExpressionTable(outputInfo$outputFileExpression, gse)

#construct and write the clinicalTable
writeClinicalTable(outputInfo$outputFileClinical, pdClinicalSubset)

#construct and write the clinical table mapping
writeClinicalColumns(outputInfo$outputFileClinicalColumns,outputInfo$outputFileClinical,pdClinicalSubset, "individual", categoryCodes )

#construct and write the expression table mapping
writeExpressionMapping(outputInfo$outputFileExpressionMapping, pdExpressionSubset, geo_id  )

#construct and write clinicalParams file
writeClinicalParams(outputInfo$outputFileClinicalParams, outputInfo$outputFileClinicalColumns )

#construct and write expressionParams file
expressionValuesType = getExpressoinValueType(gse)  #try to detect expressionValuesType
writeExpressionParams(outputInfo$outputFileExpressionParams, geo_id, outputInfo$outputFileExpressionMapping, expressionValuesType)





























