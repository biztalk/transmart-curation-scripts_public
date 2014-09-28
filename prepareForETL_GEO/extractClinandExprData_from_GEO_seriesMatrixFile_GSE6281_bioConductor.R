#for  GEOquery installation instruction and documentatoin see the links below
#http://www.bioconductor.org/packages/release/bioc/html/GEOquery.html
#http://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.pdf

#load the GEOQUERY and R.utils library
library(GEOquery)
library("R.utils")
source("/home/wim/Projects/GEOQuery/TransmartGEOAdapter.R")


#set the GEO study to load and the base output directory
geo_id <- "GSE6281"
outputDir = "/home/wim//Projects/compare_GEO_R_Scripts/bioConductor/"

#create the output subdirectories and the outputFilePaths
outputInfo <- createOutputDirsAndFiles(outputDir, geo_id)

#load the GSE
gseList <- getGEO(geo_id)

#get the gse from the list by index or name
gse <- gseList[[1]]

# now get the phenotypic data (covariates etc.) using pData()
pd <- pData(gse)

#show the available phenotypic data
names(pd)
head(pd)

#select a subset of the phenotypic data
pdSubset = pd[c( "title","title", "title", "geo_accession", "source_name_ch1", "characteristics_ch1", "characteristics_ch1", "organism_ch1"  )]


#set the labels to use for the phenotypic data. This is how the items will be named in transmart
colnames(pdSubset) <- c(  "individual","allergic","time", "sample_id", "tissue", "sex", "age",  "organism"   )

#individual, source and treatment are all in the source_name_ch1 column. Split them to seperate columns
#get the individual information from the title
pdSubset <- splitLevelsOfDataFrameColumsTail(pdSubset, c( "individual"), "_")

#get the treatment information from the title
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c("allergic"), "_", 3)
levels(pdSubset$allergic) <- c("control", "nickel")

#get the time information from the title
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c("time"), "_", 2)
levels(pdSubset$time) <- c("00h", "00h", "48h", "07h", "96h")

#get the tissue type from the source_name_ch1
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c("tissue"), " ", 5)
levels(pdSubset$tissue) <- c("upper nates")

levels(pdSubset$sex)
levels(pdSubset$sex) <- c("Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female", "Female")

levels(pdSubset$age)
levels(pdSubset$age) <- c("33-49","33-49","33-49","33-49","33-49","33-49","33-49","33-49","33-49","33-49","33-49","33-49", "33-49")


#select a further subset of the phenotypic data that describe the individuals. One column should uniquely identify the individuals
pdClinicalSubset = pdSubset[c( "individual", "allergic", "tissue", "sex", "age", "organism" )]

#collapse the clinicalTable to individual (remove the timepoint and lesional dimension)
pdClinicalSubset <- unique(pdClinicalSubset)
#set the individual ids as rowNames
row.names(pdClinicalSubset) <- pdClinicalSubset$individual

#set category codes for every phenotypic data as a named list. Under these category codes the phenotypic data will be listed in transmart
categoryCodes <- list()
categoryCodes$individual     <- "Subjects"
categoryCodes$allergic       <- "Subjects"
categoryCodes$organism       <- "Subjects"
categoryCodes$tissue         <- "Subjects"
categoryCodes$sex            <- "Subjects"
categoryCodes$age            <- "Subjects"


#select the individual sample mapping from the first subset
#add the platform_id and tissueType
#add time and lesional columns as the  attributes attributes
pdExpressionSubset = pdSubset[c( "individual", "sample_id" )]
pdExpressionSubset["PLATFORM"] <- attr(gse, "annotation")
pdExpressionSubset["TISSUETYPE"] <- as.character(pdSubset$tissue)
pdExpressionSubset["ATTR1"] <- pdSubset$time
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






























