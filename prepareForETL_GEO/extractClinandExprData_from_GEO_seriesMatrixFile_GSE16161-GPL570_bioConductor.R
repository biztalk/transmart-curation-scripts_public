#for  GEOquery installation instruction and documentatoin see the links below
#http://www.bioconductor.org/packages/release/bioc/html/GEOquery.html
#http://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.pdf

#load the GEOQUERY and R.utils library
library(GEOquery)
library("R.utils")
source("/home/wim/Projects/GEOQuery/TransmartGEOAdapter.R")


#set the GEO study to load and the base output directory
geo_id <- "GSE16161"
outputDir = "/home/wim//Projects/compare_GEO_R_Scripts/bioConductor/"

#load the GSE
gseList <- getGEO(geo_id)

#concate the platform_id to the geo_id because this geo study uses two different platforms
geo_id <-paste(geo_id, "-GPL570",sep="")

#create the output subdirectories and the outputFilePaths
outputInfo <- createOutputDirsAndFiles(outputDir,geo_id )



#get the gse from the list by index or name
gse <- gseList$`GSE16161-GPL570_series_matrix.txt.gz`

# now get the phenotypic data (covariates etc.) using pData()
pd <- pData(gse)

#show the available phenotypic data
names(pd)

#select a subset of the phenotypic data
pdSubset = pd[c( "geo_accession", "source_name_ch1",  "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2"   )]

#set the labels to use for the phenotypic data. This is how the items will be named in transmart
colnames(pdSubset) <- c( "sample_id", "source", "organism", "tissue", "individual", "condition" )


#remove the keys from the phenotypic columns by splitting the levels for the column on ": " and selecting index 2
#ie "time: day 1, time: day 2" becomes "day1, day 2"
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c("tissue", "individual", "condition"), ": ", 2)
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c("source"), ", ", 2)

#change leasional to lesional
levels(pdSubset$source)
levels(pdSubset$source) <- c("lesional", "ma")

pdClinicalSubset = pdSubset[c( "source", "organism", "tissue", "individual", "condition" )]

#set category codes. Under these category codes the phenotypic data will be listed in transmart
categoryCodes <- list()
categoryCodes$organism       <- "Subjects"
categoryCodes$tissue         <- "Subjects"
categoryCodes$condition      <- "Subjects"
categoryCodes$source         <- "Subjects"
categoryCodes$individual     <- "Subjects"

#select the individual to sample mapping from the first subset
#add the platform_id and tissueType
#if required set extra attributes ATTR1 and ATTR2 columns for every sample. Attributes mighs be timePoint or tissue location / type
pdExpressionSubset <- pdSubset[c( "individual", "sample_id" )]
pdExpressionSubset["PLATFORM"] <- attr(gse, "annotation")
pdExpressionSubset["TISSUETYPE"] <- as.character(pdClinicalSubset$tissue)
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





























