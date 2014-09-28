#for  GEOquery installation instruction and documentatoin see the links below
#http://www.bioconductor.org/packages/release/bioc/html/GEOquery.html
#http://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.pdf

#load the GEOQUERY and R.utils library
library(GEOquery)
library("R.utils")
source("/home/wim/Projects/GEOQuery/TransmartGEOAdapter.R")


#set the GEO study to load and the base output directory
geo_id <- "GSE39278"
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
pdSubset = pd[c( "geo_accession","source_name_ch1","source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.4", "characteristics_ch1.5"  )]


#set the labels to use for the phenotypic data. This is how the items will be named in transmart
colnames(pdSubset) <- c( "geo_accession","lesional","condition", "organism", "gender", "age", "breed", "disease_severity", "biopsy_site", "treatment"    )

#get lesional information from the source_name_ch1
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c( "lesional"), " ", 2)

#get condition information from the source_name_ch1
pdSubset <- splitLevelsOfDataFrameColumsTail(pdSubset, c( "condition"), " ")
levels(pdSubset$condition) <- c("atopic dermatitis", "healthy")

#remove the keys from the phenotypic columns by splitting the levels for the column on ": " and selecting index 2
#ie "time: day 1, time: day 2" becomes "day1, day 2"
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c("gender", "age", "breed", "disease_severity", "biopsy_site", "treatment" ), ": ", 2)




#set category codes for every phenotypic data as a named list. Under these category codes the phenotypic data will be listed in transmart
categoryCodes <- list()
categoryCodes$lesional       <- "Measurements"
categoryCodes$condition      <- "Subjects"
categoryCodes$organism       <- "Subjects"
categoryCodes$gender         <- "Subjects"
categoryCodes$age            <- "Subjects"
categoryCodes$breed          <- "Subjects"
categoryCodes$disease_severity  <- "Measurements"
categoryCodes$biopsy_site       <- "Measurements"
categoryCodes$treatment         <- "Measurements"




#write the expressionTable
writeExpressionTable(outputInfo$outputFileExpression, gse)

#construct and write the clinicalTable
writeClinicalTable(outputInfo$outputFileClinical, pdSubset)

#construct and write the clinical table mapping
writeClinicalColumns(outputInfo$outputFileClinicalColumns,outputInfo$outputFileClinical,pdSubset, "geo_accession", categoryCodes )

#construct and write the expression table mapping
platform_id <- attr(gse, "annotation")
sampleNames = row.names(pd)
tissueType <- "skin"
categoryCd <- "Biomarker_data+PLATFORM+TISSUETYPE"
writeExpressionMapping(outputInfo$outputFileExpressionMapping, sampleNames, geo_id, platform_id, tissueType, categoryCd )

#construct and write clinicalParams file
writeClinicalParams(outputInfo$outputFileClinicalParams, outputInfo$outputFileClinicalColumns )

#construct and write expressionParams file
writeExpressionParams(outputInfo$outputFileExpressionParams, geo_id, outputInfo$outputFileExpressionMapping)





























