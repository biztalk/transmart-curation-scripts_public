#for  GEOquery installation instruction and documentatoin see the links below
#http://www.bioconductor.org/packages/release/bioc/html/GEOquery.html
#http://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.pdf

#load the GEOQUERY and R.utils library
library(GEOquery)
library("R.utils")
source("/home/wim/Projects/GEOQuery/TransmartGEOAdapter.R")


#set the GEO study to load and the base output directory
geo_id <- "GSE36842"
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
pdSubset = pd[c( "geo_accession", "source_name_ch1",  "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2", "characteristics_ch1.3","characteristics_ch1.4", "characteristics_ch1.5", "characteristics_ch1.6", "characteristics_ch1.7"   )]

#set the labels to use for the phenotypic data. This is how the items will be named in transmart
colnames(pdSubset) <- c( "sample_id","source",  "organism", "individual", "condition", "scorad" ,"age", "family_history_of_atopy", "felevated_ige","eosinophil_count", "filaggrin_mutation" )

#remove the keys from the phenotypic columns by splitting the levels for the column on ": " and selecting index 2
#ie "time: day 1, time: day 2" becomes "day1, day 2"
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c("individual", "condition", "scorad", "age", "family_history_of_atopy", "felevated_ige", "eosinophil_count", "filaggrin_mutation"), ": ", 2)

#normal samples don't have all characteristics_ch(4,5,6,7,) . Somehow this means that unused characteristics_ch  and 9 are now listed under characteristics_ch 4 and 5 for normal samples.
#remove these here
levels(pdSubset$age)
levels(pdSubset$age)[9] <- NA
levels(pdSubset$family_history_of_atopy)
levels(pdSubset$family_history_of_atopy)[1] <- NA
levels(pdSubset$family_history_of_atopy)[1] <- NA       #order of levels changes after overwriting a level

#select a further subset of the phenotypic data that describe the individuals. One column should uniquely identify the individuals
pdClinicalSubset = pdSubset[c( "source", "organism", "individual", "condition", "age", "family_history_of_atopy", "felevated_ige","eosinophil_count", "filaggrin_mutation" )]

#overwrite the NL, CLS and ANL sample information for patient so we can later collapse the clinical table
levels(pdClinicalSubset$condition) <- c("Patient" , "Patient","Patient", "Normal")

#fix one patient that has as age 40, 40 and 41 for 3 samples by overwriting 41 to 40
levels(pdClinicalSubset$age)[5] <- 40

#add the scorad columns for acute and chronic skin. Will be filled later
pdClinicalSubset["acute_scorad"] <- NA
pdClinicalSubset["chronic_scorad"] <- NA

#collapse the clinicalTable to individual level ( ie remove the sample (CLS,NL ANL) dimension )
pdClinicalSubset <- unique(pdClinicalSubset)

#set the individual ids as rowNames
row.names(pdClinicalSubset) <- pdClinicalSubset$individual

# get the individual, condition and scorad score from the first subset
# use it to set the scorad scores for the different time points in the clinicalSubset
pdScoradSubset = pdSubset[c( "individual", "condition", "scorad" )]

#convert the dataframe with factors to a matrix
pdScoradSubset <- as.matrix(pdScoradSubset)

#for the unique individuals in the scoradSubset
for (ind in unique(pdScoradSubset[,1])) {

    #get the timepoints records availble for the individual
    individualRecordsIndexes <- which(pdScoradSubset[,1] == ind)

    #get the individual, timepoint and scorad score from the record and set the socrad timepoint combination in the clinical subset
    for(i in individualRecordsIndexes)
    {
        individual <- pdScoradSubset[[i, 1]]
        condition <- pdScoradSubset[[i, 2]]
        scorad <- pdScoradSubset[[i, 3]]

        if(condition == "ALS"){ pdClinicalSubset[individual, "acute_scorad"] <- scorad }
        if(condition == "CLS"){ pdClinicalSubset[individual, "chronic_scorad"] <- scorad }
    }
}


#set category codes for every phenotypic data as a named list. Under these category codes the phenotypic data will be listed in transmart
categoryCodes <- list()
categoryCodes$source         <- "Subjects"
categoryCodes$organism       <- "Subjects"
categoryCodes$individual     <- "Subjects"
categoryCodes$condition      <- "Subjects"
categoryCodes$acute_scorad         <- "Subjects"
categoryCodes$chronic_scorad       <- "Subjects"
categoryCodes$age            <- "Subjects"
categoryCodes$family_history_of_atopy  <- "Subjects"
categoryCodes$felevated_ige            <- "Subjects"
categoryCodes$eosinophil_count         <- "Subjects"
categoryCodes$filaggrin_mutation       <- "Subjects"

#select the individual sample mapping from the first subset
#add the platform_id and tissueType
#add time and lesional columns as the  attributes attributes
pdExpressionSubset = pdSubset[c( "individual", "sample_id" )]
pdExpressionSubset["PLATFORM"] <- attr(gse, "annotation")
pdExpressionSubset["TISSUETYPE"] <- as.character(pd$source_name_ch1)
pdExpressionSubset["ATTR1"] <- pdSubset$condition
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





























