library(GEOquery)
library("R.utils")
source("/home/wim/Projects/GEOQuery/TransmartGEOAdapter.R")


#set the GEO study to load and the base output directory
geo_id <- "GSE58558"
#geo_id <- "GSE45016"
outputDir = "/home/wim//Projects/compare_GEO_R_Scripts/bioConductor/"

#create the output subdirectories and the outputFilePaths
outputInfo <- createOutputDirsAndFiles(outputDir, geo_id)

#load the GSE
gseList <- getGEO(geo_id)

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

#remove the sample information from the title
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c("individual"), " ", 1)

#remove the keys from the phenotypic columns by splitting the levels for the column on ": " and selecting index 2
#ie "time: day 1, time: day 2" becomes "day1, day 2"
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c("time","lesional",  "responder", "age", "gender", "scorad" ), ": ", 2)

levels(pdSubset$time) <- c("week_00","week_12", "week_02")

#select a further subset of the phenotypic data that describe the individuals. One column should uniquely identify the individuals
pdClinicalSubset = pdSubset[c( "individual", "organism", "responder", "age", "gender" )]

#add empty columns to store the scorad scores for the different time points
pdClinicalSubset["scorad_week00"] <- NA
pdClinicalSubset["scorad_week02"] <- NA
pdClinicalSubset["scorad_week12"] <- NA

#collapse the clinicalTable to individual (remove the timepoint and lesional dimension)
pdClinicalSubset <- unique(pdClinicalSubset)
#set the individual ids as rowNames
row.names(pdClinicalSubset) <- pdClinicalSubset$individual


# get the individual, timepoint and scorad score from the first subset
# use it to set the scorad scores for the different time points in the clinicalSubset
pdScoradSubset = pdSubset[c( "individual", "time", "scorad" )]
#collapse to timepoints (remove the lesional dimension)
pdScoradSubset <- unique(pdScoradSubset)

#change the levels of time in the scoradSubset so we can use them to acces the right scorad column in the clinicalSubset
levels(pdScoradSubset$time)
levels(pdScoradSubset$time) <- c("scorad_week00",  "scorad_week12" , "scorad_week02")

#convert the dataframe with factors to a matrix
pdScoradSubset <- as.matrix(pdScoradSubset)

#for the unique individuals in the scoradSubset
for (ind in unique(pdScoradSubset[,1])) {

    #get the timepoints records availble for the individual
    individualTimePointRecordsIndexes <- which(pdScoradSubset[,1] == ind)

    #get the individual, timepoint and scorad score from the record and set the socrad timepoint combination in the clinical subset
    for(i in individualTimePointRecordsIndexes)
    {
        individual <- pdScoradSubset[[i, 1]]
        timepoint <- pdScoradSubset[[i, 2]]
        scorad <- pdScoradSubset[[i, 3]]
        if(scorad == "NA"){scorad <- ""}    #replace NA values with empty string

        pdClinicalSubset[individual, timepoint] <- scorad
    }
}

#set category codes. Under these category codes the phenotypic data will be listed in transmart
categoryCodes <- list()
categoryCodes$organism       <- "Subjects"
categoryCodes$responder      <- "Subjects"
categoryCodes$age            <- "Subjects"
categoryCodes$gender         <- "Subjects"
categoryCodes$scorad_week00         <- "Subjects"
categoryCodes$scorad_week02         <- "Subjects"
categoryCodes$scorad_week12         <- "Subjects"


#select the individual sample mapping from the first subset
#add the platform_id and tissueType
#add time and lesional columns as the  attributes attributes
pdExpressionSubset = pdSubset[c( "individual", "sample_id" )]
pdExpressionSubset["PLATFORM"] <- attr(gse, "annotation")
pdExpressionSubset["TISSUETYPE"] <- as.character(pd$source_name_ch1)
pdExpressionSubset["ATTR1"] <- pdSubset$time
pdExpressionSubset["ATTR2"] <- pdSubset$lesional



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





























