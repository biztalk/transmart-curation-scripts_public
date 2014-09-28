#for  GEOquery installation instruction and documentatoin see the links below
#http://www.bioconductor.org/packages/release/bioc/html/GEOquery.html
#http://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.pdf

#load the GEOQUERY and R.utils library
library(GEOquery)
library("R.utils")
source("/home/wim/Projects/GEOQuery/TransmartGEOAdapter.R")


#set the GEO study to load and the base output directory
geo_id <- "GSE27887"
outputDir = "/home/wim/Projects/compare_GEO_R_Scripts/bioConductor/"

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
pdSubset = pd[c( "geo_accession", "source_name_ch1",  "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2", "characteristics_ch1.3","characteristics_ch1.4"   )]


#set the labels to use for the phenotypic data. This is how the items will be named in transmart
colnames(pdSubset) <- c( "sample_id", "condition", "organism", "treatment", "individual", "time", "tissue","scorad"   )

#remove the keys from the phenotypic columns by splitting the levels for the column on ": " and selecting index 2
#ie "time: day 1, time: day 2" becomes "day1, day 2"
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c( "treatment", "individual", "time", "tissue","scorad"), ": ", 2)
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c( "condition"), ", ", 3)

#select a further subset of the phenotypic data that describe the individuals. One column should uniquely identify the individuals
pdClinicalSubset = pdSubset[c( "individual", "organism","treatment"  )]
pdClinicalSubset["scorad_pre_treatment"] <- NA
pdClinicalSubset["scorad_post_treatment"] <- NA

pdClinicalSubset <- unique(pdClinicalSubset)
#set the individual ids as rowNames
row.names(pdClinicalSubset) <- pdClinicalSubset$individual



# get the individual, timepoint and scorad score from the first subset
# use it to set the scorad scores for the different time points in the clinicalSubset
pdScoradSubset = pdSubset[c( "individual", "time","tissue", "scorad" )]

#change the levels of time in the scoradSubset so we can use them to acces the right scorad column in the clinicalSubset
levels(pdScoradSubset$time)
levels(pdScoradSubset$time) <- c("scorad_post_treatment", "scorad_pre_treatment")


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
        tissue <- pdScoradSubset[[i, 3]]
        scorad <- pdScoradSubset[[i, 4]]

        if(tissue == "lesional skin"){

            pdClinicalSubset[individual, timepoint] <- scorad
        }
    }
}

# not all scorad scores wre set under lesional skins, sometimes under non lesional skin.
# Copy these manual to the clinicalTable
# P9 does not seem to have post treatment samples records
pdClinicalSubset["P5","scorad_pre_treatment" ] <-pdScoradSubset["GSM688769", "scorad"]

#set category codes for every phenotypic data as a named list. Under these category codes the phenotypic data will be listed in transmart
categoryCodes <- list()
categoryCodes$organism       <- "Subjects"
categoryCodes$treatment      <- "Subjects"
categoryCodes$individual     <- "Subjects"
categoryCodes$scorad_pre_treatment         <- "Subjects"
categoryCodes$scorad_post_treatment        <- "Subjects"


#select the individual sample mapping from the first subset
#add the platform_id and tissueType
#add time and lesional columns as the  attributes attributes
pdExpressionSubset = pdSubset[c( "individual", "sample_id" )]
pdExpressionSubset["PLATFORM"] <- attr(gse, "annotation")
pdExpressionSubset["TISSUETYPE"] <- "skin"
pdExpressionSubset["ATTR1"] <- pdSubset$time
pdExpressionSubset["ATTR2"] <- pdSubset$tissue




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































