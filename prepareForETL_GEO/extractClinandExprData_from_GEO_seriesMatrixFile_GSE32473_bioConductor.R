#for  GEOquery installation instruction and documentatoin see the links below
#http://www.bioconductor.org/packages/release/bioc/html/GEOquery.html
#http://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.pdf

#load the GEOQUERY and R.utils library
library(GEOquery)
library("R.utils")
source("/home/wim/Projects/GEOQuery/TransmartGEOAdapter.R")


#set the GEO study to load and the base output directory
geo_id <- "GSE32473"
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
pdSubset = pd[c( "geo_accession", "source_name_ch1", "source_name_ch1", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2"    )]

#set the labels to use for the phenotypic data. This is how the items will be named in transmart
colnames(pdSubset) <- c(  "sample_id", "individual", "source", "treatment", "organism", "condition", "tissue", "peasi_reduction_perc"   )

#individual, source and treatment are all in the source_name_ch1 column. Split them to seperate columns
#get the patient information from source_name_ch1
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c( "individual"), ", ", 1)
#get the source
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c( "source"), " ", 3)
#make the labels more clear, ad arm to lef and right
levels(pdSubset$source)
levels(pdSubset$source) <- c("baseline", "left arm", "right arm")

pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c( "treatment"), " ", 7)


#change the levels for a certain type of data
#ie gender(f, m) to gender(female, male)

#not all samples have the same characteristics. This results in some samples having wrong values under a column
#These lines fix this manually
pdSubset$peasi_reduction_perc[6] <- "peasi reduction [%]: 0"
pdSubset$peasi_reduction_perc[7] <- "peasi reduction [%]: 0"
pdSubset$peasi_reduction_perc[8] <- "peasi reduction [%]: 0"


#remove the keys from the phenotypic columns by splitting the levels for the column on ": " and selecting index 2
#ie "time: day 1, time: day 2" becomes "day1, day 2"
pdSubset <- splitLevelsOfDataFrameColums(pdSubset, c("condition", "tissue", "peasi_reduction_perc" ), ": ", 2)


#select a further subset of the phenotypic data that describe the individuals. One column should uniquely identify the individuals
pdClinicalSubset = pdSubset[c( "individual", "organism", "condition", "tissue" )]

#add empty columns to store the scorad scores for the different time points
pdClinicalSubset["betamethasone_peasi_reduction_perc"] <- NA
pdClinicalSubset["pimecrolimus_peasi_reduction_perc"] <- NA

#collapse the clinicalTable to individual (remove the timepoint and lesional dimension)
pdClinicalSubset <- unique(pdClinicalSubset)
#set the individual ids as rowNames
row.names(pdClinicalSubset) <- pdClinicalSubset$individual


# get the individual, timepoint and scorad score from the first subset
# use it to set the scorad scores for the different time points in the clinicalSubset
pdPeasi_reductionSubset = pdSubset[c( "individual", "treatment", "peasi_reduction_perc" )]

#change the levels of time in the scoradSubset so we can use them to acces the right scorad column in the clinicalSubset
levels(pdPeasi_reductionSubset$treatment)
levels(pdPeasi_reductionSubset$treatment) <- c("betamethasone_peasi_reduction_perc",  "pimecrolimus_peasi_reduction_perc" )

#convert the dataframe with factors to a matrix
pdPeasi_reductionSubset <- as.matrix(pdPeasi_reductionSubset)

#for the unique individuals in the scoradSubset
for (ind in unique(pdPeasi_reductionSubset[,1])) {

    #get the timepoints records availble for the individual
    individualTimePointRecordsIndexes <- which(pdPeasi_reductionSubset[,1] == ind)

    #get the individual, timepoint and scorad score from the record and set the socrad timepoint combination in the clinical subset
    for(i in individualTimePointRecordsIndexes)
    {
        individual <- pdPeasi_reductionSubset[[i, 1]]
        treatment <- pdPeasi_reductionSubset[[i, 2]]
        peasi_reduction_perc <- pdPeasi_reductionSubset[[i, 3]]

        if(!is.na(treatment) ){
            pdClinicalSubset[individual, treatment] <- peasi_reduction_perc
        }
    }
}


#set category codes for every phenotypic data as a named list. Under these category codes the phenotypic data will be listed in transmart
categoryCodes <- list()
categoryCodes$individual     <- "Subjects"
categoryCodes$condition      <- "Subjects"
categoryCodes$organism       <- "Subjects"
categoryCodes$condition      <- "Subjects"
categoryCodes$tissue         <- "Subjects"
categoryCodes$betamethasone_peasi_reduction_perc   <- "Subjects"
categoryCodes$pimecrolimus_peasi_reduction_perc    <- "Subjects"

#overwrite NA in treatment column with none
pdSubset$treatment <- factor(pdSubset$treatment, exclude = NULL )
levels(pdSubset$treatment) <- c("betamethasone","pimecrolimus",  "none")

#select the individual sample mapping from the first subset
#add the platform_id and tissueType
#add time and lesional columns as the  attributes attributes
pdExpressionSubset = pdSubset[c( "individual", "sample_id" )]
pdExpressionSubset["PLATFORM"] <- attr(gse, "annotation")
pdExpressionSubset["TISSUETYPE"] <- as.character(pdSubset$tissue)
pdExpressionSubset["ATTR1"] <- pdSubset$source
pdExpressionSubset["ATTR2"] <- pdSubset$treatment




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






























