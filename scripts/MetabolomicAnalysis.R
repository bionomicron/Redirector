#---------------------------------------------
# Analysis of correlation between genomic
# and metabolic data for engineered strains
# Updated 2014.09.01
#---------------------------------------------

#Examples/Users/bionomicron/Projects/Redirector/scripts/MetabolomicAnalysis.R
#var_file_name = 'FC_00678_seq_analysis/Variance_FC_00678_Report.txt'
#coverage_file_name = 'FC_00678_seq_analysis/Coverage_FC_00678_Report.txt'

var_file_name = 'FC_00707_seq_analysis/Variance_FC_00707_Report.txt'
coverage_file_name = 'FC_00707_seq_analysis/Coverage_FC_00707_Report.txt'

var_data = read.table(var_file_name,sep="\t",na.strings=c(''),header=T)
coverage_data = read.table(coverage_file_name,sep="\t",header=T)

met_data = read.table('20140429_FA_GCMS_T40_Data.txt',sep="\t",header=T)
metabolite_cols =
c("OD","Total","Total.OD",
"c08","c10","c12","c12.1","c14","c16","c18",
"c08.OD","c10.OD","c12.OD","c12.1.OD","c14.OD","c16.OD","c18.OD")

#seqIDCol= 'Seq.1.ID'
seq_id= 'Seq.2.ID'

#>result <- MetaboliteAnalysis(varData,coverage,metData,metaboliteCols,minCoverage,seqIDCol="Seq.1.ID",varID="Row.Name")
MetaboliteAnalysis <- function(varData,coverData,metData,metaboliteCols,minCoverage=10,seqIDCol=seq_id,varID="Row_Name"){
                            
print("Starting Metabolic Analysis Version 1.0")
result = list()

#------------------------
# 1.Prepare Data Matrixes
#------------------------

#----Variance Data----
print ("Preparing variance data")
observation_ids = colnames(varData)[varID != colnames(varData)]

#vColID <- dimnames(varData)[[2]] == varID
varDataRowIDs <- as.vector(varData[,varID])
#varDataRowIDs <- as.vector(varData[,varID])
varData <- varData[,observation_ids]

#Order variance data
dimnames(varData)[[1]] <- varDataRowIDs
vOrder <- order(varDataRowIDs)
varData <- varData[vOrder,]
result$var_data = varData

#-----Coverage Data-------
print("Preparing coverage data")
coverDataRowIDs <- as.vector(coverData[,varID])
coverData <- coverData[,observation_ids]

#Order coverage data
coverData <- coverData[vOrder,]

#------Metabolic Data-------
varRowNames = dimnames(varData)[[1]]
varColumnNames = dimnames(varData)[[2]]

#Set strain IDs as row names for metabolic data
print("Preparing metabolic data")
metIDs = as.vector(metData[,seqIDCol])
#print(c("metIDs",seqIDCol,metIDs))
metData <- metData[,metabolite_cols]
dimnames(metData)[[1]] = metIDs

#-----Filter Data-------
print("Filtering data")
data_sum <- apply(metData,1,sum)
filter <- !is.na(data_sum) & data_sum > 0
filter_names <- names(data_sum[filter])
#print(filter_names)
obs_ids <- intersect(filter_names,dimnames(varData)[[2]])
varData <- varData[,obs_ids]
coverData <- coverData[,obs_ids]
metData <- metData[obs_ids,]

#-----Data Mask-----
print("Making data masks")
variant_mask <- varData != "NA"
coverage_mask <- coverData > minCoverage
data_mask <- variant_mask
data_mask[variant_mask] <- 1
data_mask[!variant_mask] <- 0
data_mask[!coverage_mask] <- 0

result$nvar_ids = dimnames(varData)[[2]]
result$obs_ids = dimnames(varData)[[1]]
result$var_data = varData
result$coverage_data = coverData
result$data_mask = t(data_mask)
result$meta_data = metData

return(result)
#---------------------------
# 2. Loop over metabolites
#---------------------------

result = c()
i = 0
for (metCol in metaboliteCols){
	print(c("Analysis of metabolite",metCol))
    variance_names = dimnames(varData)[[1]]
    coverage_names = dimnames(coverData)[[1]]
    
    stats = c()
    
	#for (r in 1:length(varRowNames)){
    for (row_name in variance_names){
        #print(c("Row",row_name))
        #print(c("Row",varRowNames[r]))
        
        #activeVariants= variant_mask[r,]
        #coveredVariants = coverage_mask[r,]
        
        ivariant = as.vector(varData[row_name,])
        icoverage = as.vector(coverData[row_name,])
        activeVariants = variant_mask[row_name,]
        coveredVariants = coverage_mask[row_name,]
        
        activeVar <- activeVariants
        activeVar[!coveredVariants] = NA
        
        stat_data <- t(rbind(icoverage,ivariant))
        x <- cbind(activeVariants,coveredVariants,activeVar)
        
        n1 = varColumnNames[activeVar==T]
		n2 = varColumnNames[activeVar==F]
		v1 = as.numeric(metData[n1,metCol])
		v1 = v1[!is.na(v1)]
		v2 = as.numeric(metData[n2,metCol])
		v2 = v2[!is.na(v2)]
        
        #print(stat_data)
        #print(x)
        #print(c("v1",v1))
		#print(c("v2",v2))
        
        if (length(v1) < 2 || length(v2) < 2){
			stats = append(stats,NA)
		}
		else{
			z = t.test(v1,v2)
			#s = as.numeric(z[[1]])
			s = z[["p.value"]]
			stats = append(stats,s)
		}
        #print(c("Stats",stats))
        i = i + 1
        #if (i >= 3){ return(stats) }
    }
	result = cbind(result,stats)
}
dimnames(result)[[1]] = variance_names
dimnames(result)[[2]] = metaboliteCols

return (result)	

}
