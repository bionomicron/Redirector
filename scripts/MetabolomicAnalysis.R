#-----------------------
# Analysis of correlation between genomic
# and metabolic data for engineered strains
#-----------------------

MetaboliteAnalysis <- function(varData,coverData,metData,metaboliteCols,minCoverage=100,
                        seqIDCol="Sequencing.ID",varID="Row Names"){
                            
print("Starting Metabolic Analysis Version 1.0")

#replace row names with data from variant ID column
vColID <- dimnames(varData)[[2] == varID
varDataRowIDs <- as.vector(varData[,vColID])
varData <- varData[,!vColID]
dimnames(varData)[[2]] <- varDataRowIDs

cColID = dimnames(coverData)[[2] == varID
coverDataIDs <- as.vector(coverData[,cColID])
coverData <- coverData[,!cColID]
dimnames(coverData)[[2]] <- coverDataIDs

varBinary <- varData != "N/A"
coverMask <- coverData > minCoverData

varRowNames = dimnames(varData)[[1]]
varColumnNames = dimnames(varData)[[2]]

#Set strain IDs as row names for metabolic data
dimnames(metData)[[1]] = metData[,seqIDCol]

result = c()
for (metCol in metaboliteCols){
	stats = c()
	print(c("Analysis of metabolite",metCol))
	for (r in 1:length(varRowNames)){
		print(c("Row",r))
        
		activeVar = varBinary[r,]
        coveredVar = coverMask[r,]
        
        print(c("varBinary",activeVar))
        print(c("coverMask",coveredVar))
        
        activeVar[!coverMask[]] = NA
        
		print(c("activeVar",activeVar))
        
        n1 = varColumnNames[activeVar==T ]
		n2 = varColumnNames[activeVar==F]
		v1 = as.numeric(metData[n1,metCol])
		v1 = v1[!is.na(v1)]
		v2 = as.numeric(metData[n2,metCol])
		v2 = v2[!is.na(v2)]
        
        print(c("v1",v1))
		print(c("v2",v2))
		
        if (length(v1) < 3 || length(v2) < 3){
			stats = append(stats,NA)
		}
		else{
			z = t.test(v1,v2)
			#s = as.numeric(z[[1]])
			s = z[["p.value"]]
			stats = append(stats,s)
		}
	}
	result = cbind(result,stats)
    
	return(result)
}

dimnames(result)[[2]] = metaboliteCols
return (result)	

}
