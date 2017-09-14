#This function calls the STAMS program
#Need a few things to be specified:
# -where the VEGAS *.out files are
# -what the phenotype is (for file naming)
# -permutation designations for search_results and normalized_results
STAMS <- function(outfiledir = "STAMS/META_PLMD",
                  pheno = "PLMD",
                  SRperm = 10000,
                  NRperm = 1000000){
    library(STAMS)
    string_db_instance <- STRINGdb$new()
    dirpattern <- "STAMS.R"
    myfunction.path <- list.files(getwd(),recursive=TRUE,full.names=TRUE,pattern=dirpattern)
    setwd(paste(gsub(dirpattern,'\\1',myfunction.path),outfiledir,sep="/"))
    outfilelist <- list.files(pattern='.out$')
    
    # LOAD PVALUES
    string_db <- STRINGdb$new()
    gene_values=read.table(file=outfilelist, header=T, sep=' ')
    # get rid of pvalue = 0 lines and replace with upper bound
    gene_values$Pvalue=pmax(gene_values$Pvalue, 1/gene_values$nSims)
    
    # remove pvalue = 1
    print("Removing genes with a pvalue = 1")
    print(paste("Total genes removed: ",length(which(gene_values$Pvalue ==1)),sep=""))
    gene_values=gene_values[-which(gene_values$Pvalue == 1),]
    
    gene_mapped = string_db$map(gene_values, "Gene", removeUnmappedRows = T)
    # check mapping
    data3=string_db$add_proteins_description(gene_mapped)
    data4=data3[which(data3$Gene==toupper(data3$preferred_name)),]
    ensg_genes=data3[grep('ENSG',data3$preferred_name),]
    data4=rbind(data4,ensg_genes)
    mapped_data=data4[,-which(names(data4)=='annotation')]
    
    # TO SAVE THIS WORK SO THAT YOU CAN LOAD IT UP NEXT TIME
    save.image(paste(pheno,"_mapped_data.Rdata",sep=''))
    
    mapped_data$weight=qnorm(1-mapped_data$Pvalue)
    network_with_weights=get_string_edges(mapped_data, "combined_score", string_db_instance)
    
    out.path <- paste(gsub(dirpattern,'\\1',myfunction.path),"STAMS/out-STAMS",sep="/")
    setwd(out.path)
    
    t1 <- Sys.time()
    print(paste("STAMS search starting at ",t1,sep = ""))
    # For a real experiment with default lambda, set lambda="default" and lambda_perm to 10,000 +
    search_results=stams_search(network_with_weights,
                                r=0.1,
                                lambda="default",
                                lambda_perm=SRperm)
    save.image(paste(pheno,"_STAMSsearch.Rdata",sep = ""))
    t2 <- Sys.time()
    
    t3 <- Sys.time()
    print(paste("Score normalization starting at ",t3,sep = ""))
    # For a real experimental, increase the number of permutations here to 100,000 +
    normalized_results=normalize_scores(network_with_weights,
                                        search_results,
                                        filename=paste("STAMS_",pheno,".txt",sep = ""),
                                        perm=NRperm)
    t4 <- Sys.time()
    
    t5 <- Sys.time()
    print(paste("Choose & Plot starting at ",t5,sep = ""))
    # This method is from dmGWAS, the plot uses tcl/tk. To plot the results this way, set plot=TRUE
    chosen=chooseModule(normalized_results,
                        top=0.01,
                        plot=FALSE)
    # This plot doesn't need tcl/tk or X
    png=string_db_instance$get_png(chosen$modules[[1]],
                                   file=paste("STAMS_",pheno,".png",sep = ""))
    # This plot uses X
    string_db_instance$plot_network(chosen$modules[[1]],
                                    add_summary=FALSE)
    t6 <- Sys.time()
    
    clock <- data.frame(Start_Stop=as.POSIXct(c(t1, t2, t3, t4, t5, t6)),
                        Process=as.character(c("STAMS search start",
                                               "STAMS search stop",
                                               "Normalization start",
                                               "Normalization stop",
                                               "Choose & Plot start",
                                               "Choose & Plot stop")))
    write.table(clock, file = "Process_speed.txt",quote = F,row.names = F)
}

STAMS(outfiledir = "STAMS/META_PLMD", pheno = "PLMD", SRperm = 10000, NRperm = 1000000)