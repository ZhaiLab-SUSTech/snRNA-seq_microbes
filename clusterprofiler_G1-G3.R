library(clusterProfiler)
library(optparse)
option_list <- list(
    make_option(c('-b','--database'),type = 'character',default = FALSE,action='store',help='GSEA database'),
    make_option(c('-f','--DEGfile'),type = 'character',default = FALSE,action='store',help='diffxpy result'),
    make_option(c('-o','--outdir'),type = 'character',default = FALSE,action='store',help='output path'),
    make_option(c('-c','--cluster'),type = 'character',default = FALSE,action='store',help='cluster'),
    make_option(c('-d','--direction'),type = 'character',default = FALSE,action='store',help='direction'))

opt = parse_args(OptionParser(option_list = option_list, usage = 'GSEA analysis from diffxpy result'))
GSEA_analysis <- function(database,DEGfile,outdir,clustername,directionname){
    plantGSEA <- clusterProfiler::read.gmt(database)
    DEG_file_df <- read.csv(DEGfile)
    DEG_file_df_cluster_direction <- subset(DEG_file_df,cluster== clustername & direction == directionname) 
    enrich_result = enricher(gene=DEG_file_df_cluster_direction$gene,TERM2GENE=plantGSEA,)
    split_path = as.vector(strsplit(DEGfile,split='/')[[1]])
    split_name = as.vector(strsplit(tail(split_path,n=1),split='_')[[1]])
    sample = head(split_name,n=1)
    #save table
    file_title = paste(sample,clustername,directionname,sep='_')
    genelistname = paste(paste(outdir,"/",file_title,sep=''),'genelist.csv',sep='_')
    write.csv(DEG_file_df_cluster_direction,genelistname)
    filename = paste(paste(outdir,"/",file_title,sep=''),'clusterprofiler_result.csv',sep='_')
    write.csv(enrich_result,filename)
    #plot
    
    plot_title = paste(sample,clustername,directionname,sep='_')
    barplot_filename = paste(paste(outdir,"/",plot_title,sep=''),'barplot.png',sep='_')
    dotplot_filename = paste(paste(outdir,"/",plot_title,sep=''),'dotplot.png',sep='_')

    #barplot
    png(file = barplot_filename)
    print(barplot(enrich_result,
        showCategory = 20,
        font.size = 8,
        title = plot_title))
    dev.off()

    #dotplot
    png(file = dotplot_filename)
    print(dotplot(enrich_result,
        showCategory = 20,
        font.size = 8,	
        title = plot_title))
    dev.off()
    return(paste(filename,barplot_filename,dotplot_filename,'finish!',sep=','))
}

tryc <- tryCatch({
result <- GSEA_analysis(opt$database,opt$DEGfile,opt$outdir,opt$cluster,opt$direction
)},error = function(e){print(paste(opt$DEGfile,'Error:',e,sep=','));return(NA)},
   warning = function(w){print(paste(opt$DEGfile,'Warning:',w,sep=','));return(NA)},
   finally = {print(paste(opt$DEGfile,'Finish!',sep=','))}
)

