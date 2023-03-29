#Xander Nuttle
#pdf_pb_mips.r
#Call: Rscript /data/talkowski/xander/MIPs/analysis_programs/pdf_pb_mips.r experiment_name GENE_NAME
#
#This R script should be invoked from a shell script
#Generates raw allele/paralog-specific count frequency plots for regions targeted by MIPs in the PB MIP pool

args<-commandArgs(TRUE)
base_name<-args[1]
gene<-args[2]
experiment_name<-paste(base_name,"_",gene,".mipcounts",sep="")
barcodekey_file<-paste(base_name,".barcodekey",sep="")
barcodekey<-read.table(barcodekey_file,header=F,sep='\t',colClasses="character")
guide_locations<-paste(gene,".guidelocs",sep="")

pdfname<-paste(base_name,"_",gene,".pdf",sep="")
pdf(pdfname,width=11,height=7,useDingbats=FALSE)
plot_program<-"/data/talkowski/xander/MIPs/analysis_programs/mipplot_haps.r"
for(indiv in 1:length(barcodekey[,1]))
{
	individual<-barcodekey[indiv,1]
	source(plot_program)
	mytitle<-individual
	p<-p+ggtitle(mytitle)
	print(p)
}
dev.off()

