#' Generator of baseline abundances based on GSE227750 data
#'
#' This object is a function returned from `stats::approxfun` which is used
#' in the simulation study to generate baseline abundances of gene-level expression.
#' The gene list used to create `baselineAbundance_function` is exported as a character
#' vector in the `baselineAbundance_genes` object.
#'
#' @docType data
#'
#' @usage 
#' data(baselineAbundance_function) 
#' data(baselineAbundance_genes)
#' 
#' @details 
#' Genes included in `baselineAbundance_genes` are protein-coding and lncRNA 
#' genes from the reference chromosomes 1, ..., 19, X, and Y from the mouse 
#' Ensembl version 112 annotation (GENCODE - Mus musculus - release M35). Only
#' genes with expected CPM>1 in at least 6 out of the 9 libraries are used.
#'
#' @format 
#' baselineAbundance_function is a function and baselineAbundance_genes is a
#' character vector
#'
#' @keywords datasets
#'
#' @references Milevskiy MJG, Coughlan HD, Kane SR, Johanson TM et al. Three-dimensional genome architecture coordinates key regulators of lineage specification in mammary epithelial cells. Cell Genom 2023 Nov 8;3(11):100424. PMID: 38020976
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227750}{NCBI}
#'
#' @examples
#' \dontrun{
#' library(AnnotationHub) #. v3.12.0
#' library(tximport) #. v1.32.0
#' library(SummarizedExperiment) #. v1.33.3
#' library(edgeR) #. v4.2.0
#' library(data.table) #. v1.15.4
#' library(Rsubread) #. v2.18.0
#' library(rtracklayer) #. v1.62.0
#' 
#' #. Loading annotation
#' gtf <- import('../../data/annotation/mm39/gencode.vM35.annotation.gtf.gz')
#' gtf.tx <- gtf[gtf$type == 'transcript']
#' tx2gene <- data.frame(tx = gtf.tx$transcript_id,gene = gtf.tx$gene_id)
#' 
#' #. Loading quantifications
#' files.quant <- list.files('../../output/mouse/salmon','quant.sf',recursive = TRUE,full.names = TRUE)
#' names(files.quant) <- basename(dirname(files.quant))
#' se.gene <- tximport(files = files.quant,type = 'salmon',tx2gene = tx2gene,countsFromAbundance = 'no',txOut = FALSE,ignoreAfterBar = TRUE)
#' 
#' #. Loading fasta information to check duplication
#' df.fa <- as.data.table(scanFasta('../../data/annotation/mm39/gencode.vM35.transcripts.fa.gz'))
#' TranscriptID.fa <- strsplit2(df.fa$TranscriptID,"\\|")[,1]
#' GeneID.fa <- strsplit2(df.fa$TranscriptID,"\\|")[,2]
#' TranscriptBioType.fa <- strsplit2(df.fa$TranscriptID,"\\|")[,8]
#' df.fa$TranscriptID <- TranscriptID.fa
#' df.fa$GeneID <- GeneID.fa
#' df.fa$TranscriptBioType <- TranscriptBioType.fa
#' df.fa[,allUnique := all(Unique),by = 'GeneID']
#' df.fa <- df.fa[order(GeneID,TranscriptID),]
#' 
#' #. Filtering transcripts from 
#' #. (1) protein coding and lncRNA genes, 
#' #. (2) from standard chromosomes, and 
#' #. (3) from genes with only non-duplicated transcripts
#' is.tx.relev <- gtf.tx$gene_type %in% c('protein_coding','lncRNA') & #. (1)
#'   as.character(seqnames(gtf.tx)) %in% paste0('chr',c(1:19,'X','Y')) & #. (2)
#'   gtf.tx$gene_id %in% df.fa[allUnique == TRUE,unique(GeneID)] #. (3)
#' gene.relev <- unique(gtf.tx$gene_id[is.tx.relev])
#' is.gene.relev <- rownames(se.gene$counts) %in% gene.relev
#' se.gene.subset <- se.gene$counts[is.gene.relev,]
#' 
#' #. Estimating baseline proportions (goodTuringProportions needs integers)
#' cts.gene <- round(se.gene.subset)
#' prop <- goodTuringProportions(cts.gene)
#' 
#' #. Choose genes with expected CPM>1 in at least 6 libraries
#' is.expr <- rowSums(prop > 1e-6) >= 6
#' prop <- prop[is.expr,]
#' baselineAbundance_genes <- rownames(prop)
#' 
#' #. Creating interpolation function and saving object
#' i <- seq(from=0,to=1,length=length(prop))
#' baselineAbundance_function <- approxfun(i,sort(prop),rule=2)
#' 
#' save(baselineAbundance_genes,file = 'data/baselineAbundance_genes.rda',compress = 'xz')
#' save(baselineAbundance_function,file = 'data/baselineAbundance_function.rda',compress = 'xz')
#'  }
"baselineAbundance_function"
"baselineAbundance_genes"