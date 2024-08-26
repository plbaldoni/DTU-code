#' @importFrom edgeR catchKallisto catchSalmon DGEList
#' @importFrom tximeta tximeta
getCounts <- function(targets,
                      tx.gene,
                      quantifier = c('salmon','kallisto'),
                      count.type = c('raw','scaled','scaledTPM','dtuScaledTPM')){

  if (!quantifier %in% c('salmon', 'kallisto')) {
    stop('quantifier must be either salmon or kallisto')
  }
  if (!count.type %in% c('raw', 'scaled', 'scaledTPM', 'dtuScaledTPM')) {
    stop("count.type must be one of 'raw','scaled','scaledTPM', or 'dtuScaledTPM'")

  }

  if (count.type %in% c('raw', 'scaled')) {
    if (quantifier == 'kallisto') {
      cts <- catchKallisto(paths = targets$path)
    } else{
      cts <- catchSalmon(paths = targets$path)
    }
    anno <- strsplit2(rownames(cts$counts),"\\|")
    TranscriptID <- anno[,1]
    cts$annotation$TranscriptID <- TranscriptID
    cts$annotation$GeneID <- tx.gene$GeneID[match(TranscriptID,tx.gene$TranscriptID)]
    rownames(cts$counts) <- rownames(cts$annotation) <- TranscriptID
    if (count.type == 'scaled'){
      cts <- DGEList(counts = cts$counts/cts$annotation$Overdispersion,
                     samples = targets,
                     group = targets$group,
                     genes = cts$annotation)
    } else{
      cts <- DGEList(counts = cts$counts,
                     samples = targets,
                     group = targets$group,
                     genes = cts$annotation)
    }
    colnames(cts) <- basename(colnames(cts))
  } else{
    if (quantifier == 'kallisto') {
      targets$files <- file.path(targets$path,'abundance.h5')
    } else{
      targets$files <- file.path(targets$path,'quant.sf')
    }
    cts <- tximeta(coldata = targets,type = quantifier,countsFromAbundance = count.type,tx2gene = tx.gene,ignoreAfterBar = TRUE,txOut = TRUE)
  }
  return(cts)
}

filterCounts <- function(x,method,design = NULL){
  OK <- method %in% c('edgeR','limma','satuRn','DRIMSeq','DEXSeq')
  if (!OK) stop('wrong method')

  if (method %in% c('edgeR', 'limma','satuRn','DEXSeq')) {
    # Use default values for filterByExpr
    keep.tx <- filterByExpr(x,min.count = 10,min.total.count = 15)
    y <- x[keep.tx, , keep.lib.sizes = FALSE]

    if (method %in% c('satuRn','DEXSeq')) {
      # For satuRn/DEXSeq we need to explicitly remove the single-transcript genes
      # This is be done automatically in limma/edgeR
      cname <- ifelse(method == 'satuRn','GeneIDNoVersion','GeneID')
      tb.gene <- table(y$genes[[cname]])
      keep.gene <- y$genes[[cname]] %in% names(tb.gene[tb.gene > 1])
      y <- y[keep.gene,,keep.lib.sizes = FALSE]
    }
  }

  if (method == 'DRIMSeq') {
    # The smallest group sample size (the minimum inverse leverage computed from the design matrix)
    n.small <- min(1/diag(design %*% solve(t(design) %*% design) %*% t(design)))
    # To match as close as possible featureByExpr, but min.total.count is not available in dmFilter
    y <- dmFilter(x,min_samps_feature_expr = n.small,min_feature_expr = 10)
  }

  return(y)
}

#' @importFrom edgeR normLibSizes glmQLFit glmQLFTest topTags
#' @importFrom edgeR topSpliceDGE diffSpliceDGE filterByExpr estimateDisp
runEdgeR <- function(targets,quantifier,simes,count.type,tx.gene,legacy){
  dge <- getCounts(targets,tx.gene,quantifier,count.type)
  design <- model.matrix(~group,data = dge$samples)

  dge <- filterCounts(dge,method = 'edgeR')
  dge <- normLibSizes(dge)
  if(legacy){
    dge <- estimateDisp(y = dge,design = design)
  }
  fit <- glmQLFit(dge,design,legacy = legacy)
  ds <- diffSpliceDGE(fit,coef = 2,geneid = 'GeneID',exonid = 'TranscriptID')

  out.transcript <- topSpliceDGE(ds, test = "exon",number = Inf)

  if (simes) {
    out.gene <- topSpliceDGE(ds, test = "Simes",number = Inf)
  } else{
    out.gene <- topSpliceDGE(ds, test = "gene",number = Inf)
  }

  return(list('transcript' = out.transcript,'gene' = out.gene))
}

#' @importFrom edgeR voomLmFit
#' @importFrom limma diffSplice topSplice
runLimma <- function(targets,quantifier,simes,count.type,tx.gene){
  dge <- getCounts(targets,tx.gene,quantifier,count.type)
  design <- model.matrix(~group,data = dge$samples)

  dge <- filterCounts(dge,method = 'limma')
  dge <- normLibSizes(dge)
  fit <- voomLmFit(counts = dge,design = design)
  ds <- diffSplice(fit, geneid = "GeneID", exonid = "TranscriptID")

  out.transcript <- topSplice(ds, test = "t",number = Inf)

  if (simes) {
    out.gene <- topSplice(ds, test = "simes", number = Inf)
  } else{
    out.gene <- topSplice(ds, test = "F", number = Inf)
  }

  return(list('transcript' = out.transcript,'gene' = out.gene))
}

#' @importFrom DRIMSeq dmDSdata dmFilter dmPrecision dmFit dmTest results samples
#' @importFrom SummarizedExperiment assay
runDRIMSeq <- function(targets,quantifier,count.type,tx.gene){
  se <- getCounts(targets,tx.gene,quantifier,count.type)

  if (count.type %in% c('raw', 'scaled')) {
    cts <- se$counts
    anno <- se$genes
  }
  if (count.type %in% c('scaledTPM', 'dtuScaledTPM')) {
    cts <- assay(se,'counts')
    anno <- data.frame(GeneID = as(rowData(se)$gene_id,'character'),
                       TranscriptID = as(rowData(se)$tx_name,'character'))
  }

  counts <- data.frame(gene_id = anno$GeneID,feature_id = anno$TranscriptID,cts,row.names = NULL)

  d <- dmDSdata(counts = counts, samples = targets)
  design <- model.matrix( ~ as.factor(condition), data = samples(d))
  colnames(design) <- c('Intercept','condition')

  d <- filterCounts(d,method = 'DRIMSeq',design = design)
  d <- dmPrecision(d, design = design)
  d <- dmFit(d, design = design)
  d <- dmTest(d, coef = "condition")

  out.gene <- DRIMSeq::results(d)
  colnames(out.gene) <- c("GeneID","lr","df","PValue","FDR")

  out.transcript <- DRIMSeq::results(d,level = "feature")
  colnames(out.transcript) <- c("GeneID",'TranscriptID',"lr","df","PValue","FDR")

  # For NA p-values, set 1
  # DRIMSeq will not estimate a dispersion parameters for features with 0 counts in one group
  # See https://f1000research.com/articles/7-952/v3
  out.gene$PValue <- out.gene$FDR <- no.na(out.gene$PValue)
  out.transcript$PValue <- out.transcript$FDR <- no.na(out.transcript$PValue)

  return(list('transcript' = out.transcript,'gene' = out.gene))
}

#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData
#' @importFrom limma makeContrasts
#' @importFrom S4Vectors metadata
#' @importFrom satuRn fitDTU testDTU
#' @importFrom BiocParallel SerialParam
runSatuRn <- function(targets,quantifier,count.type,tx.gene){
  se <- getCounts(targets,tx.gene,quantifier,count.type)
  if (count.type %in% c('raw', 'scaled')) {
    TranscriptIDNoVersion <- sub("\\..*", "",se$genes$TranscriptID)
    GeneIDNoVersion <- sub("\\..*", "",se$genes$GeneID)
    cts <- se$counts
    anno <- se$genes
    samp.anno <- se$samples
  }
  if (count.type %in% c('scaledTPM', 'dtuScaledTPM')) {
    TranscriptIDNoVersion <- sub("\\..*", "",as(rowData(se)$tx_name,'character'))
    GeneIDNoVersion <- sub("\\..*", "",as(rowData(se)$gene_id,'character'))
    cts <- assay(se,'counts')
    anno <- data.frame(GeneID = as(rowData(se)$gene_id,'character'),
                       TranscriptID = as(rowData(se)$tx_name,'character'))
    samp.anno <- colData(se) |> as.data.frame()
  }

  dge <- DGEList(counts = cts,samples = samp.anno,genes = anno)
  rownames(dge) <- TranscriptIDNoVersion
  dge$genes$GeneIDNoVersion <- GeneIDNoVersion
  dge$genes$TranscriptIDNoVersion <- TranscriptIDNoVersion

  dge.filtr <- filterCounts(dge,method = 'satuRn')

  # Creating SE object
  sumExp <- SummarizedExperiment(
    assays = list(counts = dge.filtr$counts),
    colData = dge.filtr$samples,
    rowData = data.frame('isoform_id' = rownames(dge.filtr),
                         'gene_id' = dge.filtr$genes$GeneIDNoVersion)
  )

  # Setting up design and contrast
  design <- model.matrix(~0 + group,data = colData(sumExp))
  colnames(design) <- levels(sumExp$group)
  L <- makeContrasts(BvsA = B - A,levels = design)
  sumExp@metadata$formula <- ~ 0 + colData(sumExp)$group

  # Model fitting  and testing (only 1 worker is used)
  BPPARAM <- SerialParam()
  register(BPPARAM)

  sumExp <- fitDTU(object = sumExp,formula = ~ 0 + group,parallel = FALSE,BPPARAM = BPPARAM,verbose = FALSE)
  sumExp <- testDTU(object = sumExp,contrasts = L,diagplot1 = FALSE,diagplot2 = FALSE,sort = FALSE)

  # Getting gene-level adjusted p-values
  pvals <- rowData(sumExp)[["fitDTUResult_BvsA"]]$empirical_pval

  geneID <- factor(rowData(sumExp)$gene_id)
  geneSplit <- split(seq(along = geneID), geneID)
  pGene <- sapply(geneSplit, function(i) min(pvals[i]))
  pGene[is.na(pGene)] <- 1
  theta <- unique(sort(pGene))

  q <- perGeneQValueExact(pGene, theta, geneSplit) # From satuRn vignette
  qScreen <- rep(NA_real_, length(pGene))
  qScreen <- q[match(pGene, theta)]
  qScreen <- pmin(1, qScreen)
  names(qScreen) <- names(geneSplit)

  # Organizing output (satuRn authors recommend 'empirical' p-values)
  out.transcript <-
    data.frame(rowData(sumExp)[,c('gene_id','isoform_id')],
               rowData(sumExp)[["fitDTUResult_BvsA"]][,c('t','df','empirical_pval','empirical_FDR')],
               row.names = NULL)
  colnames(out.transcript) <- c("GeneID",'TranscriptID',"t","df","PValue","FDR")

  # satuRn does not return raw p-values for gene-level tests
  out.gene <- data.frame(GeneID = names(qScreen),PValue = NA,FDR = qScreen)

  # Bringing GeneID and TranscriptID versions
  out.transcript$GeneID <- dge$genes$GeneID[match(out.transcript$GeneID,dge$genes$GeneIDNoVersion)]
  out.transcript$TranscriptID <- dge$genes$TranscriptID[match(out.transcript$TranscriptID,dge$genes$TranscriptIDNoVersion)]

  out.gene$GeneID <- dge$genes$GeneID[match(out.gene$GeneID,dge$genes$GeneIDNoVersion)]
  rownames(out.gene) <- dge$genes$GeneID[match(rownames(out.gene),dge$genes$GeneIDNoVersion)]

  return(list('transcript' = out.transcript,'gene' = out.gene))
}

#' @importFrom DRIMSeq samples
#' @importFrom DEXSeq DEXSeqDataSet estimateSizeFactors estimateDispersions testForDEU perGeneQValue DEXSeqResults
runDEXSeq <- function(targets,quantifier,count.type,tx.gene){
  se <- getCounts(targets,tx.gene,quantifier,count.type)

  if (count.type %in% c('raw', 'scaled')) {
    cts <- se$counts
    anno <- se$genes
    samp.anno <- se$samples
  }
  if (count.type %in% c('scaledTPM', 'dtuScaledTPM')) {
    cts <- assay(se,'counts')
    anno <- data.frame(GeneID = as(rowData(se)$gene_id,'character'),
                       TranscriptID = as(rowData(se)$tx_name,'character'))
    samp.anno <- colData(se) |> as.data.frame()
  }

  dge <- DGEList(counts = cts,samples = samp.anno,genes = anno)

  dge.filtr <- filterCounts(dge,method = 'DEXSeq')

  # DEXSeq pipeline
  cts.filtr <- round(dge.filtr$counts)
  dxd <- DEXSeqDataSet(countData = round(cts.filtr),
                       sampleData = dge$samples,
                       design =  ~ sample + exon + group:exon,
                       featureID = dge.filtr$genes$TranscriptID,
                       groupID = dge.filtr$genes$GeneID)

  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet=TRUE)
  dxd <- testForDEU(dxd, reducedModel=~sample + exon)

  dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)

  out.transcript <- data.frame(GeneID = dxr$groupID,
                               TranscriptID = dxr$featureID,
                               PValue = dxr$pvalue,
                               FDR = dxr$padj,
                               row.names = NULL)

  # DEXSeq does not return raw p-values for gene-level tests
  out.gene <- perGeneQValue(dxr)
  out.gene <- data.frame(GeneID = names(out.gene), PValue = NA, FDR = out.gene,row.names = NULL)

  return(list('transcript' = out.transcript,'gene' = out.gene))
}

callMethods <- function(targets,quantifier,tx.gene){

  res <- list()
  time <- list()

  # edgeR-v4-diffSpliceDGE methods

  time[['edger.v3-scaled-simes']] <-
    system.time({res[['edger.v3-scaled-simes']] <- runEdgeR(targets = targets, quantifier = quantifier,legacy = TRUE,simes = TRUE,count.type = 'scaled', tx.gene = tx.gene)})

  time[['edger.v3-raw-simes']] <-
    system.time({res[['edger.v3-raw-simes']] <- runEdgeR(targets = targets, quantifier = quantifier,legacy = TRUE,simes = TRUE,count.type = 'raw', tx.gene = tx.gene)})

  time[['edger.v3-scaled-ftest']] <-
    system.time({res[['edger.v3-scaled-ftest']] <- runEdgeR(targets = targets, quantifier = quantifier,legacy = TRUE,simes = FALSE,count.type = 'scaled', tx.gene = tx.gene)})

  time[['edger.v3-raw-ftest']] <-
    system.time({res[['edger.v3-raw-ftest']] <- runEdgeR(targets = targets, quantifier = quantifier,legacy = TRUE,simes = FALSE,count.type = 'raw', tx.gene = tx.gene)})

  # edgeR-v4-diffSpliceDGE methods

  time[['edger.v4-scaled-simes']] <-
    system.time({res[['edger.v4-scaled-simes']] <- runEdgeR(targets = targets, quantifier = quantifier,legacy = FALSE,simes = TRUE,count.type = 'scaled', tx.gene = tx.gene)})

  time[['edger.v4-raw-simes']] <-
    system.time({res[['edger.v4-raw-simes']] <- runEdgeR(targets = targets, quantifier = quantifier,legacy = FALSE,simes = TRUE,count.type = 'raw', tx.gene = tx.gene)})

  time[['edger.v4-scaled-ftest']] <-
    system.time({res[['edger.v4-scaled-ftest']] <- runEdgeR(targets = targets, quantifier = quantifier,legacy = FALSE,simes = FALSE,count.type = 'scaled', tx.gene = tx.gene)})

  time[['edger.v4-raw-ftest']] <-
    system.time({res[['edger.v4-raw-ftest']] <- runEdgeR(targets = targets, quantifier = quantifier,legacy = FALSE,simes = FALSE,count.type = 'raw', tx.gene = tx.gene)})

  # limma-diffSplice methods

  time[['limma-scaled-simes']] <-
    system.time({res[['limma-scaled-simes']] <- runLimma(targets = targets, quantifier = quantifier,simes = TRUE,count.type = 'scaled', tx.gene = tx.gene)})

  time[['limma-raw-simes']] <-
    system.time({res[['limma-raw-simes']] <- runLimma(targets = targets, quantifier = quantifier,simes = TRUE,count.type = 'raw', tx.gene = tx.gene)})

  time[['limma-scaled-ftest']] <-
    system.time({res[['limma-scaled-ftest']] <- runLimma(targets = targets, quantifier = quantifier,simes = FALSE,count.type = 'scaled', tx.gene = tx.gene)})

  time[['limma-raw-ftest']] <-
    system.time({res[['limma-raw-ftest']] <- runLimma(targets = targets, quantifier = quantifier,simes = FALSE,count.type = 'raw', tx.gene = tx.gene)})

  # DRIMSeq methods

  time[['drimseq-raw']] <-
    system.time({res[['drimseq-raw']] <- runDRIMSeq(targets = targets, quantifier = quantifier,count.type = 'raw', tx.gene = tx.gene)})

  time[['drimseq-scaledTPM']] <-
    system.time({res[['drimseq-scaledTPM']] <- runDRIMSeq(targets = targets, quantifier = quantifier,count.type = 'scaledTPM', tx.gene = tx.gene)})

  time[['drimseq-dtuScaledTPM']] <-
    system.time({res[['drimseq-dtuScaledTPM']] <- runDRIMSeq(targets = targets, quantifier = quantifier,count.type = 'dtuScaledTPM', tx.gene = tx.gene)})

  # satuRn methods

  time[['saturn-raw']] <-
    system.time({res[['saturn-raw']] <- runSatuRn(targets = targets, quantifier = quantifier,count.type = 'raw', tx.gene = tx.gene)})

  time[['saturn-scaledTPM']] <-
    system.time({res[['saturn-scaledTPM']] <- runSatuRn(targets = targets, quantifier = quantifier,count.type = 'scaledTPM', tx.gene = tx.gene)})

  time[['saturn-dtuScaledTPM']] <-
    system.time({res[['saturn-dtuScaledTPM']] <- runSatuRn(targets = targets, quantifier = quantifier,count.type = 'dtuScaledTPM', tx.gene = tx.gene)})

  # DEXSeq methods

  time[['dexseq-raw']] <-
    system.time({res[['dexseq-raw']] <- runDEXSeq(targets = targets, quantifier = quantifier,count.type = 'raw', tx.gene = tx.gene)})

  time[['dexseq-scaledTPM']] <-
    system.time({res[['dexseq-scaledTPM']] <- runDEXSeq(targets = targets, quantifier = quantifier,count.type = 'scaledTPM', tx.gene = tx.gene)})

  time[['dexseq-dtuScaledTPM']] <-
    system.time({res[['dexseq-dtuScaledTPM']] <- runDEXSeq(targets = targets, quantifier = quantifier,count.type = 'dtuScaledTPM', tx.gene = tx.gene)})

  # Computing time

  time <- as.data.frame(do.call(rbind,time))
  time <- cbind('method' = rownames(time),time)

  return(list('res' = res, 'time' = time))
}

runMethods <- function(meta.path,quant.path,dest,quantifier){

  # Getting tx2gene info
  meta.path <- normalizePath(meta.path)
  tx.gene <- fread(file.path(meta.path,'counts.tsv.gz'),select = c('TranscriptID','GeneID')) |> as.data.frame()

  quant.path <- normalizePath(quant.path)
  sample.names <- basename(list.dirs(quant.path,full.names = TRUE,recursive = FALSE))
  sample.names.split <- strsplit2(sample.names,'_')

  dir.create(path = dest,recursive = TRUE,showWarnings = FALSE)
  dest <- normalizePath(dest)
  if (file.exists(file.path(dest, 'targets.tsv'))) {
    message('Files are already processed. Skipping them and moving on...')
    return(invisible(NULL))
  }

  targets <- data.frame(group = gsub('group','',sample.names.split[,1]),
                        replicate = gsub('rep','',sample.names.split[,2]))

  targets$group <- relevel(as.factor(targets$group),ref = 'A')
  targets$names <- sample.names
  targets$sample <- targets$names
  targets$path <- file.path(quant.path,targets$names)
  targets$sample_id <- targets$sample
  targets$condition <- as.numeric(targets$group)

  out <- callMethods(targets,quantifier,tx.gene)

  for (meth.name in names(out$res)) {
    saveRDS(object = out$res[[meth.name]],
            file = file.path(dest,paste0(meth.name,'.rds')),
            compress = 'xz')
  }
  write_tsv(x = out$time,file = file.path(dest,'time.tsv'),
            col_names = TRUE,quote = 'none')
  write_tsv(x = targets,file = file.path(dest,'targets.tsv'),
            col_names = TRUE,quote = 'none')

  return(invisible(NULL))
}

getGeneStatus <- function(GeneStatus){
  gs <- unique(GeneStatus)
  if (all(is.na(GeneStatus))) {
    out <- 'Not Expressed'
  } else{
    out <- gs[!is.na(gs)]
  }
  return(out)
}

#' @importFrom Biobase listLen
perGeneQValueExact <- function(pGene, theta, geneSplit) {
  # From DEXSeq package. This function is not exported.
  # This function is used in the satuRn pipeline.
  stopifnot(length(pGene) == length(geneSplit))
  numExons = listLen(geneSplit)
  tab = tabulate(numExons)
  notZero = (tab > 0)
  numerator = mapply(function(m, n) m * (1 - (1 - theta)^n),
                     m = tab[notZero], n = which(notZero))
  numerator = rowSums(numerator)
  bins = cut(pGene, breaks = c(-Inf, as.vector(theta)), right = TRUE,
             include.lowest = TRUE)
  counts = tabulate(bins, nbins = nlevels(bins))
  denom = cumsum(counts)
  stopifnot(denom[length(denom)] == length(pGene))
  return(numerator/denom)
}

no.na <- function(x){
  x[is.na(x)] <- 1
  x
}

#' @importFrom fishpond scaleInfReps labelKeep swish isoformProportions
#' @importFrom S4Vectors mcols
runSwish <- function(targets,quantifier,count.type,tx.gene){
  # Not run
  se <- getCounts(targets,tx.gene,quantifier,count.type)
  se <- scaleInfReps(se)
  se <- labelKeep(se)
  se <- se[mcols(se)$keep,]
  se <- isoformProportions(se)
  se <- swish(y = se, x = "group")
  out <- as.data.frame(mcols(se))
  out <- data.frame(TranscriptID = rownames(out),GeneID = unlist(out$gene_id),logFC = out$log2FC,P.Value = out$pvalue,adj.P.Value = out$qvalue)
  return(list('transcript' = out,'gene' = NULL))
}
