gammaToLogNormal <- function(x,shape,scale){
  # This function transforms a vector of gamma random variables to an equivalent
  # vector of log-normal random variables using quantile transformation. We 
  # implemented this function aiming to compare downstream DTE results between
  # scenarios in which the expression levels follow either a gamma or log-normal
  # distribution
  
  sigma2 <- log((shape + 1) / shape)
  mu <- 3 / 2 * log(shape) + log(scale) - 0.5 * log(shape + 1)
  
  ylnorm <- x
  pupper <- pgamma(x,shape = shape,scale = scale,lower.tail = FALSE,log.p = TRUE)
  plower <- pgamma(x,shape = shape,scale = scale,lower.tail = TRUE,log.p = TRUE)
  up <- pupper < plower
  if(any(up)) ylnorm[up] <- qlnorm(pupper[up],meanlog=mu,sdlog=sqrt(sigma2),lower.tail=FALSE,log.p=TRUE)
  if(any(!up)) ylnorm[!up] <- qlnorm(plower[!up],meanlog=mu,sdlog=sqrt(sigma2),lower.tail=TRUE,log.p=TRUE)
  
  return(ylnorm)
}

createDE <- function(Baseline,GeneStatus,fc){
  n <- length(Baseline)
  u <- runif(1)
  NewFC <- (fc * (u <= 0.5) + (u > 0.5) / fc)
  if(all(GeneStatus == 'Null')){
    NewBaseline <-  Baseline
  }
  if(all(GeneStatus == 'DGE')){
    NewBaseline <-  NewFC * Baseline
  }
  if(all(GeneStatus == 'DGE/DTU')){
    w <- sample(n,1,replace = FALSE)
    NewBaseline <- Baseline
    NewBaseline[w] <- NewFC * Baseline[w]
  }
  if(all(GeneStatus == 'DTU')){
    z <- sample(n,2,replace = FALSE)
    s <- sum(Baseline[z])
    NewBaseline <- Baseline
    NewBaseline[z] <- Baseline[z] <- c(NewFC*s/(1+NewFC),s/(1+NewFC))
    NewBaseline[z[2:1]] <- NewBaseline[z]
  }
  
  v <- runif(1)
  b1 <- (v < 0.5)*Baseline + (v >= 0.5)*NewBaseline
  b2 <- (v >= 0.5)*Baseline + (v < 0.5)*NewBaseline
  
  TranscriptStatus <- vector('numeric',n)
  TranscriptStatus[b2/b1 > 1] <- 1L # Up in group 2
  TranscriptStatus[b2/b1 < 1] <- -1L # Down in group 2
  
  list(Baseline.G1 = b1, Baseline.G2 = b2,TranscriptStatus = TranscriptStatus)
}

simulateExpr <- function(x,n.feat,n.libs,lib.sizes,num.DE,fc,lognormal,df.bcv = 40, bcv.true = 0.2){
  # This function was written with the goal of mimicking the simulation setup
  # used in the voom paper. Specifically, we use goodTuringProportions to
  # estimate baseline abundances, different asymptotic BCV value as well as
  # Chi-square degrees of freedom, and a different dispersion trend.
  # See ?baselineAbundance_function for goodTuringProportions usage in this
  # simulation.
  
  y <- as.data.table(copy(x))
  n.groups <- length(n.libs)
  n.samples <- sum(n.libs)
  if (n.groups > 2) stop('This function does not support more than 2 groups')
  
  # Generating baseline proportions
  goodTuring <- get('baselineAbundance_function')
  baselineAbundance <- goodTuring(seq_len(n.feat) / (n.feat + 1))
  baselineAbundance <- baselineAbundance/sum(baselineAbundance)
  y$Baseline[sample(nrow(y),nrow(y),replace = F)] <- baselineAbundance/sum(baselineAbundance)
  
  # Random allocating gene status
  y.gene <- y[,.(NTranscripts = .N), by = 'GeneID']
  y.gene[,Status := "Null"]
  y.gene[sample(.N,num.DE['DGE'],replace = FALSE), Status := "DGE"]
  y.gene[sample(which(Status == "Null" & NTranscripts > 1),num.DE['DTU'],replace = FALSE),Status := "DTU"]
  y.gene[sample(which(Status == "Null" & NTranscripts > 1),num.DE['DGE/DTU'],replace = FALSE),Status := "DGE/DTU"]
  y$GeneStatus <- y.gene$Status[match(y$GeneID,y.gene$GeneID)]
  
  y[,c('Baseline.G1','Baseline.G2','TranscriptStatus') := createDE(Baseline,GeneStatus,fc = fc),by = 'GeneID']
  baselineAbundancePerGroup <- as.matrix(y[,c('Baseline.G1','Baseline.G2')])
  de <- y$TranscriptStatus
  de.gene <- y$GeneStatus
  
  # Generating expected counts
  mu0 <- lapply(seq_len(n.groups),function(x){
    size <- lib.sizes[seq(1 + n.libs[x] * (x - 1), n.libs[x] * x)]
    matrix(baselineAbundancePerGroup[,x],n.feat,1) %*% matrix(size,1,length(size))
  })
  mu0 <- do.call(cbind,mu0)
  
  # Generating random noise around dispersion trend. I am generating 1 RV per
  # transcript per group. In the voom simulation, we had 1 RV per gene per sample.
  # Here, I am arguing that the expression level and dispersion should be exactly
  # the same among libraries of the same group. Since we already generated DE
  # status in the steps above, it makes sense to have 1 random shift around the
  # dispersion trend per group and a fixed resulting dispersion per group.
  chisq <- lapply(seq_len(n.groups),function(x){
    rv <- df.bcv / rchisq(n.feat, df = df.bcv)
    matrix(rv,ncol = n.libs[x],nrow = n.feat)
  })
  chisq <- do.call(cbind,chisq)
  
  # Biological variation and Dispersion trend
  bcv0 <- bcv.true + 1/sqrt(mu0)
  disp <- bcv0 ^ 2 * chisq
  
  # Biological variation
  shape <- 1/disp
  scale <- mu0/shape
  
  expr <- rgamma(n.feat * n.samples, shape = shape, scale = scale)
  
  if (isTRUE(lognormal)) {
    gammaToLogNormalVectorized <- Vectorize(gammaToLogNormal)
    expr <- gammaToLogNormalVectorized(expr,shape = c(shape),scale = c(scale))
  }
  
  expr <- matrix(expr,nrow =  n.feat, ncol =  n.samples)
  
  # Keying back
  key <- match(x$TranscriptID,y$TranscriptID)
  expr <- expr[key,]
  mu0 <- mu0[key,]
  disp <- disp[key,]
  de <- de[key]
  de.gene <- de.gene[key]
  
  return(list('expr' = expr,'mu' = mu0,'disp' = disp,'de' = de,'de.gene' = de.gene))
}

#' @importFrom data.table setkey copy setnames
simulateTPM <- function(contigs,contigs.subset,
                        n.libs,lib.sizes,
                        num.DE,fc,lognormal){
  
  # Generating sample labels
  group <- rep(LETTERS[seq_len(length(n.libs))],times = n.libs)
  rep <- unlist(lapply(n.libs,seq_len))
  group.name <- paste(paste0('group', group), paste0('rep',rep), sep = '_')
  
  # Simulating transcript-wise expression
  trExpr <- simulateExpr(x = contigs.subset,n.feat = nrow(contigs.subset),
                         n.libs = n.libs,lib.sizes = lib.sizes,
                         num.DE = num.DE,fc = fc,lognormal = lognormal)
  
  # Generating TPM values
  tpm <- trExpr$expr / contigs.subset$Length
  tpm <- 1e6 * t(t(tpm) / colSums(tpm))
  
  # Organizing contigs.subset
  dimnames(tpm) <- dimnames(trExpr$expr) <- dimnames(trExpr$mu) <- dimnames(trExpr$disp) <- list(contigs.subset$TranscriptID,group.name)
  names(trExpr$de) <-  names(trExpr$de.gene) <- contigs.subset$TranscriptID
  
  # Merging contigs.subset values to contigs
  tpm.contigs <- expr.contigs <- mu.contigs <- disp.contigs <- matrix(NA, nrow(contigs), sum(n.libs))
  de.contigs <- de.gene.contigs <- rep(NA,nrow(contigs))
  
  dimnames(tpm.contigs) <- dimnames(expr.contigs) <- dimnames(mu.contigs) <- dimnames(disp.contigs) <- list(contigs$TranscriptID,group.name)
  names(de.contigs) <- names(de.gene.contigs) <- contigs$TranscriptID
  
  key <- match(contigs.subset$TranscriptID, contigs$TranscriptID)
  tpm.contigs[key,] <- tpm
  expr.contigs[key,] <- trExpr$expr
  mu.contigs[key,] <- trExpr$mu
  disp.contigs[key,] <- trExpr$disp
  de.contigs[key] <- trExpr$de
  de.gene.contigs[key] <- trExpr$de.gene
  
  return(list('tpm' = tpm.contigs,'expr' = expr.contigs,'mu' = mu.contigs,'disp' = disp.contigs,'de' = de.contigs, 'de.gene' = de.gene.contigs))
}

#' @importFrom data.table as.data.table
selectTx <- function(contigs,genome){
  DT <- as.data.table(contigs)
  
  # Selecting reference ranking
  if (is.character(genome)) {
    if (genome == 'mm39') DT.ref <- get('GSE227750')
    if (genome == 'hg38') stop('simulation for human-like data not yet done')
  } else {
    DT.ref <- genome
  }
  
  # Subsetting transcripts
  DT.sub <- DT[TranscriptID %in% DT.ref$TranscriptID,]
  as.data.frame(DT.sub)
}

#' @importFrom limma strsplit2
readFasta <- function(fasta){
  # Reading input and getting tx information from fasta
  contigs <- scanFasta(fasta,quiet = TRUE)
  tx.info <- strsplit2(contigs$TranscriptID,"\\|")
  contigs[,c('TranscriptID','GeneID','TranscriptType')] <- tx.info[,c(1,2,8)]
  contigs
}

#' @importFrom Rsubread scanFasta simReads
#' @importFrom BiocParallel bplapply
#' @importFrom readr write_tsv
simulateFASTQ <- function(fasta,n.libs,lib.sizes,dest,tmpdir,paired.end,
                          fc,num.DE,genome,BPPARAM,read.length,
                          fragment.length.min,lognormal){
  
  # Checking if simulation has already been run
  dir.create(dest,showWarnings = FALSE,recursive = TRUE)
  out.files <- file.path(dest,c('counts.tsv.gz','targets.tsv.gz','tpm.tsv.gz','transcript.status.tsv.gz','gene.status.tsv.gz'))
  names(out.files) <- c('counts','targets','tpm','transcript','gene')
  if (all(file.exists(out.files))) {
    message('Reads already simulated!')
    return(invisible(NULL))
  }
  
  # Scanning fasta
  message('Preparing fasta file for simulation...')
  contigs <- readFasta(fasta)
  
  # Subsetting transcripts according to curated data and bringing metadata
  message('Subsetting transcripts...')
  contigs.subset <- selectTx(contigs,genome)
  
  # Simulating tx-wise TPMs from genes
  message('Simulating transcript-wise TPM...')
  txTPM <- simulateTPM(contigs = contigs, contigs.subset = contigs.subset,
                       lib.sizes = lib.sizes,n.libs = n.libs, 
                       num.DE = num.DE,fc = fc,lognormal = lognormal)
  
  # Getting quality reference
  quality.source <- ifelse(read.length %in% c(75, 100),'Rsubread','rfun')
  quality.reference <- list.files(system.file(package = quality.source,'qualf'),
                                  paste0('-',read.length,'bp'),full.names = TRUE)
  
  # Running simReads
  message('Simulating reads...')
  curwd <- getwd()
  setwd(tmpdir)
  out <- bplapply(seq_len(sum(n.libs)), FUN = function(i){
    simReads(transcript.file = fasta, expression.levels = txTPM$tpm[, i], 
             output.prefix = colnames(txTPM$tpm)[i],library.size = lib.sizes[i],
             paired.end = paired.end, simplify.transcript.names = TRUE,
             fragment.length.min = fragment.length.min,read.length = read.length,
             quality.reference = quality.reference)
  },BPPARAM = BPPARAM)
  out.fastq <- normalizePath(list.files('.',"group*.*.fastq.gz",full.names = TRUE))
  setwd(curwd)
  
  # Saving metadata
  message('Saving metadata...')
  mat <- do.call(cbind,lapply(out,function(x){x$NReads}))
  colnames(mat) <- colnames(txTPM$tpm)
  
  out <- cbind(contigs[,c("TranscriptID","Length",'GeneID')],mat)
  rownames(out) <- NULL
  write_tsv(x = out,file = out.files['counts'],col_names = TRUE,quote = 'none')
  
  out.tpm <- txTPM$tpm
  out.tpm[is.na(out.tpm)] <- 0L
  out.tpm <- cbind(contigs[,c("TranscriptID","Length",'GeneID')],out.tpm)
  rownames(out.tpm) <- NULL
  write_tsv(x = out.tpm,file = out.files['tpm'],col_names = TRUE,quote = 'none')
  
  out.transcript <- data.table(TranscriptID = out$TranscriptID,TranscriptStatus = txTPM$de)
  rownames(out.transcript) <- NULL
  write_tsv(x = out.transcript,file = out.files['transcript'],col_names = TRUE,quote = 'none')
  
  out.gene <- data.table(GeneID = out$GeneID,GeneStatus = txTPM$de.gene)
  out.gene <- out.gene[,.(GeneStatus = getGeneStatus(GeneStatus)),by = 'GeneID']
  out.gene[GeneStatus == 'Not Expressed', GeneStatus := NA]
  rownames(out.gene) <- NULL
  write_tsv(x = out.gene,file = out.files['gene'],col_names = TRUE,quote = 'none')
  
  targets <- data.frame('R1' = out.fastq[grepl("R1.fastq.gz",out.fastq)])
  if (isTRUE(paired.end)) {
    targets$R2 <- out.fastq[grepl("R2.fastq.gz",out.fastq)]
  }
  write_tsv(x = targets,file = out.files['targets'],col_names = TRUE,quote = 'none')
  
  message('Simulation of sequencing reads completed!')
}

#' @importFrom BiocParallel MulticoreParam register
simulateExperiment <- function(dest,
                               fasta,
                               genome,
                               bin.salmon,
                               index.salmon,
                               bin.kallisto,
                               index.kallisto,
                               tmpdir = tempdir(),
                               workers = 1,
                               num.DE = c('DGE' = 1500,'DTU' = 1500, 'DGE/DTU' = 1500),
                               fc = 2,
                               n.libs = c(3,3),
                               lib.sizes = rep(50e6,sum(n.libs)),
                               paired.end = FALSE,
                               keep.fastq = FALSE,
                               read.length = 75,
                               fragment.length.min = 150,
                               lognormal = FALSE,
                               opts.salmon =  paste('-p',workers,'-l A --numGibbsSamples 100 --validateMappings'),
                               opts.kallisto = paste0('--bootstrap-samples=100 --threads=',workers),
                               run.salmon = TRUE,
                               run.kallisto = FALSE,
                               run.dtu = FALSE,
                               seed = NULL){
  
  # Setting up parallel computing
  if(is.null(seed)){
    BPPARAM <- MulticoreParam(workers = workers, progressbar = TRUE)
  } else{
    set.seed(seed)
    BPPARAM <- MulticoreParam(workers = workers, progressbar = TRUE,RNGseed = seed)
  }
  register(BPPARAM = BPPARAM)
  
  # Preparing directories
  dir.create(tmpdir,recursive = TRUE,showWarnings = FALSE)
  dir.create(dest,recursive = TRUE,showWarnings = FALSE)
  
  tmpdir <- normalizePath(tmpdir)
  dest <- normalizePath(dest)
  
  # Simulating FASTQs
  simulateFASTQ(fasta = fasta,n.libs = n.libs,tmpdir = tmpdir,
                lib.sizes = lib.sizes,dest = file.path(dest,'meta'), fc = fc,
                paired.end = paired.end,num.DE = num.DE,
                genome = genome,BPPARAM = BPPARAM,
                read.length = read.length,fragment.length.min = fragment.length.min,
                lognormal = lognormal)
  
  # Quantifying FASTQs
  path.targets <- file.path(dest,'meta/targets.tsv.gz')
  quantifyReads(targets = path.targets,dest = dest,
                genome = genome,workers = workers, keep.fastq = keep.fastq,
                bin.salmon = bin.salmon,index.salmon = index.salmon,opts.salmon = opts.salmon,
                bin.kallisto = bin.kallisto,index.kallisto = index.kallisto,opts.kallisto = opts.kallisto,
                run.salmon = run.salmon, run.kallisto = run.kallisto)
  
  # Running methods
  runDTUMethods(dest = file.path(dest),run.salmon = run.salmon, run.kallisto = run.kallisto,run.dtu = run.dtu)
  
  # Organizing FASTQ files
  path.fastq <- read.delim(path.targets,header = TRUE)
  if (isTRUE(keep.fastq)) {
    message('Copying FASTQ files...')
    dir.create(file.path(dest,'fastq'))
    file.copy(from = unlist(path.fastq),to = file.path(dest,'fastq'))
  }
  unlink(tmpdir,recursive = TRUE)
}

runDTUMethods <- function(dest,run.salmon,run.kallisto,run.dtu){
  if (run.dtu) {
    if (run.salmon) {
      message('Running DTU methods with Salmon quantification...')
      dir.create(file.path(dest,'dtu-salmon'),recursive = TRUE,showWarnings = FALSE)
      runMethods(meta.path = file.path(dest,'meta'),quant.path = file.path(dest,'quant-salmon'),dest = file.path(dest,'dtu-salmon'),quantifier = 'salmon') 
      if (file.exists(file.path(dest, 'dtu-salmon', 'time.tsv'))) {
        message('DTU analysis w/ Salmon completed!')
      } else{
        stop('DTU analysis w/ Salmon failed!')
      }
    }
    if (run.kallisto) {
      message('Running DTU methods with kallisto quantification...')
      dir.create(file.path(dest,'dtu-kallisto'),recursive = TRUE,showWarnings = FALSE)
      runMethods(meta.path = file.path(dest,'meta'),quant.path = file.path(dest,'quant-kallisto'),dest = file.path(dest,'dtu-kallisto'),quantifier = 'kallisto') 
      if (file.exists(file.path(dest, 'dtu-kallisto', 'time.tsv'))) {
        message('DTU analysis w/ kallisto completed!')
      } else{
        stop('DTU analysis w/ kallisto failed!')
      }
    } 
  }
  return(invisible())
}

quantifyReads <- function(targets,
                          dest,
                          genome,
                          workers,
                          keep.fastq,
                          bin.salmon,
                          index.salmon,
                          opts.salmon,
                          bin.kallisto,
                          index.kallisto,
                          opts.kallisto,
                          run.salmon,
                          run.kallisto){
  
  dir.salmon <- file.path(dest,'quant-salmon')
  dir.kallisto <- file.path(dest,'quant-kallisto')
  
  df.targets <- read.delim(targets,header = TRUE)
  
  if(run.kallisto){
    message('Quantifying reads with kallisto...')
    dir.create(dir.kallisto,recursive = TRUE,showWarnings = FALSE)
    runKallisto(bin = bin.kallisto,
                index = index.kallisto,
                options = opts.kallisto,
                targets = df.targets, 
                dest = dir.kallisto)
    message('Quantification with kallisto completed!')
  }
  
  if(run.salmon){
    message('Quantifying reads with Salmon...')
    dir.create(dir.salmon,recursive = TRUE,showWarnings = FALSE)
    runSalmon(bin = bin.salmon,
              index = index.salmon,
              options = opts.salmon,
              targets = df.targets,
              dest = dir.salmon)
    message('Quantification with Salmon completed!')
  }
}