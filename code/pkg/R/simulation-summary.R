#' @importFrom grDevices rgb colorRamp
colramp <- function(n,x){
  y <- colorRamp(x)
  y <- y(n)
  apply(y,1,function(z){rgb(z[1],z[2],z[3],maxColorValue = 255)})
}

#' @importFrom RColorBrewer brewer.pal
methodsNames <- function(){
  method <- c('edger.v3-scaled-simes',
              'edger.v3-raw-ftest',
              'edger.v3-raw-simes',
              'edger.v3-scaled-ftest',
              'edger.v4-scaled-simes',
              'edger.v4-raw-ftest',
              'edger.v4-raw-simes',
              'edger.v4-scaled-ftest',
              'limma-scaled-simes',
              'limma-raw-ftest',
              'limma-raw-simes',
              'limma-scaled-ftest',
              'drimseq-raw',
              'drimseq-scaledTPM',
              'drimseq-dtuScaledTPM',
              'saturn-raw',
              'saturn-scaledTPM',
              'saturn-dtuScaledTPM',
              'dexseq-raw',
              'dexseq-scaledTPM',
              'dexseq-dtuScaledTPM')

  labels <- c('edgeR.v3-scaled-Simes',
              'edgeR.v3-raw-F',
              'edgeR.v3-raw-Simes',
              'edgeR.v3-scaled-F',
              'edgeR-scaled-Simes',
              'edgeR-raw-F',
              'edgeR-raw-Simes',
              'edgeR-scaled-F',
              'limma-scaled-Simes',
              'limma-raw-F',
              'limma-raw-Simes',
              'limma-scaled-F',
              'DRIMSeq-raw',
              'DRIMSeq-scaledTPM',
              'DRIMSeq-dtuScaledTPM',
              'satuRn-raw',
              'satuRn-scaledTPM',
              'satuRn-dtuScaledTPM',
              'DEXSeq-raw',
              'DEXSeq-scaledTPM',
              'DEXSeq-dtuScaledTPM')

  c4 <- c(0,0.25,0.75,1)
  c3 <- c(0,0.5,1)

  color <- c(colramp(c4,c("orange", "lightgoldenrod")),
             colramp(c4,c("red", "lightsalmon")),
             colramp(c4,c("blue", "lightblue")),
             colramp(c3,c("green3", "lightgreen")),
             colramp(c3,c("purple", "lightpink")),
             colramp(c3,c("black", "lightgray")))

  names(method) <- names(color) <- labels
  return(list(labels = labels,method = method,color = color))
}

roundPretty <- function(x,digits = 1){
  formatC(round(x,digits),digits = digits,format = 'f')
}

loadRDS <- function(name,type,path){
  rds <- readRDS(path[name])
  dt <- data.table('Method' = names(path[name]), as.data.table(rds[[type]]))
  setnames(x = dt,
           old = c('pvalue','pval','P.Value','qvalue','qval'),
           new = c('PValue','PValue','PValue','FDR','FDR'),skip_absent = TRUE)
  feature.name <- ifelse(type == 'gene','GeneID','TranscriptID')
  dt <- dt[,c('Method',feature.name,'PValue','FDR'),with = FALSE]

  return(dt)
}

#' @importFrom data.table data.table fread
loadResults <- function(path,genome,len,fc,read,scenario,libs.per.group,simulation,quantifier){

  meth <- methodsNames()
  path.time <- file.path(path,'time.tsv')
  path.method <- file.path(path,paste0(meth$method,'.rds'))
  names(path.method) <- meth$labels

  dt.scenario <- data.table('Genome' = genome,
                            'Length' = len,
                            'FC' = fc,
                            'Reads' = read,
                            'Scenario' = scenario,
                            'LibsPerGroup' = libs.per.group,
                            'Simulation' = simulation,
                            'Quantifier' = quantifier)

  dt.results.gene <- lapply(names(path.method),loadRDS,type = 'gene',path = path.method)
  dt.results.transcript <- lapply(names(path.method),loadRDS,type = 'transcript',path = path.method)

  dt.results.gene <- cbind(dt.scenario,do.call(rbind,dt.results.gene))
  dt.results.transcript <- cbind(dt.scenario,do.call(rbind,dt.results.transcript))

  # For transcript-level tests, drop one of F or Simes rows, because they both give
  # the same result as there is only one transcript-level test (t-test in limma
  # and 'exon' test in edgeR)
  dt.results.transcript <- dt.results.transcript[!grepl("Simes",Method,ignore.case = TRUE),]

  dt.time <- fread(file = path.time,header = TRUE)
  dt.time <- cbind(dt.scenario,dt.time[,c('method','elapsed')])
  setnames(x = dt.time,old = c('method','elapsed'),new = c('Method','Time'))
  dt.time$Method <- names(meth$method)[match(dt.time$Method,meth$method)]

  out <- list('results.gene' = dt.results.gene,'results.transcript' = dt.results.transcript,'time' = dt.time)

  return(out)
}

loadMetadata <- function(path,genome,len,fc,read,scenario,libs.per.group,simulation){
  path.gene.status <- file.path(path,'gene.status.tsv.gz')
  path.transcript.status <- file.path(path,'transcript.status.tsv.gz')
  path.counts <- file.path(path,'counts.tsv.gz')

  dt.scenario <- data.table('Genome' = genome,
                            'Length' = len,
                            'FC' = fc,
                            'Reads' = read,
                            'Scenario' = scenario,
                            'LibsPerGroup' = libs.per.group,
                            'Simulation' = simulation)

  dt.gene.metadata <- fread(path.gene.status,header = TRUE)
  dt.transcript.metadata <- fread(path.transcript.status,header = TRUE)
  dt.counts <- fread(path.counts,select = c('TranscriptID','GeneID'),header = TRUE)

  dt.gene <- dt.gene.metadata[,c('GeneID')]
  dt.transcript <- dt.transcript.metadata[,c('TranscriptID')]

  dt.gene.sim <- cbind(dt.scenario,dt.gene.metadata[dt.gene.metadata$GeneStatus %in% c("DGE/DTU","DTU"), 'GeneID'])

  is.dtu.transcript <- dt.transcript.metadata$TranscriptStatus %in% c(-1,1) &
    dt.transcript.metadata$TranscriptID %in% dt.counts[GeneID %in% dt.gene.sim$GeneID,TranscriptID]
  dt.transcript.sim <- cbind(dt.scenario,dt.transcript.metadata[is.dtu.transcript,])

  # If fold-change = 1, status should be 0
  if (fc == 'fc1') {
    dt.gene.sim$GeneStatus <- 0L
    dt.transcript.sim$TranscriptStatus <- 0L
  } else{
    dt.gene.sim$GeneStatus <- 1L
    dt.transcript.sim$TranscriptStatus <- 1L
  }

  out <- list('simulation.gene' = dt.gene.sim,
              'simulation.transcript' = dt.transcript.sim,
              'gene' = dt.gene,
              'transcript' = dt.transcript)

  return(out)
}

aggregateScenario <- function(path,genome,len,fc,read,scenario,libs.per.group,quantifier,nsim){

  subpath <- paste0('simulation-',seq_len(nsim))

  ls.results <- lapply(seq_len(nsim),function(x){
    res.path <- file.path(path,subpath[x],paste0('dtu-',quantifier))
    loadResults(res.path,genome,len,fc,read,scenario,libs.per.group,x,quantifier)
  })

  ls.metadata <- lapply(seq_len(nsim),function(x){
    meta.path <- file.path(path,subpath[x],'meta')
    loadMetadata(meta.path,genome,len,fc,read,scenario,libs.per.group,x)
  })

  dt.results.gene <- do.call(rbind,lapply(ls.results,function(x){x[['results.gene']]}))
  dt.results.transcript <- do.call(rbind,lapply(ls.results,function(x){x[['results.transcript']]}))
  dt.time <- as.data.table(do.call(rbind,lapply(ls.results,function(x){x[['time']]})))
  dt.simulation.gene <- do.call(rbind,lapply(ls.metadata,function(x){x[['simulation.gene']]}))
  dt.simulation.transcript <- do.call(rbind,lapply(ls.metadata,function(x){x[['simulation.transcript']]}))
  dt.gene <- ls.metadata[[1]]$gene
  dt.transcript <- ls.metadata[[1]]$transcript

  out <- list('results.gene' = dt.results.gene,'results.transcript' = dt.results.transcript,
              'simulation.gene' = dt.simulation.gene,'simulation.transcript' = dt.simulation.transcript,
              'features.gene' = dt.gene,'features.transcript' = dt.transcript,
              'time' = dt.time)
  return(out)
}

computeGeneMetrics <- function(x,simulation,features,fdr,alpha){

  GeneID.DE <- x$GeneID[x$FDR < fdr]
  n <- length(x$GeneID)
  n.lt.alpha <- sum(x$PValue < alpha)

  call.DE <- data.table(GeneID = GeneID.DE,call = 1)

  truth.DE <- simulation[Genome == x$Genome & Length == x$Length & FC == x$FC &
                           Reads == x$Reads &
                           Scenario == x$Scenario & LibsPerGroup == x$LibsPerGroup &
                           Simulation == x$Simulation, c('GeneID','GeneStatus')]

  tb.DE <- merge(features,truth.DE,by = 'GeneID',all.x = TRUE)
  tb.DE <- merge(tb.DE,call.DE,by = 'GeneID',all.x = TRUE)

  tb.DE[is.na(GeneStatus),GeneStatus := 0]
  tb.DE[, GeneStatus := abs(GeneStatus)]
  tb.DE[is.na(call),call := 0]

  tb.DE$GeneStatus <- factor(tb.DE$GeneStatus,levels = c(0,1))
  tb.DE$call <- factor(tb.DE$call,levels = c(0,1))

  tb.results <- tb.DE[,table(GeneStatus,call)]

  total.call <- sum(tb.results[,"1"])
  total.de <- sum(tb.results["1",])

  out <- list('N' = n,
              'N.ALPHA' = n.lt.alpha,
              'TP' = tb.results["1","1"],
              'FP' = tb.results["0","1"],
              'FPR' = tb.results["0","1"]/sum(tb.results["0",]),
              'FDR' = ifelse(total.call == 0,NA,tb.results["0","1"]/total.call),
              'TPR' = ifelse(total.de == 0,NA,tb.results["1","1"]/total.de))

  return(lapply(out,as.double))
}

computeGeneROCCurve <- function(x,simulation,features,fdr,seq.fdr){
  truth.DE <- simulation[Genome == x$Genome & Length == x$Length & FC == x$FC &
                           Reads == x$Reads &
                           Scenario == x$Scenario & LibsPerGroup == x$LibsPerGroup &
                           Simulation == x$Simulation, c('GeneID','GeneStatus')]
  tb.DE <- merge(features,truth.DE,by = 'GeneID',all.x = TRUE)
  tb.DE[is.na(GeneStatus),GeneStatus := 0]
  tb.DE[, GeneStatus := abs(GeneStatus)]

  feature.DE <- data.table(GeneID = x$GeneID,FDR = x$FDR)

  tb.DE <- merge(tb.DE,feature.DE,by = 'GeneID',all.x = TRUE)
  tb.DE[,call := 0]
  tb.DE$GeneStatus <- factor(tb.DE$GeneStatus,levels = c(0,1))
  tb.DE$call <- factor(tb.DE$call,levels = c(0,1))
  tb.DE[is.na(FDR),FDR := 1]

  out <- lapply(seq.fdr,function(w){
    tb.results <- copy(tb.DE)
    tb.results[FDR < w/100, call := "1"]
    tb.results <- tb.results[,table(GeneStatus,call)]

    num.de <- sum(tb.results[2,])
    num.call <- sum(tb.results[,2])

    return(c('TPR' = ifelse(num.de == 0,NA,tb.results[2,2]/sum(tb.results[2,])),
             'FDR' = ifelse(num.call == 0,NA,tb.results[1,2]/sum(tb.results[,2]))))
    # o <- c('TPR' = tb.results[2,2]/num.de,'FDR' = tb.results[1,2]/num.call)
    #
    # if(num.de > 0 & num.call == 0){
    #   o['TPR'] <- o['FDR'] <- 0L
    # }
    # if(num.de == 0){
    #   if(num.call > 0){
    #     o['TPR'] <- 0
    #     o['FDR'] <- 1L
    #   } else{
    #     o['TPR'] <- o['FDR'] <- 0L
    #   }
    # }
    # return(o)
  })

  out.tpr <- lapply(out,function(x){as.numeric(x['TPR'])})
  out.fdr <- lapply(out,function(x){as.numeric(x['FDR'])})

  names(out.tpr) <- paste0('tpr.',seq.fdr)
  names(out.fdr) <- paste0('fdr.',seq.fdr)
  return(c(out.tpr,out.fdr))
}

computeTranscriptROCCurve <- function(x,simulation,features,fdr,seq.fdr){
  truth.DE <- simulation[Genome == x$Genome & Length == x$Length & FC == x$FC &
                           Reads == x$Reads &
                           Scenario == x$Scenario & LibsPerGroup == x$LibsPerGroup &
                           Simulation == x$Simulation, c('TranscriptID','TranscriptStatus')]

  tb.DE <- merge(features,truth.DE,by = 'TranscriptID',all.x = TRUE)
  tb.DE[is.na(TranscriptStatus),TranscriptStatus := 0]
  tb.DE[, TranscriptStatus := abs(TranscriptStatus)]

  feature.DE <- data.table(TranscriptID = x$TranscriptID,FDR = x$FDR)

  tb.DE <- merge(tb.DE,feature.DE,by = 'TranscriptID',all.x = TRUE)
  tb.DE[,call := 0]
  tb.DE$TranscriptStatus <- factor(tb.DE$TranscriptStatus,levels = c(0,1))
  tb.DE$call <- factor(tb.DE$call,levels = c(0,1))
  tb.DE[is.na(FDR),FDR := 1]

  out <- lapply(seq.fdr,function(w){
    tb.results <- copy(tb.DE)
    tb.results[FDR < w/100, call := "1"]
    tb.results <- tb.results[,table(TranscriptStatus,call)]

    num.call <- sum(tb.results[,2])
    num.de <- sum(tb.results[2,])

    return(c('TPR' = ifelse(num.de == 0,NA,tb.results[2,2]/sum(tb.results[2,])),
             'FDR' = ifelse(num.call == 0,NA,tb.results[1,2]/sum(tb.results[,2]))))
    # o <- c('TPR' = tb.results[2,2]/num.de,'FDR' = tb.results[1,2]/num.call)
    #
    # if(num.de > 0 & num.call == 0){
    #   o['TPR'] <- o['FDR'] <- 0L
    # }
    # if(num.de == 0){
    #   if(num.call > 0){
    #     o['TPR'] <- 0
    #     o['FDR'] <- 1L
    #   } else{
    #     o['TPR'] <- o['FDR'] <- 0L
    #   }
    # }
    # return(o)
  })

  out.tpr <- lapply(out,function(x){as.numeric(x['TPR'])})
  out.fdr <- lapply(out,function(x){as.numeric(x['FDR'])})

  names(out.tpr) <- paste0('tpr.',seq.fdr)
  names(out.fdr) <- paste0('fdr.',seq.fdr)

  return(c(out.tpr,out.fdr))
}

computeTranscriptMetrics <- function(x,simulation,features,fdr,alpha){

  TranscriptID.DE <- x$TranscriptID[x$FDR < fdr]
  n <- length(x$TranscriptID)
  n.lt.alpha <- sum(x$PValue < alpha)

  call.DE <- data.table(TranscriptID = TranscriptID.DE,call = 1)

  truth.DE <- simulation[Genome == x$Genome & Length == x$Length & FC == x$FC &
                           Reads == x$Reads &
                           Scenario == x$Scenario & LibsPerGroup == x$LibsPerGroup &
                           Simulation == x$Simulation, c('TranscriptID','TranscriptStatus')]

  tb.DE <- merge(features,truth.DE,by = 'TranscriptID',all.x = TRUE)
  tb.DE <- merge(tb.DE,call.DE,by = 'TranscriptID',all.x = TRUE)

  tb.DE[is.na(TranscriptStatus),TranscriptStatus := 0]
  tb.DE[, TranscriptStatus := abs(TranscriptStatus)]
  tb.DE[is.na(call),call := 0]

  tb.DE$TranscriptStatus <- factor(tb.DE$TranscriptStatus,levels = c(0,1))
  tb.DE$call <- factor(tb.DE$call,levels = c(0,1))

  tb.results <- tb.DE[,table(TranscriptStatus,call)]

  total.call <- sum(tb.results[,"1"])
  total.de <- sum(tb.results["1",])

  out <- list('N' = n,
              'N.ALPHA' = n.lt.alpha,
              'TP' = tb.results["1","1"],
              'FP' = tb.results["0","1"],
              'FPR' = tb.results["0","1"]/sum(tb.results["0",]),
              'FDR' = ifelse(total.call == 0,NA,tb.results["0","1"]/total.call),
              'TPR' = ifelse(total.de == 0,NA,tb.results["1","1"]/total.de))

  return(lapply(out,as.double))
}

computeGeneFDRCurve <- function(x,simulation,features,fdr,seq.n){

  truth.DE <- simulation[Genome == x$Genome & Length == x$Length & FC == x$FC &
                           Reads == x$Reads &
                           Scenario == x$Scenario & LibsPerGroup == x$LibsPerGroup &
                           Simulation == x$Simulation, c('GeneID','GeneStatus')]

  tb.DE <- merge(features,truth.DE,by = 'GeneID',all.x = TRUE)
  tb.DE[is.na(GeneStatus),GeneStatus := 0]
  tb.DE[, GeneStatus := abs(GeneStatus)]
  
  if(grepl("edgeR|limma|DRIMSeq",x$Method)){
    ranking.variable <- x$PValue
  } else{
    ranking.variable <- x$FDR
  }

  feature.DE <- data.table(GeneID = x$GeneID,RankingVar = ranking.variable,call = 1)

  tb.DE <- merge(tb.DE,feature.DE,by = 'GeneID',all.x = TRUE)
  tb.DE[is.na(call),call := 0]
  tb.DE$GeneStatus <- factor(tb.DE$GeneStatus,levels = c(0,1))
  tb.DE$call <- factor(tb.DE$call,levels = c(0,1))
  tb.DE <- tb.DE[order(RankingVar),]

  out <- lapply(seq.n,function(w){
    tb.results <- tb.DE[seq(1,w),][,table(GeneStatus,call)]
    return(tb.results["0","1"])
  })

  names(out) <- paste0('n.',seq.n)

  return(out)
}

computeTranscriptFDRCurve <- function(x,simulation,features,fdr,seq.n){

  truth.DE <- simulation[Genome == x$Genome & Length == x$Length & FC == x$FC &
                           Reads == x$Reads &
                           Scenario == x$Scenario & LibsPerGroup == x$LibsPerGroup &
                           Simulation == x$Simulation, c('TranscriptID','TranscriptStatus')]

  tb.DE <- merge(features,truth.DE,by = 'TranscriptID',all.x = TRUE)
  tb.DE[is.na(TranscriptStatus),TranscriptStatus := 0]
  tb.DE[, TranscriptStatus := abs(TranscriptStatus)]
  
  if(grepl("edgeR|limma|DRIMSeq",x$Method)){
    ranking.variable <- x$PValue
  } else{
    ranking.variable <- x$FDR
  }

  feature.DE <- data.table(TranscriptID = x$TranscriptID,RankingVar = ranking.variable,call = 1)

  tb.DE <- merge(tb.DE,feature.DE,by = 'TranscriptID',all.x = TRUE)
  tb.DE[is.na(call),call := 0]
  tb.DE$TranscriptStatus <- factor(tb.DE$TranscriptStatus,levels = c(0,1))
  tb.DE$call <- factor(tb.DE$call,levels = c(0,1))
  tb.DE <- tb.DE[order(RankingVar),]

  out <- lapply(seq.n,function(w){
    tb.results <- tb.DE[seq(1,w),][,table(TranscriptStatus,call)]
    return(tb.results["0","1"])
  })

  names(out) <- paste0('n.',seq.n)

  return(out)
}

#' @importFrom thematic okabe_ito
#' @importFrom ggplot2 ggplot geom_line theme_bw scale_color_manual element_rect
#' @importFrom ggplot2 scale_x_continuous theme element_blank labs aes alpha geom_point
#' @importFrom ggplot2 scale_y_continuous geom_abline facet_grid vars unit coord_cartesian
plotFDRCurve <- function(x,max.n,base_size = 8,xlab = 'Genes chosen'){

  meth <- methodsNames()

  plot <- ggplot(x,aes(x = N,y = FDR,color = Method,group = Method)) +
    facet_grid(rows = vars(LibsPerGroup),scales = 'free_y') +
    geom_line(size = 0.75) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    scale_color_manual(values = meth$color) +
    coord_cartesian(xlim = c(0,max.n)) +
    theme(strip.text = element_blank(),
          strip.background = element_blank(),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"line"),
          legend.position = 'top',
          panel.grid = element_blank(),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size)) +
    labs(y = 'False discoveries',x = xlab)

  return(plot)
}

#' @importFrom ggplot2 geom_vline scale_shape_manual
plotROCCurve <- function(x,seq.fdr,base_size = 8,max.x = 0.35){
  seq.fdr <- seq.fdr/100
  meth <- methodsNames()

  shape <- 21:25
  names(shape) <- as.character(seq.fdr)

  y <- copy(x)

  y$fill <- meth$color[match(y$Method,names(meth$color))]
  y[oFDR>nFDR,fill := "#FFFFFF"]

  ggplot(data = y,aes(x = oFDR,y = oTPR,group = Method,color = Method)) +
    geom_vline(xintercept = seq.fdr,linetype = 'dashed',colour = 'gray') +
    geom_line() +
    geom_point(aes(shape = as.character(nFDR)),fill = y$fill,size = 3) +
    scale_color_manual(values = meth$color) +
    scale_shape_manual(values = shape) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_continuous(limits = c(0,max.x),breaks = seq(0,max.x,0.05)) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          # legend.title = element_blank(),
          legend.key.size = unit(0.75,"line"),
          legend.position = 'top',
          panel.grid = element_blank(),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size)) +
    labs(y = 'Detection power',x = 'Observed FDR',shape = 'FDR')

}

#' @importFrom ggplot2 geom_col geom_text scale_fill_brewer .pt scale_fill_manual
plotPowerBars <- function(x,fdr,max.n,base_size = 8){

  sub.byvar <- colnames(x)[-which(colnames(x) %in% c('P.SIG','TP','FP'))]

  gap <- 0.05*max(x$TP + x$FP)

  x.melt <- melt(x,id.vars = sub.byvar,
                 measure.vars = c('TP','FP'),
                 variable.name = 'Type',
                 value.name = 'Value')
  x.melt$Type <- factor(x.melt$Type,levels = c('FP','TP'),labels = c('False Positive','True Positive'))

  plot <- ggplot(x.melt,aes(x = Method,y = Value,fill = Type)) +
    facet_grid(rows = vars(LibsPerGroup)) +
    geom_col() +
    geom_text(aes(x = Method,y = (TP + FP) + gap,
                  label = roundPretty(ifelse((FP+TP) == 0,NA,100*FP/(FP+TP)),1)),
              inherit.aes = FALSE,data = x,vjust = 0,size = base_size/.pt) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    scale_fill_manual(values = c('False Positive' = 'salmon','True Positive' = 'lightblue')) +
    labs(x = NULL,y = paste('# DE Transcripts at FDR <',roundPretty(fdr,2))) +
    scale_y_continuous(limits = c(0,max.n)) +
    theme(strip.text = element_blank(),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'top',
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 90,colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.key.size = unit(0.75,"line"))

  return(plot)
}

plotType1Error <- function(x,alpha,base_size = 8){

  sub.byvar <- colnames(x)[-which(colnames(x) %in% c('P.SIG','TP','FP'))]

  x.melt <- melt(x,id.vars = sub.byvar,
                 measure.vars = c('P.SIG'),variable.name = 'Type',value.name = 'Value')

  plot <- ggplot(x.melt,aes(x = Method,y = Value)) +
    facet_grid(rows = vars(LibsPerGroup)) +
    geom_col(fill = 'grey') +
    theme_bw(base_size = base_size,base_family = 'sans') +
    geom_hline(yintercept = alpha,color = 'red',linetype = 'dashed') +
    labs(x = NULL,y = paste('Proportion of genes with p-value <',roundPretty(alpha,2))) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),#element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"line"),
          legend.position = 'top',
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 90),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size))

  return(plot)
}

summarizePValue <- function(x,byvar,step = 0.05){
  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  cut.sq <- seq(step,1 - step,by = step)
  cut.match <- c(0,cut.sq) + step/2
  names(cut.match) <- paste0('(',roundPretty(c(0,cut.sq),2),'-',roundPretty(c(cut.sq,1),2),']')
  n.groups <- length(cut.match)

  x.sub <- copy(x)
  x.sub[,PValue := cut(x = PValue,breaks = c(-Inf,cut.sq,Inf),labels = names(cut.match))]

  x.sub.method <- x.sub[,list(N = .N),by = byvar]

  table <- x.sub[,list(N.cat = .N),by = c(byvar,'PValue')]
  table <- merge(table,x.sub.method,by = byvar,all.x = TRUE)
  table$PValue <- factor(table$PValue,levels = names(cut.match))

  table <- table[,list(Density.Avg = mean(n.groups*N.cat/N)),by = c(sub.byvar,'PValue')]

  table$PValue.Midpoint <- cut.match[match(table$PValue,names(cut.match))]

  return(table)
}

#' @importFrom ggplot2 facet_wrap geom_histogram geom_hline scale_x_discrete rel
plotPValues <- function(x,base_size = 8){
  plot <- ggplot(data = x,aes(x = PValue,y = Density.Avg)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(Method)) +
    geom_col(col = 'black',fill = 'grey',position = position_dodge(0.9),width = 0.8) +
    geom_hline(yintercept = 1,col = 'red',linetype = 'dashed') +
    theme_bw(base_size = base_size,base_family = 'sans') +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          strip.text.y = element_blank(),
          strip.background.y = element_blank(),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"line"),
          legend.position = 'top',
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size)) +
    labs(x = 'P-values',y = 'Density')
  return(plot)
}

#' @importFrom ggplot2 geom_bar position_dodge element_text
plotTime <- function(x,base_size = 8){
  plot <- ggplot(data = x,aes(x = Method,y = Time)) +
    facet_grid(rows = vars(LibsPerGroup)) +
    geom_bar(stat = 'identity',position = position_dodge()) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    scale_y_continuous(limits = c(0,max(5,max(x$Time)))) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = 'top',
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 90,colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size),
          strip.text = element_text(colour = 'black',size = base_size)) +
    labs(y = 'Time (min)',x = NULL)
  return(plot)
}

#' @importFrom data.table melt
summarizeFDRCurve <- function(x,byvar){

  cnames <- colnames(x)

  sub.cnames <- cnames[grepl('n\\.',cnames)]
  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  x.mean <- x[,lapply(.SD,mean),by = sub.byvar,.SDcols = sub.cnames]

  table <- melt(data = x.mean,id.vars = sub.byvar,variable.name = 'N',value.name = 'FDR')

  table[,N := as.numeric(gsub('n\\.','',N))]

  return(table)
}

summarizeMetrics <- function(x,byvar){

  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  table <- x[,.(P.SIG = mean(N.ALPHA/N),TP = mean(TP),FP = mean(FP)),sub.byvar]

  return(table)
}

summarizeTime <- function(x,byvar){

  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  table <- x[,.(Time = mean(Time/60)),sub.byvar]

  return(table)
}

summarizeQQ <- function(x,byvar,step = 0.001){

  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  # Only DEXSeq and satuRn have NA p-values: they do not output raw gene-level
  # p-values and NAs have been explicitly set for those
  if(x[!grepl("DEXSeq|satuRn",Method),any(is.na(PValue))]) stop('NA p-values')

  x.quant <- x[, list(q.sample = quantile(PValue,probs = seq(0,1,length.out = .N),na.rm = TRUE),
                      q.theory = seq(0,1,length.out = .N)),by = byvar]

  quant.sq <- seq(step,1 - step,by = step)
  quant.match <- c(0,quant.sq) + step/2
  names(quant.match) <- paste0('(',c(0,quant.sq),'-',c(quant.sq,1),']')

  x.quant[,Q.Theory.Cat := cut(x = q.theory,breaks = c(-Inf,quant.sq,Inf),labels = names(quant.match))]

  table <- x.quant[,.(Q.Sample.Avg = mean(q.sample,na.rm = TRUE)),by = c(sub.byvar,'Q.Theory.Cat')]

  table$Q.Theory.Midpoint <- quant.match[match(table$Q.Theory.Cat,names(quant.match))]

  return(table)
}

plotQQPlot <- function(x,base_size = 8){
  meth <- methodsNames()

  plot <- ggplot(x,
                 aes(x = Q.Theory.Midpoint,y = Q.Sample.Avg,color = Method,group = Method)) +
    facet_grid(rows = vars(LibsPerGroup)) +
    # geom_abline(intercept = 0,slope = 1,colour = 'black',linetype = 'dashed') +
    geom_line(inherit.aes = FALSE,aes(x=0,y = 0,color = Method,group = Method)) +
    geom_point(pch = '.',size = 2) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    scale_color_manual(values = meth$color) +
    scale_x_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1)) +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          # legend.key.size = unit(2,"line"),
          legend.position = 'top',
          panel.grid = element_blank(),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size)) +
    labs(y = 'Sample Quantiles',x = 'Theoretical Quantiles')

  return(plot)
}

summarizeOverdispersion <- function(path, genome, len, fc, read,
                                    scenario, libs.per.group, quantifier, nsim){

  subpath <- paste0('simulation-',seq_len(nsim))
  catchFunction <- get(ifelse(quantifier == 'salmon','catchSalmon','catchKallisto'))

  ls.results <- lapply(seq_len(nsim),function(x){
    res.path <- file.path(path,subpath[x],paste0('quant-',quantifier))
    meta.path <- file.path(path,subpath[x],'meta/transcript.status.tsv.gz')

    meta <- fread(meta.path)

    catch <- catchFunction(list.dirs(res.path,recursive = FALSE))
    rownames(catch$annotation) <- strsplit2(rownames(catch$annotation),"\\|")[,1]

    keep <- rownames(catch$annotation) %in% meta$TranscriptID[!is.na(meta$TranscriptStatus)]

    out <- catch$annotation$Overdispersion[keep]
    return(out)
  })

  res <- log10(unlist(ls.results))

  out <- data.table('Genome' = genome,
                    'Length' = len,
                    'FC' = fc,
                    'Reads' = read,
                    'Scenario' = scenario,
                    'LibsPerGroup' = libs.per.group,
                    'Quantifier' = quantifier,
                    'Mean' = mean(res),
                    'SD' = sd(res),
                    '2.5Pct' = quantile(res,0.025),
                    '25Pct' = quantile(res,0.25),
                    '50Pct' = quantile(res,0.5),
                    '75Pct' = quantile(res,0.75),
                    '97.5Pct' = quantile(res,0.975))

  return(out)
}


#' @importFrom kableExtra kbl add_header_above kable_styling pack_rows cell_spec
#' @importFrom kableExtra landscape collapse_rows
tabulateMetrics <- function(x,cap,
                            seq.len = seq(50,150,25),
                            lib.group = c(3,5,10),
                            lib.size = c('50M','25/100M'),
                            color = TRUE,
                            font_size = NULL,
                            color.fdr = 0.05,
                            format = 'latex',...){

  dt <- copy(x)

  methods <- methodsNames()$labels

  dt$Length %<>% factor(levels = paste0(seq.len,'bp'))
  dt$LibsPerGroup %<>% mapvalues(from = paste0("#Lib/Group = ",lib.group),to = lib.group)
  dt$Scenario %<>% mapvalues(from = c('balanced','unbalanced'),to = lib.size) %<>% factor(levels = lib.size)

  dt[,A.Power := TP/3000]
  dt[,B.FDR := ifelse((FP+TP) == 0,NA,FP/(FP+TP))]

  dt.dcast <- dcast(dt, Reads + LibsPerGroup + Scenario + Length ~ Method,value.var = c('A.Power','B.FDR'))

  dt.dcast <- dt.dcast[order(Reads,LibsPerGroup,Scenario,Length),]
  dt.dcast <- dt.dcast[,c('Reads','LibsPerGroup','Scenario','Length',paste0('A.Power_',methods),paste0('B.FDR_',methods)),with = FALSE]

  dt.dcast[,c(paste0('A.Power_',methods),paste0('B.FDR_',methods)) := lapply(.SD,roundPretty,digits = 3),.SDcols = c(paste0('A.Power_',methods),paste0('B.FDR_',methods))]

  dt.dcast.color <- copy(dt.dcast)

  if(color == TRUE){
    mat.power <- as.matrix(dt.dcast[,paste0('A.Power_',methods),with = FALSE])
    class(mat.power) <- 'numeric'
    mat.fdr <- as.matrix(dt.dcast[,paste0('B.FDR_',methods),with = FALSE])
    class(mat.fdr) <- 'numeric'

    col.power <- t(sapply(1:nrow(mat.power),FUN = function(i){
      ifelse(mat.power[i,] == max(mat.power[i,which(mat.fdr[i,] < color.fdr)]) &
               mat.fdr[i,] < color.fdr,'blue','black')
    }))
    col.power[is.na(col.power)] <- 'black'
    col.fdr <- t(apply(mat.fdr,1,function(x){ifelse(x > color.fdr,'red','black')}))
    col.fdr[is.na(col.fdr)] <- 'black'

    for(imethod in methods){
      dt.dcast.color[[paste0('A.Power_',imethod)]] <-
        cell_spec(x = dt.dcast.color[[paste0('A.Power_',imethod)]],color = col.power[,paste0('A.Power_',imethod)],format = format)
      dt.dcast.color[[paste0('A.Power_',imethod)]] <- gsub("NA","-",dt.dcast.color[[paste0('A.Power_',imethod)]])

      dt.dcast.color[[paste0('B.FDR_',imethod)]] <-
        cell_spec(x = dt.dcast.color[[paste0('B.FDR_',imethod)]],color = col.fdr[,paste0('B.FDR_',imethod)],format = format)
      dt.dcast.color[[paste0('B.FDR_',imethod)]] <- gsub("NA","-",dt.dcast.color[[paste0('B.FDR_',imethod)]])
    }
  }

  kb <- kbl(dt.dcast.color,
            escape = FALSE,
            format = format,
            booktabs = TRUE,
            align = c('l',rep('r',13)),
            caption = cap,
            col.names = c('Read','Samples/Group','Library Size','Read Length',rep(methods,2)),...) %>%
    add_header_above(c(" " = 4, "Power" = length(methods), "False Discovery Rate" = length(methods))) %>%
    { if(format == 'latex'){
      dotDotDot <- match.call(expand.dots = FALSE)
      if('longtable' %in% names(dotDotDot$...)){
        getLongTable <- get('longtable')
      } else{
        getLongTable <- FALSE
      }
      if(isTRUE(getLongTable)){
        opts <- c("repeat_header")
      } else{
        opts <- c("scale_down")
      }
      kable_styling(kable_input = .,latex_options = opts,font_size = font_size) %>%
        collapse_rows(kable_input = .,columns = 1, latex_hline = "major", valign = "top")
    } else{
      kable_styling(kable_input = .,bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
        collapse_rows(kable_input = .,columns = 1, valign = "top") %>%
        landscape()
    }
    }

  return(kb)
}

summarizeROCCurve <- function(x,byvar){

  cnames <- colnames(x)

  sub.cnames1 <- cnames[grepl('tpr\\.',cnames)]
  sub.cnames2 <- cnames[grepl('fdr\\.',cnames)]

  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  x.mean <- x[,lapply(.SD,mean),by = sub.byvar,.SDcols = c(sub.cnames1,sub.cnames2)]

  table1 <- melt(data = x.mean,id.vars = sub.byvar,measure.vars = sub.cnames1,variable.name = 'nFDR',value.name = 'oTPR')
  table2 <- melt(data = x.mean,id.vars = sub.byvar,measure.vars = sub.cnames2,variable.name = 'nFDR',value.name = 'oFDR')

  table1$nFDR <- as.numeric(gsub("tpr.","",table1$nFDR))/100
  table2$nFDR <- as.numeric(gsub("fdr.","",table2$nFDR))/100

  table <- merge(table1,table2,by = c(sub.byvar,'nFDR'),all.x = TRUE)

  return(table)
}

#' @importFrom data.table fwrite
summarizeQuantification <- function(path,dest,genome,fc,read,len,
                                    scenario,libs.per.group,quantifier,
                                    nsim = 20, fdr = 0.05, seq.n.gene = seq(100,3000,100),seq.n.transcript = seq(100,4500,100),seq.fdr = c(1,5,10,15,20),alpha = fdr){

  byvar <- c('Genome','Length','FC','Reads','Scenario','LibsPerGroup','Quantifier','Method','Simulation')

  res <- aggregateScenario(path = path, genome = genome, len = len, fc = fc , read = read, scenario = scenario, libs.per.group = libs.per.group, quantifier = quantifier, nsim = nsim)

  table.gene.metrics <- res$results.gene[,computeGeneMetrics(c(.BY,.SD),simulation = res$simulation.gene,features = res$features.gene,fdr = fdr,alpha = alpha),by = byvar]
  table.transcript.metrics <- res$results.transcript[,computeTranscriptMetrics(c(.BY,.SD),simulation = res$simulation.transcript,features = res$features.transcript,fdr = fdr,alpha = alpha),by = byvar]

  table.gene.fdr <- res$results.gene[,computeGeneFDRCurve(c(.BY,.SD),simulation = res$simulation.gene,features = res$features.gene,fdr = fdr,seq.n = seq.n.gene),by = byvar]
  table.transcript.fdr <- res$results.transcript[,computeTranscriptFDRCurve(c(.BY,.SD),simulation = res$simulation.transcript,features = res$features.transcript,fdr = fdr,seq.n = seq.n.transcript),by = byvar]

  table.gene.roc <- res$results.gene[,computeGeneROCCurve(c(.BY,.SD),simulation = res$simulation.gene,features = res$features.gene,fdr = fdr,seq.fdr = seq.fdr),by = byvar]
  table.transcript.roc <- res$results.transcript[,computeTranscriptROCCurve(c(.BY,.SD),simulation = res$simulation.transcript,features = res$features.transcript,fdr = fdr,seq.fdr = seq.fdr),by = byvar]

  out <- list('time' = summarizeTime(res$time,byvar),
              'metrics.gene' = summarizeMetrics(table.gene.metrics,byvar),
              'metrics.transcript' = summarizeMetrics(table.transcript.metrics,byvar),
              'fdr.gene' = summarizeFDRCurve(table.gene.fdr,byvar),
              'fdr.transcript' = summarizeFDRCurve(table.transcript.fdr,byvar),
              'roc.gene' = summarizeROCCurve(table.gene.roc,byvar),
              'roc.transcript' = summarizeROCCurve(table.transcript.roc,byvar))

  if(fc == 'fc1'){
    out[['quantile.gene']] = summarizeQQ(res$results.gene,byvar)
    out[['quantile.transcript']] = summarizeQQ(res$results.transcript,byvar)

    out[['pvalue.gene']] = summarizePValue(res$results.gene,byvar)
    out[['pvalue.transcript']] = summarizePValue(res$results.transcript,byvar)
  }

  # plotFDRCurve(out$fdr,max.n = max(seq.n))
  # plotPowerBars(out$metrics,fdr,max.n)
  # plotType1Error(out$metrics,alpha)
  # plotTime(out$time)
  # plotQQPlot(out$quantile)
  # plotPValues(out$pvalue)

  dir.create(dest,showWarnings = FALSE,recursive = TRUE)
  lapply(names(out),function(x){fwrite(x = out[[x]], file = file.path(dest,paste0(x,'.tsv.gz')),quote = FALSE,sep = '\t')})

  return(invisible())
}

summarizeScenario <- function(x,table,path,dest){
  dt <- as.character(table[x,])
  names(dt) <- colnames(table)
  table.names = c('fdr','metrics','time','quantile','pvalue','overdispersion') # Need to update this with the actual output table names

  in.path <- file.path(path,do.call(file.path,as.list(dt)))
  out.path <- file.path(dest,do.call(file.path,as.list(dt)))

  # Check is simulation directory exists
  if (!dir.exists(in.path)) return(invisible())
  message("Reading simulations from:")
  message(in.path)

  # Summarizing results with Salmon
  if (!all(file.exists(file.path(out.path,'dtu-salmon',paste0(table.names,'.tsv.gz'))))){
    # Verbose
    summarizeQuantification(path = in.path,quantifier = 'salmon',
                            dest = file.path(out.path,'dtu-salmon'),
                            genome = dt['genome'],fc = dt['fc'],read = dt['read'],
                            scenario = dt['scenario'],
                            libs.per.group = dt['libs.per.group'],len = dt['len'])
  }
}

summarizeSimulation <- function(path,
                                dest,
                                genome = c('mm39'),
                                len = c(50,75,100,125,150),
                                fc = c(1,2),
                                read = c('single-end','paired-end'),
                                scenario = c('balanced','unbalanced'),
                                libs.per.group = c(3,5,10), ...){

  path <- normalizePath(path)
  dir.create(dest,showWarnings = FALSE,recursive = TRUE)
  dest <- normalizePath(dest)

  dt.scenario <- expand.grid('genome' = genome,
                             'len' = paste0('readlen-',len),
                             'fc' = paste0('fc',fc),
                             'read' = read,
                             'scenario' = scenario,
                             'libs.per.group' = paste0(libs.per.group,'libsPerGroup'),
                             stringsAsFactors = FALSE)

  bplapply(seq_len(nrow(dt.scenario)),summarizeScenario,table = dt.scenario,dest = dest,path = path,...)

  return(invisible())
}
