#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/mofapipeline
========================================================================================
 nf-core/mofapipeline Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/mofapipeline
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     nf-core/mofapipeline v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/mofapipeline --omics '*.csv'

    Mandatory arguments:
      --omics                       Path to input data (must be surrounded with quotes)

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}


//SET UP CONFIGURATION VARIABLES

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables


// Create a channel for input omics files

params.omics = "/home/gust/Desktop/project/nf-core-mofapipeline/CLL_data.RData"


// Different analysis


// Implementation of MOFA


process mofa {

  output:
  stdout mofa

    """
    #!/usr/bin/env Rscript

  rm(list=ls())
  library(MOFAtools)
  library(magrittr)

  GEdata <- read.csv("/home/gust/Desktop/project/test_MOFA/NormCountsElisabeth_OT2AVG.csv", row.names=1, comment.char="" , colClasses=c("character",rep("numeric",9)), strip.white=FALSE)
  MetabData <- read.csv("/home/gust/Desktop/project/test_MOFA/Metabolomics_transposed.csv", row.names=1, comment.char="" , colClasses=c("character",rep("numeric",9)), strip.white=FALSE)
  ProteomicsData <- read.csv("/home/gust/Desktop/project/test_MOFA/Proteomics_massaged_reordered.csv", row.names=1, comment.char="", colClasses=c("character",rep("numeric",9)), strip.white=FALSE)

  MATS <- list (GEdata, MetabData, ProteomicsData)
  names(MATS) <- c("GEdata", "MetabData", "ProteomicsData")

  #filtering out features without variance
  for (i in names(MATS)) {
    MATS[[i]] <- MATS[[i]][0<apply(MATS[[i]], 1, var),]
  }

  MATS2 <- list(GEdata=MATS\$GEdata %>% as.matrix, MetabData=MATS\$MetabData %>% as.matrix, ProteomicsData=MATS\$ProteomicsData %>% as.matrix)
  
  str(MATS2)

  MOFAobject <- createMOFAobject(MATS2)
  MOFAobject

  pdf("TilesData.pdf")
  plotTilesData(MOFAobject)
  dev.off()

  DataOptions <- getDefaultDataOptions()
  
  ModelOptions <- getDefaultModelOptions(MOFAobject)
  

  TrainOptions <- getDefaultTrainOptions()
  TrainOptions\$maxiter <- 5

  MOFAobject <- prepareMOFA(
      MOFAobject, 
      DataOptions = DataOptions,
      ModelOptions = ModelOptions,
      TrainOptions = TrainOptions
      )

    MOFAobject <- runMOFA(MOFAobject)
    
    MOFAobject
    r2 <- calculateVarianceExplained(MOFAobject)
  r2
  
  pdf("model.pdf")
  plotWeightsHeatmap(MOFAobject, "GEdata", factors=1:5, show_colnames=F)
  plotWeightsHeatmap(MOFAobject, "ProteomicsData", factors=1:5, show_colnames=F)
  plotWeightsHeatmap(MOFAobject, "MetabData", factors=1:5, show_colnames=F)
  plotWeights(MOFAobject, view = "MetabData", factor = 1)
  plotWeights(MOFAobject, view = "GEdata", factor = 1)
  plotWeights(MOFAobject, view = "ProteomicsData", factor = 1)
  plotTopWeights(MOFAobject, "ProteomicsData", 1)
  plotTopWeights(MOFAobject, "GEdata", 1)
  plotTopWeights(MOFAobject, "MetabData", 1)
  dev.off()
  """
}

mofa.subscribe {
  println it.trim()
}