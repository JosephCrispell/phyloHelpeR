# Tutorials
#https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
#http://r-pkgs.had.co.nz/description.html
#https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html

## Packages to install
#install.packages("devtools")
#install.packages("digest")
#devtools::install_github("klutometis/roxygen")

## Packages to load
#library("devtools")
#library("roxygen2")

## Creating package
#packageDirectory <- "/home/josephcrispell/phyloHelpeR/"
#usethis::create_package(packageDirectory)

## Documenting changes
#setwd(packageDirectory)
#document()

## Install
#setwd("..")
#install("phyloHelpeR")

#' Run command line RAxML tool from within R with bootstrapping
#'
#' Runs RAxML from within R with bootstrapping. Note that it requires raxmlHPC to available in the path.
#' @param fastaFile A character string containing the full path to the FASTA file
#' @param date A character string corresponding to date e.g. 06-12-19 for 6th December 2019
#' @param nBootstraps The number of bootstraps for RAxML to run. Defaults to 100.
#' @param nThreads The number of threads for RAxML to use. Defaults to 6.
#' @param outgroup A character string identifying a tip in the phylogeny to use as an outgroup. Defaults to \code{NULL}
#' @param model A character string identifying the substitution model to used by RAxML
#' @keywords phylogeny raxml maximum likelihood
#' @export
runRAXML <- function(fastaFile, date, path, nBootstraps=100, nThreads=6, outgroup=NULL, model="GTRCAT"){

  # Create a directory for the output file
  directory <- file.path(path, paste0("RAxML_", date))

  print(directory)

  # Check if analyses already run
  alreadyRun <- dir.exists(directory)

  # If not already run, create output directory for RAxML
  if(alreadyRun == FALSE){
    suppressWarnings(dir.create(directory))
  }

  # Get and set the Working directory - this will be where the output files are dumped
  currentDirectory <- getwd()
  setwd(directory)

  # Build analysis name
  analysisName <- paste("RaxML-R_", date, sep="")

  # Check if already Run and just want to retrieve tree
  if(alreadyRun == FALSE){

    # Build the command
    seeds <- sample(1:100000000, size=2, replace=FALSE) # For parsimony tree and boostrapping

    if(is.null(outgroup)){
      command <- paste("raxmlHPC",
                       " -f a", # Algorithm: Rapid boostrap inference
                       " -N ", nBootstraps,
                       " -T ", nThreads,
                       " -m ", model, " -V", # -V means no rate heterogenity
                       " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                       " -n ", analysisName,
                       " -s ", fastaFile, sep="")
    }else{
      command <- paste("raxmlHPC",
                       " -f a", # Algorithm: Rapid boostrap inference
                       " -N ", nBootstraps,
                       " -T ", nThreads,
                       " -m ", model, " -V", # -V means no rate heterogenity
                       " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                       " -n ", analysisName,
                       " -s ", fastaFile,
                       " -o ", outgroup, sep="")
    }

    system(command, intern=TRUE)
  }

  # Get the tree and read it in
  treeBS <- getTreeFileWithSupportValues(analysisName)

  # Reset working directory
  setwd(currentDirectory)

  return(treeBS)
}

#' Retrieves the tree file with bootstrap values from RAxML output
#'
#' Function used by \code{runRAXML()}
#' @param analysisName A character string indicating the name of the analysis ran in RAXML
#' @keywords internal
#' @return Returns a phylo tree object
getTreeFileWithSupportValues <- function(analysisName){

  # Get files in current working directory
  files <- list.files()

  # Select the tree file with BS support values
  treeBSFile <- files[grepl(files, pattern=paste("RAxML_bipartitions[.]", analysisName, sep="")) == TRUE]

  # Open the file
  treeBS <- ape::read.tree(treeBSFile)

  return(treeBS)
}
