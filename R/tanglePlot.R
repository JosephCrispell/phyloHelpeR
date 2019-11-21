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



#### Preparation ####

# Load libraries
library(ape)
library(grid)
library(phytools)

# Set the seed
set.seed(76263)

#### Create a random phylogeny and rotate a node ####

# Create a random phylogeny
tree <- rtree(n=50, rooted=TRUE)

# Rotate one of the nodes in the phylogeny
rotated <- rotateNodes(tree, nodes=c(56))

#### Compare the random and rotated phylogeny ####

tanglePlot(tree, rotated, connectingLine.col="red",
           connectingLine.lty=2, show.tip.label=FALSE)

#### FUNCTIONS ####

tanglePlot <- function(tree1, tree2, onlyIfDifferent=TRUE,
                       connectingLine.col="black",
                       connectingLine.lwd=1,
                       connectingLine.lty=1,
                       ...){

  # Get the current plotting margins
  currentMar <- par()$mar

  # Set the number of plots in window
  currentMfrow <- par()$mfrow
  par(mfrow=c(1,2))

  # Set the plotting margin - leave space on right
  par(mar=c(0,0,0,1))

  # Plot phylogeny 1
  plot.phylo(tree1, direction="rightwards", ...)
  tree1TipCoordinates <- getTipCoordinates(tree1$tip.label)

  # Set the plotting margin - leave space on right
  par(mar=c(1,0,0,0))

  # Plot phylogeny 2
  plot.phylo(tree2, direction="leftwards", ...)
  tree2TipCoordinates <- getTipCoordinates(tree2$tip.label)

  # Plot lines between the phylogenies
  plotLinesBetweenTips(tree1TipCoordinates, tree2TipCoordinates,
                       col=connectingLine.col,
                       lwd=connectingLine.lwd,
                       lty=connectingLine.lty,
                       onlyIfDifferent=onlyIfDifferent)

  # Reset the plotting margins
  par(mar=currentMar)

  # Reset the number of plots in window
  par(mfrow=currentMfrow)
}

plotLinesBetweenTips <- function(tipLocationsA, tipLocationsB,
                                 col="black", lwd=1, lty=1,
                                 onlyIfDifferent=TRUE){

  pushViewport(viewport())
  popViewport()

  for(key in names(tipLocationsA)){

    if(is.null(tipLocationsB[[key]])){
      next
    }

    if(abs(tipLocationsA[[key]][4] - tipLocationsB[[key]][4]) > 1){

      pushViewport(viewport())
      grid.lines(x = c(tipLocationsA[[key]][1], tipLocationsB[[key]][1]),
                 y = c(tipLocationsA[[key]][2], tipLocationsB[[key]][2]),
                 gp = gpar(col=col, lty=lty, lwd=lwd))
      popViewport()
    }else if(onlyIfDifferent == FALSE){
      pushViewport(viewport())
      grid.lines(x = c(tipLocationsA[[key]][1], tipLocationsB[[key]][1]),
                 y = c(tipLocationsA[[key]][2], tipLocationsB[[key]][2]),
                 gp = gpar(col=col, lty=lty, lwd=lwd))
      popViewport()
    }
  }
}

getTipCoordinates <- function(tipLabels){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tips <- list()
  for(i in 1:length(tipLabels)){
    tips[[as.character(tipLabels[i])]] <- c(grconvertX(lastPP$xx[i], "user", "ndc"),
                                            grconvertY(lastPP$yy[i], "user", "ndc"),
                                            lastPP$xx[i], lastPP$yy[i])
  }

  return(tips)
}
