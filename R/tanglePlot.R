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

#' Plot the different nucleotides present in an alignment of nucleotide sequences
#'
#' Function uses coloured polygons to illustrate which nucleotides are present in each sequence at each position
#' @param tree1 an object of class "phylo" to compared to \code{tree2}
#' @param tree2 an object of class "phylo" to compared to \code{tree1}
#' @param onlyIfDifferent a boolean variable, if \code{TRUE} only plots lines between tips that are in different locations on the phylogenies, otherwise lines between all tips are plotted. Defaults to \code{TRUE}
#' @param connectingLine.col colour of the lines connecting tips. Defaults to black
#' @param connectingLine.lwd width of the lines connecting tips. Defaults to 1
#' @param connectingLine.lty type of the lines connecting tips. Defaults to 1 (solid)
#' @param offsetProp offset value to displace starts of connecting lines from tips, set as the proportion of plotting window. Defaults to \code{NULL} (off)
#' @param leftMargins margin parameters for the left tree plot. Defaults to \code{c(0,0,0,1)} to leave small space on right
#' @param rightMargins margin parameters for the right tree plot. Defaults to \code{c(0,1,0,0)} to leave small space on left
#' @keywords pylogeny tree tangle plot compare
#' @export
#' @examples
#' # Set the seed
#' set.seed(76263)
#'
#' # Create a random phylogeny
#' tree <- ape::rtree(n=50, rooted=TRUE)
#'
#' # Rotate one of the nodes in the phylogeny - with tip labels it
#' rotated <- ape::rotate(tree, node=c("t35", "t34"))
#'
#' # Compare the random and rotated phylogeny
#' tanglePlot(tree, rotated, connectingLine.col="red",
#'            connectingLine.lty=2, show.tip.label=TRUE, offsetProp=0.02)
tanglePlot <- function(tree1, tree2, onlyIfDifferent=TRUE,
                       connectingLine.col="black",
                       connectingLine.lwd=1,
                       connectingLine.lty=1,
                       offsetProp=NULL,
                       leftMargins=c(0,0,0,1),
                       rightMargins=c(0,1,0,0), ...){

  # Get the current plotting margins
  currentMar <- par()$mar

  # Set the number of plots in window
  currentMfrow <- par()$mfrow
  par(mfrow=c(1,2))

  # Set the plotting margin - leave space on right
  par(mar=leftMargins)

  # Plot phylogeny 1
  ape::plot.phylo(tree1, direction="rightwards", ...)
  tree1TipCoordinates <- getTipCoordinates(tree1$tip.label)

  # Set the plotting margin - leave space on right
  par(mar=rightMargins)

  # Plot phylogeny 2
  ape::plot.phylo(tree2, direction="leftwards", ...)
  tree2TipCoordinates <- getTipCoordinates(tree2$tip.label)

  # Plot lines between the phylogenies
  plotLinesBetweenAllTips(tree1TipCoordinates, tree2TipCoordinates,
                          col=connectingLine.col,
                          lwd=connectingLine.lwd,
                          lty=connectingLine.lty,
                          onlyIfDifferent=onlyIfDifferent,
                          offsetProp=offsetProp)

  # Reset the plotting margins
  par(mar=currentMar)

  # Reset the number of plots in window
  par(mfrow=currentMfrow)
}

#' Plot lines between coordinates of tips on phylogenies
#'
#' Function used by \code{tanglePlot()}
#' @param tipLocationsA a list of coordinates of tips on a phylogenetic tree plot
#' @param tipLocationsB a list of coordinates of tips on a phylogenetic tree plot
#' @param col colour of the lines connecting tips. Defaults to black
#' @param lwd width of the lines connecting tips. Defaults to 1
#' @param lty type of the lines connecting tips. Defaults to 1 (solid)
#' @param onlyIfDifferent a boolean variable, if \code{TRUE} only plots lines between tips that are in different locations on the phylogenies, otherwise lines between all tips are plotted. Defaults to \code{TRUE}
#' @param offsetProp offset value to displace starts of connecting lines from tips, set as the proportion of plotting window. Defaults to \code{NULL} (off)
#' @keywords internal
plotLinesBetweenAllTips <- function(tipLocationsA, tipLocationsB,
                                    col="black", lwd=1, lty=1,
                                    onlyIfDifferent=TRUE, offsetProp=NULL){

  # Open up the port to allow lines to be added
  grid::pushViewport(grid::viewport())
  grid::popViewport()

  # Calculate the offset for the X positions
  xOffset <- 0
  if(is.null(offsetProp) == FALSE){

    xOffset <- offsetProp * (tipLocationsA$height + tipLocationsB$height)
  }

  # Examines each of the tips on the left phylogeny
  for(key in names(tipLocationsA)){

    # Skip the tips that aren't in the right phylogeny or the absolute phylogeny height
    if(is.null(tipLocationsB[[key]]) | key == "height"){
      next
    }

    # Check if the tip location (Y coordinate) is different in the two phylogenies (using relative coordinates)
    if(abs(tipLocationsA[[key]][4] - tipLocationsB[[key]][4]) > 1){

      # Plot a connecting line between tip on left and same tip on right
      plotLineBetweenPairOfTips(coordsA=tipLocationsA[[key]], coordsB=tipLocationsB[[key]],
                                col=col, lty=lty, lwd=lwd, xOffset=xOffset)

    # Check if want to plot lines between tips regardless of if different
    }else if(onlyIfDifferent == FALSE){

      # Plot a connecting line between tip on left and same tip on right
      plotLineBetweenPairOfTips(coordsA=tipLocationsA[[key]], coordsB=tipLocationsB[[key]],
                                col=col, lty=lty, lwd=lwd, xOffset=xOffset)
    }
  }
}

#' Plot lines between sets of coordinates in the plotting window
#'
#' Function used by \code{tanglePlot()}
#' @param coordsA a list of coordinates on an existing plot
#' @param coordsB a list of coordinates on an existing plot
#' @param col colour of the lines connecting tips. Defaults to black
#' @param lwd width of the lines connecting tips. Defaults to 1
#' @param lty type of the lines connecting tips. Defaults to 1 (solid)
#' @param xOffset The distance on X axis to displace the start and ends of connecting lines. Defaults to 0
#' @keywords internal
plotLineBetweenPairOfTips <- function(coordsA, coordsB,
                                      col="black", lwd=1, lty=1, xOffset=0){

  # Open up the port to plot the current line
  grid::pushViewport(grid::viewport())

  # Plot a connecting line between tip on left and same tip on right
  grid::grid.lines(x = c(coordsA[1]+xOffset, coordsB[1]-xOffset),
                   y = c(coordsA[2], coordsB[2]),
                   gp = grid::gpar(col=col, lty=lty, lwd=lwd))

  # Close the port
  grid::popViewport()
}

#' Get the relative and absolute coordinates of tips on a plotted phylogeny
#'
#' Function used by \code{tanglePlot()}
#' @param tipLabels a vector of the tip labels of a plotted phylogeny, which can be retrieved from a "phylo" object using: \code{tree$tip.label}
#' @keywords internal
getTipCoordinates <- function(tipLabels){

  # Get the tip coordinates in the last phylogeny plot
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

  # Initialise a list to store the absolute coordinates on the plotting window
  tips <- list()

  # Store the absolute tree height
  tips[["height"]] <- grconvertX(lastPP$x.lim[2], "user", "ndc")

  # Examine each of the tips
  for(i in 1:length(tipLabels)){

    # Store the aboslute and relative coordinates for the current tip
    tips[[as.character(tipLabels[i])]] <- c(grconvertX(lastPP$xx[i], "user", "ndc"),
                                            grconvertY(lastPP$yy[i], "user", "ndc"),
                                            lastPP$xx[i], lastPP$yy[i])
  }

  return(tips)
}
