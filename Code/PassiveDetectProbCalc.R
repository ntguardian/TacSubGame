#!/usr/bin/Rscript
# PassiveDetectProbCalc.R
# 2023-10-05
# curtis
# Compute probabilities for passive sonar detection in the tactical sub game

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

# LIBRARIES --------------------------------------------------------------------

suppressPackageStartupMessages(library(discreteRV))

# EXECUTABLE SCRIPT MAIN FUNCTIONALITY -----------------------------------------

main <- function(threshold,
                 submod = 0,
                 dpmod = c(0, 0)) {
  sonar_die <- SofIID(RV(1:6), 2)
  dpmod <- as.numeric(strsplit(dpmod, ",")[[1]])
  if (length(dpmod) >= 1) {
    val <- prod(sapply(dpmod, FUN = function(dp) {
        P(sonar_die < threshold - dp)})) * P(sonar_die < threshold - submod)
  } else {
    val <- P(sonar_die < threshold - submod)
  }
  cat(1 - val,  "\n", sep = "")
}

# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION -------------------------

if (sys.nframe() == 0) {
  p <- arg_parser("Compute the probability a ship is alerted by the submarine")
  p <- add_argument(p, "threshold", type = "integer", nargs = 1,
                    help = "The threshold for detection")
  p <- add_argument(p, "--submod", type = "integer", default = 0,
                    nargs = 1, help = "Modifier for the submarine")
  p <- add_argument(p, "--dpmod", type = "character", default = "0,0",
                    nargs = Inf,
                    help = "Modifiers for individual detection points, as a comma-separated list")

  cl_args <- parse_args(p)
  cl_args <- cl_args[!(names(cl_args) %in% c("help", "opts"))]
  if (any(sapply(cl_args, is.na))) {
    # User did not specify all inputs; print help message
    print(p)
    cat("\n\nNot all needed inputs were given.\n")
    quit()
  }

  do.call(main, cl_args[2:length(cl_args)])
}

