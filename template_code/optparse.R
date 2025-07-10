#!/usr/bin/env Rscript

# Sets my "packages" directory 
.libPaths("~/packages/")

library(optparse, quietly = T)

###############
# OPTIONS.....#
###############
# These are options for a chip-seq pipeline.
# Change as needed

option_list = list(
  make_option(c("-d", "--dir"), type="character", default=".", 
              help="directory", metavar="character"),
  make_option(c("-n", "--number"), type="double", default=1000, 
              help="number of peaks to select [default= %default]", metavar="double")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check null options
if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

###############
# PROGRAM.....#
###############
