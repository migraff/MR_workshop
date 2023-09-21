
###--- install Mendelian randomization package ---###

# install.packages("MendelianRandomization", dependencies = TRUE)
library(MendelianRandomization)

###--- install Meta-Analysis Package ---###

# install.packages("metafor", dependencies = TRUE)
library(metafor)

###--- install MRbase items ---###

# install.packages("devtools")
library(devtools)

# install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)

###--- install MRPRESSO ---###

# if (!require("devtools")) { install.packages("devtools") } else {}

# devtools::install_github("rondolab/MR-PRESSO")

library(MRPRESSO)

library(tidyverse)