# Load in raw data from .csv file
SampleRamanCurves <- read.csv("data-raw/SampleRamanCurves.csv")

# No preprocessing. Just saving as an .rda file in the data dir.
usethis::use_data(SampleRamanCurves)
