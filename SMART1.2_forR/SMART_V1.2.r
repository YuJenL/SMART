
cat("Welcome to SMART software (version 1.2): \n")
cat("Please provide the full path to the SMART programs:\n")
wd <- scan(n=1, what="character")
setwd(wd)
source("SMART_V1.2_Gui.r")
