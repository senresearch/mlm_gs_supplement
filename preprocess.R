# Script to format KEIO genetic screening data into a format suitable for matrix linear models. 
# Raw data is read-only

# Load libraries
library(data.table)

# Read in the keys for each of the 6 plates. These contain mutant and spatial information. 
KEIO1_key = read.csv("./processed/dataforJane/KEIO1_KEY.csv", sep="\t")
KEIO2_key = read.csv("./processed/dataforJane/KEIO2_KEY.csv", sep="\t")
KEIO3_key = read.csv("./processed/dataforJane/KEIO3_KEY.csv", sep="\t")
KEIO4_key = read.csv("./processed/dataforJane/KEIO4_KEY.csv", sep="\t")
KEIO5_key = read.csv("./processed/dataforJane/KEIO5_KEY.csv", sep="\t")
KEIO6_key = read.csv("./processed/dataforJane/KEIO6_KEY.csv", sep="\t")

# Read in the krit condition information. 
p1_krit = read.csv("./processed/dataforJane/krit_condition_files/p1_4120krit.csv", sep="\t")
p2_krit = read.csv("./processed/dataforJane/krit_condition_files/p2_4120krit.csv", sep="\t")
p3_krit = read.csv("./processed/dataforJane/krit_condition_files/p3_4120krit.csv", sep="\t")
p4_krit = read.csv("./processed/dataforJane/krit_condition_files/p4_4120krit.csv", sep="\t")
p5_krit = read.csv("./processed/dataforJane/krit_condition_files/p5_4120krit.csv", sep="\t")
p6_krit = read.csv("./processed/dataforJane/krit_condition_files/p6_4120krit.csv", sep="\t")

# Read in the new condition information. 
k1b1_new = read.csv("./processed/dataforJane/new_condition_files/k1b1_new.csv", sep="\t")
k1b4_new = read.csv("./processed/dataforJane/new_condition_files/k1b4_new.csv", sep="\t")
k2b1_new = read.csv("./processed/dataforJane/new_condition_files/k2b1_new.csv", sep="\t")
k2b4_new = read.csv("./processed/dataforJane/new_condition_files/k2b4_new.csv", sep="\t")
k3b1_new = read.csv("./processed/dataforJane/new_condition_files/k3b1_new.csv", sep="\t")
k3b4_new = read.csv("./processed/dataforJane/new_condition_files/k3b4_new.csv", sep="\t")
k4b1_new = read.csv("./processed/dataforJane/new_condition_files/k4b1_new.csv", sep="\t")
k4b4_new = read.csv("./processed/dataforJane/new_condition_files/k4b4_new.csv", sep="\t")
k5b1_new = read.csv("./processed/dataforJane/new_condition_files/k5b1_new.csv", sep="\t")
k5b4_new = read.csv("./processed/dataforJane/new_condition_files/k5b4_new.csv", sep="\t")
k6b1_new = read.csv("./processed/dataforJane/new_condition_files/k6b1_new.csv", sep="\t")
k6b4_new = read.csv("./processed/dataforJane/new_condition_files/k6b4_new.csv", sep="\t")

# Function to read in the iris/dat files and pull out the opacity column.
# Assumes that rows and columns in the data file are in the same order as in the mutant key. 
# dat_file is a character string naming the .iris or .dat file to read. 
# mutant_key is the mutant key data frame to use
# directory is the directory name holding the data files
# ... additional arguments to be passed into read.csv
read.dat = function(dat_file, mutant_key, directory, ...) {
  dat_temp = read.csv(paste("./processed/dataforJane/", directory, dat_file, sep=""), 
                    sep="\t", ...)
  if(all(mutant_key$row == dat_temp$row) && all(mutant_key$column == dat_temp$column)) {
    return(dat_temp$opacity)
  } else {
    stop("Check your rows and columns.")
  }
}

# Get the "Y" matrix of colony opacity for all batches. Each takes about 1-2 seconds to run.
# Krit conditions
dat_p1_krit = t(sapply(p1_krit$KRIT.FILE, read.dat, KEIO1_key, "krit_dat/"))
dat_p2_krit = t(sapply(p2_krit$KRIT.FILE, read.dat, KEIO2_key, "krit_dat/"))
dat_p3_krit = t(sapply(p3_krit$KRIT.FILE, read.dat, KEIO3_key, "krit_dat/"))
dat_p4_krit = t(sapply(p4_krit$KRIT.FILE, read.dat, KEIO4_key, "krit_dat/"))
dat_p5_krit = t(sapply(p5_krit$KRIT.FILE, read.dat, KEIO5_key, "krit_dat/"))
dat_p6_krit = t(sapply(p6_krit$KRIT.FILE, read.dat, KEIO6_key, "krit_dat/"))
# New conditions
dat_k1b1_new = t(sapply(k1b1_new$IRIS.0.9.4, read.dat, KEIO1_key, "new_dat/", skip=6))
dat_k1b4_new = t(sapply(k1b4_new$IRIS.0.9.4, read.dat, KEIO1_key, "new_dat/", skip=6))
dat_k2b1_new = t(sapply(k2b1_new$IRIS.0.9.4, read.dat, KEIO2_key, "new_dat/", skip=6))
dat_k2b4_new = t(sapply(k2b4_new$IRIS.0.9.4, read.dat, KEIO2_key, "new_dat/", skip=6))
dat_k3b1_new = t(sapply(k3b1_new$IRIS.0.9.4, read.dat, KEIO3_key, "new_dat/", skip=6))
dat_k3b4_new = t(sapply(k3b4_new$IRIS.0.9.4, read.dat, KEIO3_key, "new_dat/", skip=6))
dat_k4b1_new = t(sapply(k4b1_new$IRIS.0.9.4, read.dat, KEIO4_key, "new_dat/", skip=6))
dat_k4b4_new = t(sapply(k4b4_new$IRIS.0.9.4, read.dat, KEIO4_key, "new_dat/", skip=6))
dat_k5b1_new = t(sapply(k5b1_new$IRIS.0.9.4, read.dat, KEIO5_key, "new_dat/", skip=6))
dat_k5b4_new = t(sapply(k5b4_new$IRIS.0.9.4, read.dat, KEIO5_key, "new_dat/", skip=6))
dat_k6b1_new = t(sapply(k6b1_new$IRIS.0.9.4, read.dat, KEIO6_key, "new_dat/", skip=6))
dat_k6b4_new = t(sapply(k6b4_new$IRIS.0.9.4, read.dat, KEIO6_key, "new_dat/", skip=6))

# Discard plates that were flagged in the condition file. 
# For now, discard entire rows (plates) even if only some wells are flagged.
# The krit files don't have flagged plates. 
dat_k1b1_new = dat_k1b1_new[k1b1_new$FLAG=="-",]
k1b1_new = k1b1_new[k1b1_new$FLAG=="-",]
dat_k1b4_new = dat_k1b4_new[k1b4_new$FLAG=="-",]
k1b4_new = k1b4_new[k1b4_new$FLAG=="-",]

dat_k2b1_new = dat_k2b1_new[k2b1_new$FLAG=="-",]
k2b1_new = k2b1_new[k2b1_new$FLAG=="-",]
dat_k2b4_new = dat_k2b4_new[k2b4_new$FLAG=="-",]
k2b4_new = k2b4_new[k2b4_new$FLAG=="-",]

dat_k3b1_new = dat_k3b1_new[k3b1_new$FLAG=="-",]
k3b1_new = k3b1_new[k3b1_new$FLAG=="-",]
dat_k3b4_new = dat_k3b4_new[k3b4_new$FLAG=="-",]
k3b4_new = k3b4_new[k3b4_new$FLAG=="-",]

dat_k4b1_new = dat_k4b1_new[k4b1_new$FLAG=="-",]
k4b1_new = k4b1_new[k4b1_new$FLAG=="-",]
dat_k4b4_new = dat_k4b4_new[k4b4_new$FLAG=="-",]
k4b4_new = k4b4_new[k4b4_new$FLAG=="-",]

dat_k5b1_new = dat_k5b1_new[k5b1_new$FLAG=="-",]
k5b1_new = k5b1_new[k5b1_new$FLAG=="-",]
dat_k5b4_new = dat_k5b4_new[k5b4_new$FLAG=="-",]
k5b4_new = k5b4_new[k5b4_new$FLAG=="-",]

dat_k6b1_new = dat_k6b1_new[k6b1_new$FLAG=="-",]
k6b1_new = k6b1_new[k6b1_new$FLAG=="-",]
dat_k6b4_new = dat_k6b4_new[k6b4_new$FLAG=="-",]
k6b4_new = k6b4_new[k6b4_new$FLAG=="-",]

# Put all the conditions together for each plate. 
# The column names are consistent with the krit condition file column names
p1_all = rbindlist(list(p1_krit, k1b1_new, k1b4_new))
p2_all = rbindlist(list(p2_krit, k2b1_new, k2b4_new))
p3_all = rbindlist(list(p3_krit, k3b1_new, k3b4_new))
p4_all = rbindlist(list(p4_krit, k4b1_new, k4b4_new))
p5_all = rbindlist(list(p5_krit, k5b1_new, k5b4_new))
p6_all = rbindlist(list(p6_krit, k6b1_new, k6b4_new))

# Put all the opacity measurements together for each plate. 
dat_p1_all = rbind(dat_p1_krit, dat_k1b1_new, dat_k1b4_new)
dat_p2_all = rbind(dat_p2_krit, dat_k2b1_new, dat_k2b4_new)
dat_p3_all = rbind(dat_p3_krit, dat_k3b1_new, dat_k3b4_new)
dat_p4_all = rbind(dat_p4_krit, dat_k4b1_new, dat_k4b4_new)
dat_p5_all = rbind(dat_p5_krit, dat_k5b1_new, dat_k5b4_new)
dat_p6_all = rbind(dat_p6_krit, dat_k6b1_new, dat_k6b4_new)

# Plate 5: Drop rows with cond "novobiocin" and conditon "null"
dat_p5_all = dat_p5_all[!(p5_all$Condition=="novobiocin" & p5_all$Concentration=="null"),]
p5_all = p5_all[!(p5_all$Condition=="novobiocin" & p5_all$Concentration=="null"), ]

# Plate 5: Drop rows with cond "novobiocin" and conditon "null"
dat_p5_krit = dat_p5_krit[!(p5_krit$Condition=="novobiocin" & p5_krit$Concentration=="null"),]
p5_krit = p5_krit[!(p5_krit$Condition=="novobiocin" & p5_krit$Concentration=="null"), ]

# Replace "20ug/mL" concentration with "20 ug/mL" for condition "vancomycin"
# Replace multi-spaces with a single space, ";" with " +", and " sec" with "sec".
# Takes conc_col, column of concentrations, and cond_col, column of conditions. 
# Returns cleaned column of concentrations
clean_concentrations = function(conc_col, cond_col) {
  conc_col[cond_col=="vancomycin" & conc_col=="20ug/mL"] = "20 ug/mL"
  return(gsub("\\s+", " ", 
              gsub(";", " +", 
                   gsub(" sec", "sec", conc_col))))
}

p1_all$Concentration = clean_concentrations(p1_all$Concentration, p1_all$Condition)
p2_all$Concentration = clean_concentrations(p2_all$Concentration, p2_all$Condition)
p3_all$Concentration = clean_concentrations(p3_all$Concentration, p3_all$Condition)
p4_all$Concentration = clean_concentrations(p4_all$Concentration, p4_all$Condition)
p5_all$Concentration = clean_concentrations(p5_all$Concentration, p5_all$Condition)
p6_all$Concentration = clean_concentrations(p6_all$Concentration, p6_all$Condition)

p1_krit$Concentration = clean_concentrations(p1_krit$Concentration, p1_krit$Condition)
p2_krit$Concentration = clean_concentrations(p2_krit$Concentration, p2_krit$Condition)
p3_krit$Concentration = clean_concentrations(p3_krit$Concentration, p3_krit$Condition)
p4_krit$Concentration = clean_concentrations(p4_krit$Concentration, p4_krit$Condition)
p5_krit$Concentration = clean_concentrations(p5_krit$Concentration, p5_krit$Condition)
p6_krit$Concentration = clean_concentrations(p6_krit$Concentration, p6_krit$Condition)

write.csv(p1_all, file="./processed/processed_KEIO_data/p1_all.csv", row.names=FALSE)
write.csv(p2_all, file="./processed/processed_KEIO_data/p2_all.csv", row.names=FALSE)
write.csv(p3_all, file="./processed/processed_KEIO_data/p3_all.csv", row.names=FALSE)
write.csv(p4_all, file="./processed/processed_KEIO_data/p4_all.csv", row.names=FALSE)
write.csv(p5_all, file="./processed/processed_KEIO_data/p5_all.csv", row.names=FALSE)
write.csv(p6_all, file="./processed/processed_KEIO_data/p6_all.csv", row.names=FALSE)

write.csv(dat_p1_all, file="./processed/processed_KEIO_data/dat_p1_all.csv", row.names=FALSE)
write.csv(dat_p2_all, file="./processed/processed_KEIO_data/dat_p2_all.csv", row.names=FALSE)
write.csv(dat_p3_all, file="./processed/processed_KEIO_data/dat_p3_all.csv", row.names=FALSE)
write.csv(dat_p4_all, file="./processed/processed_KEIO_data/dat_p4_all.csv", row.names=FALSE)
write.csv(dat_p5_all, file="./processed/processed_KEIO_data/dat_p5_all.csv", row.names=FALSE)
write.csv(dat_p6_all, file="./processed/processed_KEIO_data/dat_p6_all.csv", row.names=FALSE)

write.csv(p1_krit, file="./processed/processed_KEIO_data/p1_krit.csv", row.names=FALSE)
write.csv(p2_krit, file="./processed/processed_KEIO_data/p2_krit.csv", row.names=FALSE)
write.csv(p3_krit, file="./processed/processed_KEIO_data/p3_krit.csv", row.names=FALSE)
write.csv(p4_krit, file="./processed/processed_KEIO_data/p4_krit.csv", row.names=FALSE)
write.csv(p5_krit, file="./processed/processed_KEIO_data/p5_krit.csv", row.names=FALSE)
write.csv(p6_krit, file="./processed/processed_KEIO_data/p6_krit.csv", row.names=FALSE)

write.csv(dat_p1_krit, file="./processed/processed_KEIO_data/dat_p1_krit.csv", row.names=FALSE)
write.csv(dat_p2_krit, file="./processed/processed_KEIO_data/dat_p2_krit.csv", row.names=FALSE)
write.csv(dat_p3_krit, file="./processed/processed_KEIO_data/dat_p3_krit.csv", row.names=FALSE)
write.csv(dat_p4_krit, file="./processed/processed_KEIO_data/dat_p4_krit.csv", row.names=FALSE)
write.csv(dat_p5_krit, file="./processed/processed_KEIO_data/dat_p5_krit.csv", row.names=FALSE)
write.csv(dat_p6_krit, file="./processed/processed_KEIO_data/dat_p6_krit.csv", row.names=FALSE)