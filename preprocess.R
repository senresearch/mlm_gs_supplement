library(data.table) # quickly read in tables

# Read in keys (with mutant and spatial information) for each of the 6 plates. 
KEIO_keys = lapply(1:6, function(i){
  read.csv(paste0("./processed/raw_KEIO_data/KEIO", i, "_KEY.csv"), sep="\t")
})

# Read in the Kritikos condition information. 
krit_cond = lapply(1:6, function(i){
  read.csv(paste0("./processed/raw_KEIO_data/krit_condition_files/p", i, 
                  "_4120krit.csv"), sep="\t")
})

###############################################################################

# Function to read in the iris/dat files and pull out the opacity column.
# Assumes that rows and columns in the data file are in the same order as in the mutant key. 
# dat_file is a character string naming the .iris or .dat file to read. 
# mutant_key is the mutant key data frame to use
# directory is the directory name holding the data files
# ... additional arguments to be passed into read.csv
read.dat = function(dat_file, mutant_key, directory, ...) {
  dat_temp = read.csv(paste(directory, dat_file, sep=""), sep="\t", ...)
  if(all(mutant_key$row == dat_temp$row) && 
     all(mutant_key$column == dat_temp$column)) {
    return(dat_temp$opacity)
  } else {
    stop("Check your rows and columns.")
  }
}

# Get the "Y" matrix of colony opacity for all batches. Each takes about 1-2 seconds to run.
# Krit conditions
krit_dat = lapply(1:6, function(i){
  t(sapply(krit_cond[[i]]$KRIT.FILE, read.dat, KEIO_keys[[i]], 
           "./processed/raw_KEIO_data/krit_dat/"))
})

# Plate 5: Drop rows with cond "novobiocin" and conditon "null"
krit_dat[[5]] = krit_dat[[5]][!(krit_cond[[5]]$Condition=="novobiocin" & 
                                  krit_cond[[5]]$Concentration=="null"),]
krit_cond[[5]] = krit_cond[[5]][!(krit_cond[[5]]$Condition=="novobiocin" & 
                                    krit_cond[[5]]$Concentration=="null"), ]

###############################################################################

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

for (i in 1:6) {
  krit_cond[[i]]$Concentration = 
    clean_concentrations(krit_cond[[i]]$Concentration, 
                         krit_cond[[i]]$Condition)
  krit_cond[[i]]$Cond_Conc = paste(krit_cond[[i]]$Concentration, 
                                   krit_cond[[i]]$Condition)
}

###############################################################################

krit_cond_conc_names = lapply(krit_cond, function(x){x$Cond_Conc})
krit_mut_names = lapply(KEIO_keys, function(x){as.character(x$name)})

length(do.call(c, krit_cond_conc_names))
length(unique(do.call(c, krit_cond_conc_names)))

length(do.call(c, krit_mut_names))
length(unique(do.call(c, krit_mut_names)))

sapply(krit_cond_conc_names, length)
sapply(krit_cond_conc_names, function(x){length(unique(x))})

sapply(krit_mut_names, length)
sapply(krit_mut_names, function(x){length(unique(x))})

###############################################################################

for (i in 1:6) {
  write.csv(krit_cond[[i]], 
            file=paste0("./processed/processed_KEIO_data/p", i, 
                        "_krit_cond.csv"), row.names=FALSE)
  write.csv(krit_dat[[i]], 
            file=paste0("./processed/processed_KEIO_data/p", i, 
                        "_krit_dat.csv"), row.names=FALSE)
  
  # Cond-conc
  write.table(krit_cond_conc_names[[i]], 
              file=paste0("./processed/processed_KEIO_data/p", i, 
                          "_krit_cond_conc_names.csv"), 
              sep=",", row.names=FALSE, col.names=FALSE)
  # Mutations
  write.table(krit_mut_names[[i]], 
              file=paste0("./processed/processed_KEIO_data/p", i, 
                          "_krit_mut_names.csv"), 
              sep=",", row.names=FALSE, col.names=FALSE)
}
