.PHONY: all compare_times dosage sim dosage_sim 

all: 
	R CMD BATCH preprocess.R
	julia compare_times.jl
	julia dosage.jl
	R CMD BATCH auxotroph.R
	R CMD BATCH dosage.R
	julia sim.jl
	R CMD BATCH sim.R
	julia sim_null.jl
	julia dosage_sim.jl
	R CMD BATCH dosage_sim.R

compare_times:
	R CMD BATCH preprocess.R
	compare_times.jl

dosage: 
	R CMD BATCH preprocess.R
	julia dosage.jl
	R CMD BATCH auxotroph.R
	R CMD BATCH dosage.R

sim: 
	R CMD BATCH preprocess.R
	julia sim.jl
	R CMD BATCH sim.R

dosage_sim: 
	R CMD BATCH preprocess.R
	julia dosage_sim.jl
	R CMD BATCH dosage_sim.R