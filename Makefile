.PHONY: all compare_times dosage sim dosage_sim 

all: 
	preprocess.R
	compare_times.jl
	dosage.jl
	auxotroph.R
	dosage.R
	sim.jl
	sim.R
	dosage_sim.jl
	dosage_sim.R

compare_times:
	preprocess.R
	compare_times.jl

dosage: 
	preprocess.R
	dosage.jl
	auxotroph.R
	dosage.R

sim: 
	preprocess.R
	sim.jl
	sim.R

dosage_sim: 
	preprocess.R
	dosage_sim.jl
	dosage_sim.R