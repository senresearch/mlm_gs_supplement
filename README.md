Repository for storing code and data to include in the supplement for the 
genetic screening paper.

---

- `preprocess.R`: Preprocess Nichols's [^fn3] data

---

- `compare_times.jl`: Compares the runtimes for matrix linear models and 
Collins's S scores [^fn1]

---

- `dosage.jl`: Run matrix linear models 
(dosage-response and condition-concentrations) and Collins's S scores 
(condition-concentrations) [^fn1] on Nichols's [^fn3] data

- `auxotroph.R`: Reproduce plots to check for auxotrophs against the lists 
provided in Supplemental Table 4 in Nichols [^fn3] and Supplemental Table 1 
in Joyce [^fn2]

- `dosage.R`: Reproduce plots of proportion of hits detected by dosage 
response approach compared to matrix linear models and Collins's S scores

---

- `sim.jl`: Run matrix linear models and Collins's S scores [^fn1] on 
simulated data

- `sim.R`: Reproduce ROC plots for comparing matrix linear models and 
Collins's S scores [^fn1]

---

- `dosage_sim.jl`: Run matrix linear models (dosage-response) and Collins's S 
scores (condition-concentrations and conditions only) on simulated data

- `dosage_sim.R`: Reproduce ROC plots for comparing matrix linear models 
(dosage-response) and Collins's S scores 
(condition-concentrations and conditions only)


[^fn1]: Collins, S. R., Schuldiner, M., Krogan, N. J., and Weissman, J. S. 
    (2006). A strategy for extracting and analyzing large-scale quantitative 
    epistatic interaction data. Genome biology, 7(7):R63. 

[^fn2]: Joyce, A. R., Reed, J. L., White, A., Edwards, R., Osterman, A., 
    Baba, T., Mori, H., Lesely, S. A., Palsson, B. Ø., and Agarwalla, S. 
    (2006). Experimental and computational assessment of conditionally 
    essential genes in Escherichia coli. Journal of bacteriology, 
    188(23):8259–8271. 

[^fn3]: Nichols, R. J., Sen, S., Choo, Y. J., Beltrao, P., Zietek, M., 
    Chaba, R., Lee, S., Kazmierczak, K. M., Lee, K. J., Wong, A., et al. 
    (2011). Phenotypic landscape of a bacterial cell. Cell, 144(1):143–156. 
