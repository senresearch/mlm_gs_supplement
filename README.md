Repository for storing code and data to include in the supplement for the 
genetic screening paper.

---

- `preprocess.R`: Preprocess Nichols [^fn3] data

---

- `compare_times.jl`: Compares the runtimes for MLM and Collins S scores 
[^fn1]

---

- `dosage.jl`: Run MLM (dosage-response and condition-concentrations) and 
Collins S scores [^fn1] on Nichols [^fn3] data

- `auxotroph.R`: Reproduce plots for auxotroph checks against Nichols [^fn3] 
and Joyce [^fn2]

- `dosage.R`: Reproduce plots of proportion of hits detected by dosage 
response approach

---

- `sim.jl`: Run MLM and Collins S scores [^fn1] on simulated data

- `sim.R`: Reproduce ROC plots for comparing MLM and Collins S scores [^fn1]

---

- `dosage_sim.jl`: Run MLM 
(dosage-response, condition-concentrations, and conditions only) on 
simulated data

- `dosage_sim.R`: Reproduce ROC plots for comparing MLM encoding approaches
(dosage-response, condition-concentrations, and conditions only)


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
