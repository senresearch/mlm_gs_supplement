# Matrix linear models for high-throughput chemical genetic screens

This repository contains code and data that can be used to reproduce
the results presented in the paper ["Matrix linear models for high-throughput chemical genetic screens"](http://dx.doi.org/10.1534/genetics.119.302299).

- Implemented in [Julia](https://julialang.org/downloads/) [^fn1] and [R](https://cran.r-project.org/mirrors.html) [^fn6]. 

- Julia functions used to run matrix linear models for genetic screening data are provided in the [GeneticScreens](https://github.com/senresearch/GeneticScreens.jl) package, which is an extension of the [MatrixLM](https://github.com/senresearch/MatrixLM.jl) package. 

- Genetic screening data from Nichols et al. [^fn5] used for analysis is available [here](https://figshare.com/s/f7da693dee83595eafd7). Once downloaded, it should be saved in the `data/raw_KEIO_data/` directory. 

- Running `auxotroph.R` additionally requires downloading the following tables and saving them as CSVs in the `data/` directory. 
    - [Supplemental Table 4](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3060659/bin/NIHMS261392-supplement-04.xls) in Nichols et al. [^fn5]
    - [Supplemental Table 1](http://systemsbiology.ucsd.edu/publications/supplemental_material/JBact2006/) in Joyce et al. [^fn3]
    - [Supplemental Table 3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5464397/bin/NIHMS72934-supplement-Supplementary_table_3.zip) in Kritikos et al. [^fn4]

---

- `code/preprocess.R`: Preprocess Nichols et al.'s data [^fn5]. Requires data contained in `data/raw_KEIO_data/` and should be run before any of the other files. 

---

- `code/compare_times.jl`: Compare the runtimes for matrix linear models and Collins et al.'s S scores [^fn2]. 

---

- `code/dosage.jl`: Run matrix linear models (dosage-response and condition-concentrations) and Collins et al.'s S scores [^fn2] (condition-concentrations) on Nichols et al.'s data [^fn5]. 

- `auxotroph.R`: Reproduce plots to check for auxotrophs against the lists provided by Supplemental Table 4 in Nichols et al. [^fn5] and Supplemental Table 1 in Joyce et al. [^fn3], as well as analysis of Kritikos et al.'s S scores (Supplemental Table 3) [^fn4]. 

- `code/dosage.R`: Reproduce plots of proportion of hits detected by dosage response approach compared to matrix linear models and Collins et al.'s S scores [^fn2]. 

---

- `code/sim.jl`: Run matrix linear models and Collins et al.'s S scores [^fn2] on simulated data. 

- `code/sim.R`: Reproduce ROC plots for comparing matrix linear models and Collins et al.'s S scores [^fn2]. 

- `code/sim_null.jl`: Permute simulated data to calculate Type I error for matrix linear models. 

---

- `code/dosage_sim.jl`: Run matrix linear models (dosage-response) and Collins et al.'s S scores [^fn2] (condition-concentrations and conditions only) on simulated data. 

- `code/dosage_sim.R`: Reproduce ROC plots for comparing matrix linear models (dosage-response) and Collins et al.'s S scores [^fn2] (condition-concentrations and conditions only). 


[^fn1]: Bezanson, J., Edelman, A., Karpinski, S., and Shah, V. B. (2017). Julia: A fresh approach to numerical computing. *SIAM review*, 59(1):65–98.

[^fn2]: Collins, S. R., Schuldiner, M., Krogan, N. J., and Weissman, J. S. (2006). A strategy for extracting and analyzing large-scale quantitative epistatic interaction data. *Genome biology*, 7(7):R63. 

[^fn3]: Joyce, A. R., Reed, J. L., White, A., Edwards, R., Osterman, A., Baba, T., Mori, H., Lesely, S. A., Palsson, B. Ø., and Agarwalla, S. (2006). Experimental and computational assessment of conditionally essential genes in Escherichia coli. *Journal of bacteriology*, 188(23):8259–8271. 
    
[^fn4]: Kritikos, G., Banzhaf, M., Herrera-Dominguez, L., Koumoutsi, A., Wartel, M., Zietek, M., and Typas, A. (2017). A tool named iris for versatile high-throughput phenotyping in microorganisms. *Nature microbiology*, 2(5):17014.

[^fn5]: Nichols, R. J., Sen, S., Choo, Y. J., Beltrao, P., Zietek, M., Chaba, R., Lee, S., Kazmierczak, K. M., Lee, K. J., Wong, A., et al. (2011). Phenotypic landscape of a bacterial cell. *Cell*, 144(1):143–156. 

[^fn6]: R Core Team (2018). *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienna, Austria.
