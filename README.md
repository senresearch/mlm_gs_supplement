Supplemental code for ["Matrix linear models for high-throughput chemical genetic screens"](https://www.biorxiv.org/content/10.1101/468140v1).

- Implemented in [Julia](https://julialang.org/downloads/) [^fn2] and [R](https://cran.r-project.org/mirrors.html) [^fn9]. 

- Julia functions used to run matrix linear models for genetic screening data are provided in the [GeneticScreen](https://github.com/janewliang/GeneticScreen.jl) package, which is an extension of the [matrixLM](https://github.com/janewliang/matrixLM.jl) package. 

- Genetic screening data from Nichols et al. [^fn8] used for analysis is available upon request. Once downloaded, it should be saved in the `data/raw_KEIO_data/` directory. 

- Running `auxotroph.R` additionally requires downloading the following tables and saving them as CSVs in the `data/` directory. 
    - [Supplemental Table 4](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3060659/bin/NIHMS261392-supplement-04.xls) in Nichols et al. [^fn8]
    - [Supplemental Table 1](http://systemsbiology.ucsd.edu/publications/supplemental_material/JBact2006/) in Joyce et al. [^fn6]
    - [Supplemental Table 3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5464397/bin/NIHMS72934-supplement-Supplementary_table_3.zip) in Kritikos et al. [^fn7]

---

- `code/Makefile`: Recipes for running each of the following components of the analysis in serial; all data dependencies should be downloaded and placed in the `data/` directory prior to running any recipes

---

- `code/preprocess.R` [^fn4]: Preprocess Nichols et al.'s data [^fn8]; requires data contained in `data/raw_KEIO_data/` and should be run before any of the other files

---

- `code/compare_times.jl`: Compare the runtimes for matrix linear models and Collins et al.'s S scores [^fn3]

---

- `code/dosage.jl`: Run matrix linear models (dosage-response and condition-concentrations) and Collins et al.'s S scores [^fn3] (condition-concentrations) on Nichols et al.'s data [^fn8]

- `auxotroph.R` [^fn11] [^fn5]: Reproduce plots to check for auxotrophs against the lists provided by Supplemental Table 4 in Nichols et al. [^fn8] and Supplemental Table 1 in Joyce et al. [^fn6], as well as analysis of Kritikos et al.'s S scores (Supplemental Table 3) [^fn7]

- `code/dosage.R` [^fn10] [^fn1]: Reproduce plots of proportion of hits detected by dosage response approach compared to matrix linear models and Collins et al.'s S scores [^fn3]

---

- `code/sim.jl`: Run matrix linear models and Collins et al.'s S scores [^fn3] on simulated data

- `code/sim.R` [^fn10] [^fn1] [^fn5]: Reproduce ROC plots for comparing matrix linear models and Collins et al.'s S scores [^fn3]

---

- `code/dosage_sim.jl`: Run matrix linear models (dosage-response) and Collins et al.'s S scores [^fn3] (condition-concentrations and conditions only) on simulated data

- `code/dosage_sim.R` [^fn10] [^fn1] [^fn5]: Reproduce ROC plots for comparing matrix linear models (dosage-response) and Collins et al.'s S scores [^fn3] (condition-concentrations and conditions only)


[^fn1]: Benjamini, Y. and Hochberg, Y. (2000). On the adaptive control of the false discovery rate in multiple testing with independent statistics. *Journal of educational and Behavioral Statistics*, 25(1):60–83.

[^fn2]: Bezanson, J., Edelman, A., Karpinski, S., and Shah, V. B. (2017). Julia: A fresh approach to numerical computing. *SIAM review*, 59(1):65–98.

[^fn3]: Collins, S. R., Schuldiner, M., Krogan, N. J., and Weissman, J. S. (2006). A strategy for extracting and analyzing large-scale quantitative epistatic interaction data. *Genome biology*, 7(7):R63. 

[^fn4]: Dowle, M. and Srinivasan, A. (2018). *data.table: Extension of ‘data.frame‘*. R package version 1.11.8.

[^fn5]: Ekstrøm, C. T. (2018). *MESS: Miscellaneous Esoteric Statistical Scripts*. R package version 0.5.2.

[^fn6]: Joyce, A. R., Reed, J. L., White, A., Edwards, R., Osterman, A., Baba, T., Mori, H., Lesely, S. A., Palsson, B. Ø., and Agarwalla, S. (2006). Experimental and computational assessment of conditionally essential genes in Escherichia coli. *Journal of bacteriology*, 188(23):8259–8271. 
    
[^fn7]: Kritikos, G., Banzhaf, M., Herrera-Dominguez, L., Koumoutsi, A., Wartel, M., Zietek, M., and Typas, A. (2017). A tool named iris for versatile high-throughput phenotyping in microorganisms. *Nature microbiology*, 2(5):17014.

[^fn8]: Nichols, R. J., Sen, S., Choo, Y. J., Beltrao, P., Zietek, M., Chaba, R., Lee, S., Kazmierczak, K. M., Lee, K. J., Wong, A., et al. (2011). Phenotypic landscape of a bacterial cell. *Cell*, 144(1):143–156. 

[^fn9]: R Core Team (2018). *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienna, Austria.

[^fn10]: Team, M. C., Blanchard, G., Dickhaus, T., Hack, N., Konietschke, F., Rohmeyer, K., Rosenblatt, J., Scheer, M., and Werft, W. (2017). *mutoss: Unified Multiple Testing Procedures*. R package version 0.1-12.

[^fn11]: Wickham, H. (2007). Reshaping data with the reshape package. *Journal of Statistical Software*, 21(12):1–20.
