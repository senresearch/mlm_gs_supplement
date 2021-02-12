# Matrix linear models for high-throughput chemical genetic screens

This repository contains code and data that can be used to reproduce
the results presented in the paper ["Matrix linear models for high-throughput chemical genetic screens"](http://dx.doi.org/10.1534/genetics.119.302299).

- Implemented in [Julia](https://julialang.org/downloads/) <sup>[1](#myfootnote1)</sup> and [R](https://cran.r-project.org/mirrors.html) <sup>[6](#myfootnote6)</sup>. 

- Julia functions used to run matrix linear models for genetic screening data are provided in the [GeneticScreens](https://github.com/senresearch/GeneticScreens.jl) package, which is an extension of the [MatrixLM](https://github.com/senresearch/MatrixLM.jl) package. 

- Genetic screening data from Nichols et al. <sup>[5](#myfootnote5)</sup> used for analysis is available [here](https://figshare.com/s/f7da693dee83595eafd7). Once downloaded, it should be saved in the `data/raw_KEIO_data/` directory. 

- Running `auxotroph.R` additionally requires downloading the following tables and saving them as CSVs in the `data/` directory. 
    - [Supplemental Table 4](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3060659/bin/NIHMS261392-supplement-04.xls) in Nichols et al. <sup>[5](#myfootnote5)</sup>
    - [Supplemental Table 1](http://systemsbiology.ucsd.edu/publications/supplemental_material/JBact2006/) in Joyce et al. <sup>[3](#myfootnote3)</sup>
    - [Supplemental Table 3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5464397/bin/NIHMS72934-supplement-Supplementary_table_3.zip) in Kritikos et al. <sup>[4](#myfootnote4)</sup>

---

- `code/preprocess.R`: Preprocess Nichols et al.'s data <sup>[5](#myfootnote5)</sup>. Requires data contained in `data/raw_KEIO_data/` and should be run before any of the other files. 

---

- `code/compare_times.jl`: Compare the runtimes for matrix linear models and Collins et al.'s S scores <sup>[2](#myfootnote2)</sup>. 

---

- `code/dosage.jl`: Run matrix linear models (dosage-response and condition-concentrations) and Collins et al.'s S scores <sup>[2](#myfootnote2)</sup> (condition-concentrations) on Nichols et al.'s data <sup>[5](#myfootnote5)</sup>. 

- `auxotroph.R`: Reproduce plots to check for auxotrophs against the lists provided by Supplemental Table 4 in Nichols et al. <sup>[5](#myfootnote5)</sup> and Supplemental Table 1 in Joyce et al. <sup>[3](#myfootnote3)</sup>, as well as analysis of Kritikos et al.'s S scores (Supplemental Table 3) <sup>[4](#myfootnote4)</sup>. 

- `code/dosage.R`: Reproduce plots of proportion of hits detected by dosage response approach compared to matrix linear models and Collins et al.'s S scores <sup>[2](#myfootnote2)</sup>. 

---

- `code/sim.jl`: Run matrix linear models and Collins et al.'s S scores <sup>[2](#myfootnote2)</sup> on simulated data. 

- `code/sim.R`: Reproduce ROC plots for comparing matrix linear models and Collins et al.'s S scores <sup>[2](#myfootnote2)</sup>. 

- `code/sim_null.jl`: Permute simulated data to calculate Type I error for matrix linear models. 

---

- `code/dosage_sim.jl`: Run matrix linear models (dosage-response) and Collins et al.'s S scores <sup>[2](#myfootnote2)</sup> (condition-concentrations and conditions only) on simulated data. 

- `code/dosage_sim.R`: Reproduce ROC plots for comparing matrix linear models (dosage-response) and Collins et al.'s S scores <sup>[2](#myfootnote2)</sup> (condition-concentrations and conditions only). 

---

<a name="myfootnote1">1</a>. Bezanson, J., Edelman, A., Karpinski, S., and Shah, V. B. (2017). Julia: A fresh approach to numerical computing. *SIAM review*, 59(1):65–98.

<a name="myfootnote2">2</a>. Collins, S. R., Schuldiner, M., Krogan, N. J., and Weissman, J. S. (2006). A strategy for extracting and analyzing large-scale quantitative epistatic interaction data. *Genome biology*, 7(7):R63. 

<a name="myfootnote3">3</a>. Joyce, A. R., Reed, J. L., White, A., Edwards, R., Osterman, A., Baba, T., Mori, H., Lesely, S. A., Palsson, B. Ø., and Agarwalla, S. (2006). Experimental and computational assessment of conditionally essential genes in Escherichia coli. *Journal of bacteriology*, 188(23):8259–8271. 
    
<a name="myfootnote4">4</a>. Kritikos, G., Banzhaf, M., Herrera-Dominguez, L., Koumoutsi, A., Wartel, M., Zietek, M., and Typas, A. (2017). A tool named iris for versatile high-throughput phenotyping in microorganisms. *Nature microbiology*, 2(5):17014.

<a name="myfootnote5">5</a>. Nichols, R. J., Sen, S., Choo, Y. J., Beltrao, P., Zietek, M., Chaba, R., Lee, S., Kazmierczak, K. M., Lee, K. J., Wong, A., et al. (2011). Phenotypic landscape of a bacterial cell. *Cell*, 144(1):143–156. 

<a name="myfootnote6">6</a>. R Core Team (2018). *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienna, Austria.
