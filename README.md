# [<img src="306030.svg" width=4.5% title="Gear" alt="Gear"/>](306030.svg) SCRAPT (Systematic CLONET ReAssignment of Ploidy and Tumour content)
SCRAPT was designed as a complimentary tool to visually and programatically reassign ploidy and tumour content in tumour sample previously analyzed using the [CLONET](https://cran.r-project.org/web/packages/CLONETv2/index.html) workflow 

## INSTALLATION
SCRAPT package can be installed using the following command: 

`devtools::install_github("TearsWillFall/SCRAPT")`


## USE
SCRAPT can be called in R as follows:  

`library(scrapt)`  
`scrapt()`

This command will invoke R Shiny server locally

## DESCRIPTION

SCRAPT is an R package built using R Shiny and Plotly to facilitate the visualization and interpretation of allele-specific copy-number (asCN) states. It integrates results from the CLONET workflow, superimposing total copy-number (tCN/Log2R) and allelic imbalance (Beta) values onto a pre-computed grid of potential asCN states across a user-defined range of tumor content and ploidy (**Figure 1**). 
<p align="center">
  <img alt="Grid" src="https://github.com/user-attachments/assets/73c5b902-377c-44ef-a8d0-e1700a995499" />
</p>

<p align="center"> 
    <i><b>Figure 1</b>: Example of SCRAPT solution grid for a solution of tumour content 0.89 and ploidy 3.86</i> 
</p>


SCRAPT systematically guides users in evaluating the fitness of proposed solutions and provides an interactive interface to semi-manually refine assignments for improved accuracy. To achieve this, SCRAPT (**Figure 2**) :

- Aligns all likely deleterious regions based on observed data (Log2R < 0 and Beta < 1)
- Estimates the most compatible tumor content and ploidy values for each loss of heterozygosity (LoH)-compatible region
- Ranks computed solutions using multiple fitness metrics, including the Subclonality Index, Subclonality Distance, and Number of Clonal Regions
  
Users can manually explore different solutions by overlaying the data with each candidate grid, allowing for a visual and systematic assessment of solution fitness.

<p align="center">
  <img alt="LoH" src="https://github.com/user-attachments/assets/d09e6db9-4017-4ddb-9cd5-2ef1cef83e1c" />
</p>

<p align="center"> 
    <i><b>Figure 2</b>: Example of programatically computed solutions for potentially deleterious regions (Log2R < 0 and Beta < 1) and their ranking based on predicted
      subclonality</i> 
</p>


Additionally, SCRAPT includes a solver to compute nearby solutions within a selected tumor content and ploidy ranges fo help further refine the pre-computed solutions (**Figure 3**).

<p align="center">
  <img alt="Solver" src="https://github.com/user-attachments/assets/f2bca6fa-2948-418a-8fda-bad3b086678a" />
</p>

<p align="center"> 
    <i><b>Figure 3</b>: Example of near-solution space for a solution of tumour content 0.89 and ploidy 3.86 with a range of +/-0.5 and +/-3, for tumour content and ploidy       respectively</i> 
</p>

## EXPLANATION
In pure tumours, a.k.a cell-lines, the observed copy-number values for the major ($f^A$) and minor alleles ($f^B$) for any given genomc region can be directly extrapolated from the observed values of allelic imbalance,  $β_i^{Obs}$, and total copy-number, $CN_i^{Obs}$, such that $β_i^{Real}=β_i^{Obs}$, and $CN_{(i,μ)}^{Real}=CN_{(i,μ)}^{Obs}$.

<p align="center"> 
  <img width="600" alt="image" src="https://github.com/user-attachments/assets/119b5067-2d3c-4d39-a0d1-fb51b1bdd46e" />
</p>

In truth, tumour biopsies are often impure, containing admixtures of tumour and non-tumour cells that can harbour multiple distinct clones. The resulting heterogeneity of these admixtures makes challenging resolving the real genomic states of bulk tumours.

In tumour and normal cell admixtures, the increased levels of normal cells, approximate the observed values for $β_i$ , and $CN_{(i,μ)}$ in the sample to those present in normal diploid cells, and the values for $f_i^A$ and $f_i^B$ to the value of 1, that are expected in the normal diploid genomes.

Under this premise, the observed total copy-number values in the admixture of cells can be modelled as follows:

<p align="center"> 
  <img width="600" alt="image" src="https://github.com/user-attachments/assets/69134a93-a284-4a82-9893-312838d8d218" />
</p>

Similarly, the observed fraction of neutral reads within the admixture of cells can be modelled as follows:

<p align="center"> 
  <img width="600" alt="image" src="https://github.com/user-attachments/assets/a37ea1c1-15c6-4342-af06-ca277287c8c9" />
</p>

Beltran et al., elegantly formulated the definition of the admixture problem for each allele independetely, and modelled it under the variables of admixture of normal cells, $\lambda$, and ploidy, $\textmu$, in the following equation: 

<p align="center">
  <img width="600" alt="image" src="https://github.com/user-attachments/assets/7eee192e-1ff8-4c3d-a42b-98d09a2c2e75" />
</p>

It was proposed that under this system of equations, values for the admixture of normal cells, $\lambda$, and ploidy, $\textmu$, can be defined such that best explain the values observed for the major, $f_i^A$ and minor alleles, $f_i^B$, in a particular bulk tumour sample.

In SCRAPT, we exploited the intrisec property of these equations to:
  1. Derive the expected values of observed proportion of neutral reads, $β_i^{Obs}$, and copy-number, $CN_{(i,μ)}^{Obs}$ for a series of known copy-number states (0/0, 0/1,  1/1, ...) and variable fraction of normal cells (0.01 ... 0.99), $\lambda$, and ploidy (2 ... 8), $\textmu$.
  2. Pre-compute a grid to define the minimum and maximum range of observed proportion of neutral reads $β_i^{Obs}$, and copy-number, $CN_{(i,μ)}^{Obs}$ that each copy-number states can be defined within for each combination of fraction of normal cells, $\lambda$ and ploidy, $\textmu$.
  3. Identify solutions of fraction of normal cells, $\lambda$ and ploidy, $\textmu$, that would align deleterious regions (Log2R < 0 and Beta < 1) to all evaluable LoH state (0/0,0/1,0/2,0/3...)




