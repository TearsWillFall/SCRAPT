# [<img src="306030.svg" width=4.5% title="Gear" alt="Gear"/>](306030.svg) SCRAPT (Systematic CLONET ReAssignment of Ploidy and Tumour content)
SCRAPT was designed as a complimentary tool to visually and programatically reassign ploidy and tumour content in tumour sample previously analyzed using the CLONET workflow (https://cran.r-project.org/web/packages/CLONETv2/index.html)

---

## DESCRIPTION

SCRAPT is an R package built using R Shiny and Plotly to facilitate the visualization and interpretation of allele-specific copy-number (asCN) states. It integrates results from the CLONET workflow, superimposing total copy-number (tCN/Log2R) and allelic imbalance (Beta) values onto a pre-computed grid of potential asCN states across a user-defined range of tumor content and ploidy (**Figure 1**). 

![Grid](https://github.com/user-attachments/assets/73c5b902-377c-44ef-a8d0-e1700a995499)  
***Figure 1**: Example of SCRAPT solution grid for a solution of tumour content 0.89 and ploidy 3.86*

SCRAPT systematically guides users in evaluating the fitness of proposed solutions and provides an interactive interface to semi-manually refine assignments for improved accuracy. To achieve this, SCRAPT (**Figure 2**) :

  -Aligns all likely deleterious regions based on observed data (Log2R < 0 and Beta < 1)
  -Estimates the most compatible tumor content and ploidy values for each loss of heterozygosity (LoH)-compatible region
  -Ranks computed solutions using multiple fitness metrics, including the Subclonality Index, Subclonality Distance, and Number of Clonal Regions
  
Users can manually explore different solutions by overlaying the data with each candidate grid, allowing for a visual and systematic assessment of solution fitness.

![Solutions](https://github.com/user-attachments/assets/d09e6db9-4017-4ddb-9cd5-2ef1cef83e1c)  
***Figure 2**: Example of programatically computed solutions for potentially deleterious regions (Log2R < 0 and Beta < 1) and their ranking based on predicted subclonality*

![Solver](https://github.com/user-attachments/assets/f2bca6fa-2948-418a-8fda-bad3b086678a)

Additionally, SCRAPT includes a solver to compute nearby solutions within a selected tumor content and ploidy ranges fo help further refine the pre-computed solutions (**Figure 3**). 
---
## EXPLANATION
SCRAPT uses the principles of the CLONET paradigm:

<img width="600" alt="image" src="https://github.com/user-attachments/assets/9cc7dbab-8b30-4b3e-a15e-f3536747c4c7" />

----
## INSTALLATION
SCRAPT package can be installed using the following command: 

`devtools::install_github("TearsWillFall/SCRAPT")`

----
## USE
SCRAPT can be called in R as follows:  

`library(scrapt)`  
`scrapt()`

This command will invoke R Shiny server locally

