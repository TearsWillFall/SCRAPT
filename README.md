# [<img src="306030.svg" width=4.5% title="Gear" alt="Gear"/>](306030.svg) SCRAPT (Systematic CLONET ReAssignment of Ploidy and Tumour content)
SCRAPT was designed as a complimentary tool to visually and programatically reassign ploidy and tumour content in tumour sample previously analyzed using the CLONET workflow (https://cran.r-project.org/web/packages/CLONETv2/index.html)

## DESCRIPTION

SCRAPT is built on R Shiny and Plotly. SCRAPT superimposes the total copy-number (tCN/Log2R) and allelic imbalance (Beta) values derived from the CLONET packages and provides a pre-computed grid of allele-specific copy-number (asCN) states for variables fractions of tumour content and ploidy. 

SCRAPT  programatically aligns all potentially deleterious regions derived from the observed data (Log2R<0 and Beta<1) and imputes plausible tumour content and ploidy for LoH compatible solutions across these regions.  Precomputed solutions are ranked using several metrics (Subclonality Index, Subclonality Distance, Number of clonal regions...) allowing the user to systematically evaluate the fitness of each solution. SCRAPT allows the user to rescale their data based on the selected solution. 

SCRAPT incorporates a solver to calculate solutions within a selected range of tumour content and ploidy. 

----
## EXPLANATION
SCRAPT uses the principles of the CLONET paradigm:

<img width="600" alt="image" src="https://github.com/user-attachments/assets/9cc7dbab-8b30-4b3e-a15e-f3536747c4c7" />

----
## INSTALLATION
Package can be installed in R using the following command:\
`devtools::install_github("TearsWillFall/SCRAPT")`

----
## USE
R Shiny application can be invoked using the following command\
`library(scrapt)\
scrapt()`
