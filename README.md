# [<img src="306030.svg" width=4.5% title="Gear" alt="Gear"/>](306030.svg) SCRAPT (Systematic CLONET ReAssignment of Ploidy and Tumour content)![Gear]
SCRAPT was designed as a complimentary tool to visually and programatically reassign ploidy and tumour content in tumour sample previously analyzed using the CLONET workflow (https://cran.r-project.org/web/packages/CLONETv2/index.html)

SCRAPT is built on R Shiny and Plotly. SCRAPT superimposes the total copy-number (tCN/Log2R) and allelic imbalance (Beta) values derived from the CLONET packages and provides a pre-computed grid of allele-specific copy-number (asCN) states for variables fractions of tumour content and ploidy. In addition, SCRAPT programatically aligns all potentially deleterious regions derived from the observed data (Log2R<0 and Beta<1) and imputes plausible tumour content and ploidy solutions for LoH compatible solutions for these regions. 


SCRAPT aims to visually and programatically guide the correct assignement 



grid of all potential grid combinations for copy-number states that could explain the observed values of total copy-number and allelic imbalance for any given genomic regions 
