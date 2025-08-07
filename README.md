# [<img src="306030.svg" width=4.5% title="Gear" alt="Gear"/>](306030.svg) SCRAPT (Systematic CLONET ReAssignment of Ploidy and Tumour content)
SCRAPT is a complementary R Shiny tool designed to visually and programmatically reassign ploidy and tumour content in tumour samples that have previously undergone analysis using the [CLONET](https://cran.r-project.org/web/packages/CLONETv2/index.html) workflow.

## INSTALLATION
SCRAPT package can be installed using the following command: 

`devtools::install_github("TearsWillFall/SCRAPT")`


## USE
### LAUNCH SCRAPT

SCRAPT can be launched in R as follows:  

`library(scrapt)`  
`scrapt()`

This command will invoke Shiny R server locally.

<p align="center">
  <img width="2832" height="1370" alt="SCRAPT" src="https://github.com/user-attachments/assets/65e0dc03-13fc-4f57-9bdc-2ce2cd9a929a" />

</p>
<p align="center"> 
    <i><b>Figure 1:</b> SCRAPT user interface. The UI is divided into two main sections: control panels on the left and visualization panels on the right. <ul> <b>Import/Export Panel</b>: Load and export sample data and select samples. <b>Main Plot Controls</b>: Adjust tumour content and ploidy parameters for grid visualization. <b>Solver Control</b>: Interactively adjust ploidy and tumour content to refine model solutions. <b>Update Plot</b>: Displays the updated grid plot based on current settings. <b>Reference Plot</b>: Displays the plot based on the default settings. <b>Solver</b>: Displays the solution space for the solver. </ul></i></p>
</p>

### LOAD DATA
Use the Import/Export panel on the left side of the interface (Figure 1) to upload your data:

- Datasets must include one or more samples, each with a unique sample ID.
  
- Currently, only datasets processed using the PCF_SELECT workflow from CLONETv2 are supported.
  
- Other formats must be converted to the compatible PCF_SELECT format before use.



### SELECT SAMPLE
Use the dropdown menu in the Import/Export panel to select one or more samples for analysis.

- SCRAPT will focus on the first selected sample to determine initial ploidy and tumour content.

- The corresponding grid will be generated and updated accordingly.

### VISUALIZE THE GRID

Overlaid on this scatter plot is a grid showing the expected positions of different integer copy-number states (e.g. 0/1, 1/2, 2/2), calculated based on the current tumour content and ploidy values. By default, the grid uses the tumour content and ploidy provided in the imported CLONETv2 data.

Each copy-number state is color-coded to reflect increasing total copy number, allowing for intuitive interpretation of the genomic landscape.

`⚠️ Note: The maximum displayed grid size is 12 (i.e., up to 12 total copies). Extremely focal amplifications beyond this limit are still resolved internally but will not display the grid.
This constraint is currently fixed to prevent excessive I/O and rendering times.`

The grid is divided into three distinct areas: 

- 🔷 Top area:
Represents autosomal regions where no informative SNPs are available, but log2 R values are estimated.
Copy-number states in this area are shown using open parentheses, e.g. )1/1(, indicating that allele-specific copy number cannot be resolved.

- 🟩 Middle area (largest):
Represents autosomal regions where both β values and log2 R values are available.
Copy-number states are shown using closed parentheses, e.g. (1/1), indicating that the copy-number state can be resolved allele-specifically.

- 🟥 Bottom area:
Represents X chromosome regions, where only log2 R values can be estimated.
Copy-number states are shown using square brackets, e.g. [0/1], indicating that only one allele is present (as expected due to X chromosome biology in male samples or loss of heterozygosity in female samples).

`📌 These visual cues help interpret the confidence and type of copy-number information available for each genomic region.`


At the top of the scatter plot, SCRAPT displays:

- Tumour content

- Ploidy

- Subclonality index: the scaled sum of normalized Euclidean distances between observed regions and their nearest predicted integer copy-number states.

- Subclonality distance: the raw (unscaled) sum of these distances.

- Number of clonal regions: by default, a region is considered clonal if it lies within a defined threshold distance from its expected copy-number state (see definition below).

`⚠️ Note: The definition of "clonal" can be adjusted by the user within the UI parameters.`

For each dot in the scatter plot, hovering reveals detailed information about the corresponding gene region, including:

X and Y coordinates: Exact values for total copy number (or log2 R) and allelic imbalance (β value)

- Gene name: The annotated gene within the genomic region

- Number of SNPs: Informative SNPs used in the β value calculation

- Predicted copy number A and B: (cnA, cnB) values based on the current grid configuration

- Nearest integer copy-number state: e.g. (1/1), (2/0), etc.

- Euclidean distance to nearest state: Distance between the observed point and the predicted integer state in copy-number space

`ℹ️ This information is useful for assessing how well individual gene regions conform to the current tumour content and ploidy model.`

<p align="center">
  <img alt="Grid" src="https://github.com/user-attachments/assets/73c5b902-377c-44ef-a8d0-e1700a995499" />
</p>

<p align="center"> 
    <i><b>Figure 1</b>: Example of SCRAPT solution grid for a solution of tumour content 0.89 and ploidy 3.86</i> 
</p>



### EXPLORE DIFFERENT SOLUTIONS

SCRAPT offers several ways to explore tumour content and ploidy solutions.

The most straightforward method is through manual adjustment of tumour content and ploidy using the Main Plot Controls panel on the left. Modifying these values will automatically update the grid overlay on the scatter plot to reflect the new solution.

To evaluate how well a particular combination of tumour content and ploidy fits the data, users can refer to the following metrics displayed above the plot:

- Subclonality Index/Subclonality Distance
- Number of Clonal Regions

An optimal solution is typically one where the majority of regions are clonal, meaning the observed values closely align with expected copy-number states. This is reflected by a low Subclonality Index and a high number of clonal regions.

`✅ Tip: Iteratively adjusting tumour content and ploidy while monitoring these metrics can help identify the most biologically plausible solution`.

Another approach SCRAPT provides for exploring tumour content and ploidy models is through the identification of gene regions compatible with loss of heterozygosity (LoH)—e.g. copy-number state 0/1.

SCRAPT selects candidate regions based on the following criteria:

- Allelic imbalance: β value < 1

- Evidence of deletion: Total copy number < 2

For each candidate region and its corresponding potential LoH copy-number states, SCRAPT estimates a tumour content and ploidy combination that minimizes the Euclidean distance between the observed copy-number values and the expected integer LoH state for that region. This allows users to evaluate and compare multiple biologically plausible solutions, grounded in focal LoH events, using global metrics of fit and clonality.

Each resulting solution is presented in a summary table, which includes:

- The gene region and the LoH state under consideration (e.g., 0/1, 0/2, 0/3 ...)

- The estimated tumour content and ploidy that best fit the selected state

- The Euclidean distance for the fit

- The corresponding Subclonality Index

- The number of clonal regions under this solution

`🧠 This method is especially useful when searching for biologically meaningful ploidy/tumour content values grounded in known or expected LoH events.`

<p align="center">
  <img alt="LoH" src="https://github.com/user-attachments/assets/d09e6db9-4017-4ddb-9cd5-2ef1cef83e1c" />
</p>

<p align="center"> 
    <i><b>Figure 2</b>: Example of programatically computed solutions for potentially deleterious regions (Log2R < 0 and Beta < 1) and their ranking based on predicted
      subclonality</i> 
</p>


Lastly, SCRAPT implements a solver to find solutions for a range of tumour content and ploidy models starting from a base tumour and content.


<p align="center">
  <img alt="Solver" src="https://github.com/user-attachments/assets/f2bca6fa-2948-418a-8fda-bad3b086678a" />
</p>

<p align="center"> 
    <i><b>Figure 3</b>: Example of near-solution space for a solution of tumour content 0.89 and ploidy 3.86 with a range of +/-0.5 and +/-3, for tumour content and ploidy       respectively</i> 
</p>



