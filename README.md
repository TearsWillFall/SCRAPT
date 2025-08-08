# [<img src="306030.svg" width=4.5% title="Gear" alt="Gear"/>](306030.svg) SCRAPT (Systematic CLONET ReAssignment of Ploidy and Tumour content)
SCRAPT is a complementary R Shiny tool designed to visually and programmatically reassign ploidy and tumour content in tumour samples that have previously undergone analysis using the [CLONET](https://cran.r-project.org/web/packages/CLONETv2/index.html) workflow.

**References**: 
- Prandi et al. (2014) [doi:10.1186/s13059-014-0439-6](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0439-6)
- Carreira et al. (2014) [doi:10.1126/scitranslmed.3009448](https://www.science.org/doi/10.1126/scitranslmed.3009448)
- Romanel et al. (2015) [doi:10.1126/scitranslmed.aac9511](https://pubmed.ncbi.nlm.nih.gov/26537258/)
- Prandi et al. (2019) [doi:10.1002/cpbi.81](https://pmc.ncbi.nlm.nih.gov/articles/PMC6778654/)
- Francesco et al. (2022) [doi:10.1093/narcan/zcac016](https://pubmed.ncbi.nlm.nih.gov/35664542/)

## Table of Contents

 - [📦 Installation](#installation)
 - [📘 Usage](#usage)
   * [🚀 Launch](#launch)
   * [📂 Importing Data](#importing-data)
   * [🔍 Selecting Samples](#selecting-samples)
   * [🧬 Copy-Number State Grid and Scatter Plot](#copy-number-state-grid-and-scatter-plot)
   * [🔧 Exploring Tumour Content and Ploidy Solutions](#exploring-tumour-content-and-ploidy-solutions)
   * [📤 Exporting Data](#exporting-data)


## 📦Installation
SCRAPT package can be installed using the following command: 


`devtools::install_github("TearsWillFall/SCRAPT")`


##  Usage
The following section describes the use of SCRAPT tool
### Launch

SCRAPT can be launched in R as follows:  

`library(scrapt)`  
`scrapt()`

This command will invoke Shiny R server locally (**Figure 1**):

- **Control Panel (Left Panel)**:
  - **Import/Export Panel**: Load and export sample data and select samples. 
  - **Main Plot Controls**: Adjust tumour content and ploidy parameters for grid visualization.
  - **Solver Control**: Interactively adjust ploidy and tumour content to refine model solutions.
 
- **Visualization (Right Plot)**:
  - **Update Plot**: Displays the updated grid plot based on current settings and LoH table.
  - **Reference Plot**: Displays the plot based on the default settings.
  - **Solver**: Displays the solution space for the solver.


<p align="center">
  <img width="2832" height="1370" alt="SCRAPT" src="https://github.com/user-attachments/assets/65e0dc03-13fc-4f57-9bdc-2ce2cd9a929a" />

</p>
<p align="center"> 
    <i><b>Figure 1:</b> SCRAPT user interface. The UI is divided into two main sections: control panels on the left and visualization panels on the right.</i></p>
</p>

---

### Importing Data
Use the Import/Export panel on the left side of the interface (**Figures  1** and  **Figure 2**) to import data into SCRAPT:

- Datasets must include one or more samples, each with a unique sample ID.
  
- Currently, only datasets processed using the PCF_SELECT workflow from CLONETv2 are supported.
  
- Other formats must be converted to the compatible PCF_SELECT format before use.

` A toy dataset with the minimal required input is provided within the package and is located in inst/data/toy_dataset.tsv`

<p align="center">
  <img alt="Loading Data" src="https://github.com/user-attachments/assets/28e570f5-83d0-45f6-8023-967e0b5c4731" />
</p>

<p align="center"> 
    <i><b>Figure 2</b>: Import/export data panel</i> 
</p>


---

### Selecting Samples

Use the dropdown menu in the Import/Export panel (**Figure 2**) to select one or more samples for analysis.

- SCRAPT will focus on the first selected sample to determine initial ploidy and tumour content.

- The corresponding grid will be generated and updated accordingly.

---

### Copy-Number State Grid and Scatter Plot

When one or more samples are selected, SCRAPT generates a scatter plot displaying allelic imbalance (β value) on the y-axis and total copy number (or log2 R) on the x-axis for each gene region (See Main Plot Controls on how to change display mode).

Overlaid on this scatter plot is a grid showing the expected positions of different integer copy-number states (e.g. 0/1, 1/2, 2/2), calculated based on the current tumour content and ploidy values. By default, the grid uses the tumour content and ploidy provided in the imported data.

Each copy-number state is color-coded to reflect increasing total copy number, allowing for intuitive interpretation of the genomic landscape.

`⚠️ Note: The maximum displayed grid size is 12 (i.e., up to 12 total copies). Extremely focal amplifications beyond this limit are still resolved internally but will not display the grid.
This constraint is currently fixed to prevent excessive I/O and rendering times.`

The grid is divided into three distinct areas (**See Figure 3**): 

- 🔷 **Top area**:
Represents autosomal regions where no informative SNPs are available, but log2 R values are estimated.
Copy-number states in this area are shown using open parentheses, e.g. )1/1(, indicating that allele-specific copy number cannot be resolved.

- 🟩 **Middle area (largest)**:
Represents autosomal regions where both β values and log2 R values are available.
Copy-number states are shown using closed parentheses, e.g. (1/1), indicating that the copy-number state can be resolved allele-specifically.

- 🟥 **Bottom area**:
Represents X chromosome regions, where only log2 R values can be estimated.
Copy-number states are shown using square brackets, e.g. [0/1], indicating that only one allele is present (as expected due to X chromosome biology in male samples or loss of heterozygosity in female samples).

`📌 These visual cues help interpret the confidence and type of copy-number information available for each genomic region.`


At the top of the scatter plot, SCRAPT displays (**See Figure 3**):

- **Tumour content**

- **Ploidy**

- **Subclonality index**: the scaled sum of normalized Euclidean distances between observed regions and their nearest predicted integer copy-number states.

- **Subclonality distance**: the raw (unscaled) sum of these distances.

- **Number of clonal regions**: by default, a region is considered clonal if it lies within a defined threshold distance from its expected copy-number state (see definition below).

`⚠️ Note: The definition of "clonal" can be adjusted by the user within the UI parameters.`

<p align="center">
  <img alt="Grid" src="https://github.com/user-attachments/assets/b4494435-43da-4d57-b197-066717965eb4" />
</p>

<p align="center"> 
    <i><b>Figure 3</b>: Example of SCRAPT solution grid for a solution of tumour content 0.89 and ploidy 3.86 </i> 
</p>

🖱️ **Data Point Hover Information**

For each dot in the scatter plot (**See Figure 3**), hovering reveals detailed information about the corresponding gene region, including:

- **X and Y coordinates**: Exact values for total copy number (or log2 R) and allelic imbalance (β value)

- **Gene name**: The annotated gene within the genomic region

- **Number of SNPs**: Informative SNPs used in the β value calculation

- **Predicted copy number A and B**: (cnA, cnB) values based on the current grid configuration

- **Nearest integer copy-number state**: e.g. (1/1), (2/0), etc.

- **Euclidean distance to nearest state**: Distance between the observed point and the predicted integer state in copy-number space

`ℹ️ This information is useful for assessing how well individual gene regions conform to the current tumour content and ploidy model.`



🎛 **Main Plot Controls**

The Main Plot Controls panel (left side of the interface) allows you to customize how the scatter plot is displayed and interpreted ( **See Figure 4**). Available options include:

- **Specify tumour content and ploidy** – Manually set values to update the grid and metrics

- **Set clonal region threshold** – Define the maximum distance for a region to be considered clonal

- **Reset to default** – Restore tumour content and ploidy to their imported default values

- **Highlight elements** – Emphasize selected genes, chromosomes, or gene types

- **Remove elements** – Hide selected genes, chromosomes, or gene types from the plot

- **Switched categories** – Change the Highlight/Remove dropdown menu mode

- **Switch display type** – Toggle between total copy number and log2 R display on the x-axis

- **Change log2 R display mode** – Options include focal all log2, left log2, or right log2

- **Display filter** – Show all gene regions or only those containing informative SNPs

<p align="center">
  <img alt="Grid" src="https://github.com/user-attachments/assets/0304f5d2-3538-498b-a9cd-0869b93aa462" />
</p>

<p align="center"> 
    <i><b>Figure 4</b>: Main Plot Controls Panel</i> 
</p>


---

### Exploring Tumour Content and Ploidy Solutions

SCRAPT offers several ways to explore tumour content and ploidy solutions.

The most straightforward method is through manual adjustment of tumour content and ploidy using the Main Plot Controls panel on the left (**See Figure 4**). Modifying these values will automatically update the grid overlay on the scatter plot to reflect the new solution.

To evaluate how well a particular combination of tumour content and ploidy fits the data, users can refer to the following metrics displayed above the plot:

- **Subclonality Index/Subclonality Distance**
- **Number of Clonal Regions**

An optimal solution is typically one where the majority of regions are clonal, meaning the observed values closely align with expected copy-number states. This is reflected by a low Subclonality Index and a high number of clonal regions.

`✅ Tip: Iteratively adjusting tumour content and ploidy while monitoring these metrics can help identify the most biologically plausible solution`.


---

📋 **LoH-Based Solution Table**

Another approach SCRAPT provides for exploring tumour content and ploidy models is through the identification of gene regions compatible with loss of heterozygosity (LoH)—e.g. copy-number state 0/1 (**See Figure 5**).

SCRAPT selects candidate regions based on the following criteria:

- **Allelic imbalance: β value < 1**

- **Evidence of deletion: Total copy number < 2**

For each candidate region and its corresponding potential LoH copy-number states, SCRAPT estimates a tumour content and ploidy combination that minimizes the Euclidean distance between the observed copy-number values and the expected integer LoH state for that region. This allows users to evaluate and compare multiple biologically plausible solutions, grounded in focal LoH events, using global metrics of fit and clonality.

Each resulting solution is presented in a summary table, which includes:

- The gene region and the LoH state under consideration (e.g., 0/1, 0/2, 0/3 ...)

- The estimated tumour content and ploidy that best fit the selected state

- The Euclidean distance for the fit (Sum of all region euclidean)

- The corresponding Subclonality Index (Scaled Euclidean Distance between 0 and 1)

- The number of clonal regions under this solution

`🧠 This method is especially useful when searching for biologically meaningful ploidy/tumour content values grounded in known or expected LoH events.`


<p align="center">
  <img alt="LoH" src="https://github.com/user-attachments/assets/b4480711-a665-4df2-a312-e10ef70db2b6" />
</p>

<p align="center"> 
    <i><b>Figure 5</b>: Example of programatically computed solutions for potentially deleterious regions (Log2R < 0 and Beta < 1) and their ranking based on predicted
      subclonality</i> 
</p>


---

🧩 **Solver for Neighbouring Solutions**

Lastly, SCRAPT includes a solver that explores a range of tumour content and ploidy models starting from a known solution (**See Figure 6**).

The solver systematically searches the surrounding parameter space to identify nearby solutions that minimize the Subclonality Index, helping refine initial estimates and improve model fit.

The search space is configured via the Solver Control panel on the left, where users can adjust:

- **Search range** – how far from the initial tumour content and ploidy values the solver should explore

- **Step size** – the resolution of the search grid

`⚠️ The number of candidate solutions increases multiplicatively with wider ranges and smaller step sizes, which can significantly increase computation time.`

`💡 This feature is particularly useful when an initial solution is close to optimal but may benefit from fine-tuning tumour content and ploidy.`

<p align="center">
  <img alt="Solver" src="https://github.com/user-attachments/assets/f2bca6fa-2948-418a-8fda-bad3b086678a" />
</p>

<p align="center"> 
    <i><b>Figure 6</b>: Example of near-solution space for a solution of tumour content 0.89 and ploidy 3.86 with a range of +/-0.5 and +/-3, for tumour content and ploidy       respectively</i> 
</p>

###  Exporting Data
Once a new tumour content and ploidy solution has been identified, the results can be exported using the Download button in the Import/Export panel (**See Figure 2**).

By default, SCRAPT will override the original tumour content, ploidy, and copy-number state values in the imported dataset.

If you wish to preserve the original values and instead save the updated results as new columns, disable the Override Calls option before exporting.






