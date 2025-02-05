
# Define UI for app that draws a histogram ----
ui <-function(input){
    shiny::fluidPage(
        # App title ----
         shiny::titlePanel("CLONET ReAssignment of tumour content and  Ploidy (CRISP)"),
        # Sidebar panel for inputs ----
         shiny::sidebarLayout(
            shiny::sidebarPanel(
             shiny::fluidRow(
                fileInput("input_files", 
                  "Choose a File:", 
                  multiple=TRUE,
                  accept = c(".csv",".txt",".tsv",".xlsx",".xlsm")
                ),
                 shiny::selectInput("sample", "Select Sample:", choices = c("No choice"),multiple=TRUE) ,
                 shiny::fluidRow(
                   shiny::column(6, shiny::downloadButton(outputId="download", class = "btn-block")),
                   shiny::column(6,shinyWidgets::materialSwitch(inputId="replace","Override Old Calls?",value=TRUE,status="info"))
                ),
                shiny::textInput(
                    inputId = "notes",
                    label = "Notes:"
              ),
              shiny::numericInput(
                inputId = "tc",
                label = "Tumour Content:",
                min = 0,
                max = 1,
                value=0.99,
                step = 0.01
              ),
              shiny::numericInput(
                inputId = "ploidy",
                label = "Ploidy:",
                min = 0,
                max = 8,
                step=0.01,
                value = 2
                ),
              shiny::numericInput(
                  inputId = "clonal_thr",
                  label = "Clonal Threshold :",
                  min = 0,
                  max = 0.5,
                  step=0.01,
                  value=0.1
                ),
                shiny::actionButton(
                  inputId = "reset",
                  label = "Reset to default"
                ),
                shiny::selectizeInput("gene_list", "Select Genes:", choices = c("No choice"),
                multiple=TRUE),
              shinyWidgets::radioGroupButtons(
                inputId="plot_type",
                label="Plot Type:",
                choices=c("CN","Log2"),
                selected=c("CN")
              ),
              shinyWidgets::radioGroupButtons(
                inputId="cn_mode",
                label="CN Mode:",
                choices=c("all_log2","all_log2_right","all_log2_left","focal_log2"),
                selected=c("all_log2")
              ),
              shinyWidgets::prettyToggle(inputId = "no_snps",label_on="Show All",label_off="Show Only with SNPs", value = TRUE),
              shiny::numericInput(
                inputId = "tc_range",
                label = "Tumour Content Range to Solve:",
                min = 0,
                max = 1,
                step=0.01,
                value=0.1
          ),
            shiny::numericInput(
                inputId = "tc_steps",
                label = "Tumour Content Steps Size:",
                min = 0,
                max = 1,
                step=0.01,
                value=0.05
          ),
            shiny::numericInput(
                inputId = "ploidy_range",
                label = "Ploidy Range to Solve:",
                min = 0,
                max = 8,
                step=0.01,
                value=1
          ),
           shiny::numericInput(
                inputId = "ploidy_steps",
                label = "Ploidy Steps Size:",
                min = 0,
                max = 1,
                step=0.01,
                value=0.1
          )
         
        )
      ),
            shiny::mainPanel(
                shiny::tabsetPanel(
                  shiny::tabPanel(title = "Update",
                    shiny::fluidRow(
                      shiny::column(12,align="center",
                        plotly::plotlyOutput(outputId = "newPlot",
                        width = "8in",
                        height = "4in")
                      )
                    ),
                    shiny::fluidRow(
                      shiny::column(12,align="center",
                        DT::DTOutput(
                          outputId = "delSol"
                        )
                      )
                    )
                  ),
              shiny::tabPanel(title = "Reference",
                    shiny::fluidRow(
                      shiny::column(12,align="center",
                        plotly::plotlyOutput(outputId = "oldPlot",
                          width = "8in",
                          height = "4in"
                        )
                      )
                )
              ),
                  shiny::tabPanel(title ="Solver",
                    shiny::fluidRow(
                      shiny::column(12,align="center",
                        plotly::plotlyOutput(outputId = "solutions",
                          width = "8in",
                          height = "4in"
                          )
                      )
                )
              )
            )
    
          )

        )
      )
}
  