
#' Start SCRAPT shiny server
#'
#' @param threads  Number of threads to launch SCRAPT
#' @import shiny
#' @import tidyverse
#' @import dplyr
#' @import ggplot2
#' @import plotly
#' @import patchwork
#' @import DT
#' @return
#' @export
#'
#' @examples





scrapt=function(idata=NULL,threads=NULL){
        options(shiny.maxRequestSize=30*1024^2)
        if(is.null(threads)){
            orange <- crayon::make_style("orange")
            cat(orange("[Warning] Number of cores was not specified\n"))
            cat(orange("[Warning] Detecting maximum number of cores available\n"))
            threads=parallel::detectCores()-1
        }

        cat(orange(paste0("[Warning] SCRAPT will use ",threads," cores\n" )))
      
        server <- function(input, output,session){

            #### Read grid solutions
            sols=read_solutions(
                tc=seq(0,1,0.01),
                threads=threads,
                loc=normalizePath(system.file("data/sols",package="SCRAPT"))
            )

            #### Read cutt off talbe
            cutoff_table<-shiny::reactive({
                cutoff_input <- system.file("data","cutoff_table.txt", package="SCRAPT")
                cutoffs=read.table(cutoff_input ,header=TRUE)
            })

            #### Delay button response
            input_tc<-shiny::reactive({
                input$tc 
            })%>% shiny::debounce(2000)

            input_ploidy<-shiny::reactive({
                input$ploidy
            })%>% shiny::debounce(2000)

            input_tc_range<-shiny::reactive({
                input$tc_range
            })%>% shiny::debounce(2000)

            input_tc_steps<-shiny::reactive({
                input$tc_steps
            })%>% shiny::debounce(2000)


            input_ploidy_range<-shiny::reactive({
                input$ploidy_range
            })%>% shiny::debounce(2000)

            
            input_ploidy_steps<-shiny::reactive({
                input$ploidy_steps
            })%>% shiny::debounce(2000)

            
            input_clonal_thr<-reactive({
                input$clonal_thr
            })%>% shiny::debounce(2000)


            #### List of genes to highlight
            input_gene_list<-shiny::reactive({
                genes_to_add=genes()%>% dplyr::filter(gene_type %in% input$gene_list)
                unique(c(genes_to_add$gene,input$gene_list))
            })%>% shiny::debounce(2000)


            #### List of genes to remove from analysis
            input_gene_list2<-shiny::reactive({
                genes_to_add=genes()%>% dplyr::filter(gene_type %in% input$gene_list2)
                unique(c(genes_to_add$gene,input$gene_list2))
            })%>% shiny::debounce(2000)


            #### Read input states 
            states<-shiny::reactive({
                read_states()
            })



            #### Read input data
            in_data <- shiny::reactive({

                if(!is.valid(input$input_files$datapath)){
                    if(is.valid(idata)){
                        ai=data.table::fread(idata)
                    }else{
                        return(NULL)
                    }
        
                }else{
                    ai=data.table::fread(input$input_files$datapath)
                }
                
                

                names(ai)<-tolower(names(ai))

                ai=suppressWarnings(ai%>% 
                    dplyr::distinct() %>%
                    dplyr::mutate(
                        id=sample
                    )%>%
                    dplyr::mutate(beta=ifelse(is.na(beta),1.1,beta))
                )

                if(!is.null(ai$patient_id)){
                ai=ai %>% dplyr::mutate(id=paste0(patient_id,"_",sample))
                }

                ai=ai%>% dplyr::arrange(id)
            
                return(ai)
            })

            

             #### Sample id list in input data
            ids<- shiny::reactive({
                shiny::req(in_data)
                unique(in_data()$id)
            })


               ###Assign current sample in all samples
            current_sample<-shiny::reactive({
                if(!is.valid(input$sample)){
                    return(NULL)
                }
                unlist(input$sample)[1]
            })



            #### Gene selection
            genes<-shiny::reactive({
                if(!is.valid(in_data())){
                    return(NULL)
                }

                if(input$gene_cat=="Type"){
                    genes=in_data() %>% 
                    dplyr::distinct(gene_type,gene)%>%
                    dplyr::arrange(order(gtools::mixedorder(gene_type)),order(gtools::mixedorder(gene)))
                    
                }else if (input$gene_cat=="Chromosome"){
                    genes=in_data() %>% 
                    dplyr::distinct(chr,gene)%>% 
                    dplyr::mutate(chr=paste0("chr",chr))%>%
                    dplyr::rename(gene_type=chr)%>% 
                    dplyr::arrange(order(gtools::mixedorder(gene_type)),order(gtools::mixedorder(gene)))
                }

                categories=data.frame(
                    gene_type="All",
                    gene=unique(genes$gene_type)
                )%>% 
                    dplyr::arrange(
                        order(gtools::mixedorder(gene))
                    )
                
                genes=dplyr::bind_rows(categories,genes)
                
                print(genes)
                return(genes)
               
            })
            
             #### Update selected genes in gene list and sample list
            shiny::observeEvent(c(genes()), {
                if(!is.valid(genes())){
                    return(NULL)
                }
                ### Conver to list of lists
                choices=unstack(genes(), gene ~ gene_type)
                ### Sort in natural order
                choices=choices[gtools::mixedsort(names(choices))]

                shiny::updateSelectizeInput(session, "gene_list", 
                choices=choices)

                shiny::updateSelectizeInput(session, "gene_list2", 
                choices=choices)
            })


            

            #### Update selected genes in gene list and sample list
            shiny::observeEvent(c(ids()), {
                if(!is.valid(ids())){
                    return(NULL)
                }
                shiny::updateSelectInput(session, "sample", choices=ids())
            })


            #### Reset sample selection when new data is loaded
            shiny::observeEvent(input$input_files$datapath,{
                shiny::updateSelectizeInput(session, "gene_list", 
                choices=NULL, selected=NULL)
                shiny::updateSelectizeInput(session, "gene_list2", 
                choices=NULL, selected=NULL)
                shiny::updateSelectInput(session,"sample",selected=NULL,choices=NULL)}
            )

            
        
            #### Process input data
            all_sample<-shiny::reactive({

                if(!is.valid(in_data(),ids(),genes(),input$sample)){
                    return(NULL)
                }
                
                rtd=in_data() %>% dplyr::filter(id %in% unlist(input$sample))
                
                if(
                    input$cn_mode=="focal_log2"
                ){
                    rtd=rtd %>% dplyr::mutate(all_log2=focal_log2)  
                }else if(
                    input$cn_mode=="all_log2_right"
                ){
                    rtd=rtd %>% dplyr::mutate(all_log2=all_log2_right) 
                }else if( input$cn_mode=="all_log2_left"){
                    rtd=rtd %>% dplyr::mutate(all_log2=all_log2_left) 
                }
                rtd=rtd %>% dplyr::mutate(obs_total=2*(2**all_log2))%>%
                mutate(
                        tc=ifelse(is.na(tc),0,tc),
                        obs_log2=as.numeric(all_log2),
                        obs_total=ifelse(chr=="X",obs_total/2,obs_total),
                        obs_beta=ifelse(chr=="X",-0.1,as.numeric(beta))
                ) %>%
                dplyr::distinct()  %>%
                dplyr::mutate(
                        state=ifelse(chr=="X",
                                paste0("[",cnb.int,"/",cna.int,"]"),
                                paste0("(",cnb.int,"/",cna.int,")")         
                        )
                )

                if(!input$no_snps){
                    rtd=rtd %>% dplyr::filter(!is.na(beta_n))
                }
                return(rtd)
                
            })


            #### Highlight/Mantain selected genes
            all_sample_gene<-shiny::reactive({
                
                if(!is.valid(all_sample())){
                    return(NULL)
                }

                this_dat=all_sample()
                this_dat$highlight=0
                this_dat$highlight2=1
                this_dat$show=1
            
                this_dat$highlight=unlist(lapply(this_dat$gene,FUN=function(x){
                    if(any(input_gene_list() %in% x)){
                    return(1)
                    }
                    return(0)
                }))

                this_dat$show=unlist(lapply(this_dat$gene,FUN=function(x){
                    if(any(input_gene_list2() %in% x)){
                    return(0)
                    }
                    return(1)
                }))


                this_dat$highlight2=unlist(lapply(this_dat$gene,FUN=function(x){
                    if(any(input_gene_list() %in% x)){
                    return(1)
                    }
                    return(0.25)
                }))


                return(this_dat %>% dplyr::filter(show==1))
            })

        

            ###Select data in current sample
            this_sample<-shiny::reactive({
                if(!is.valid(all_sample_gene(),current_sample())){
                    return(NULL)
                }
                all_sample_gene() %>% 
                dplyr::filter(id==current_sample()) 
            })

            ###Assign original tc in sample
            tc<-shiny::reactive({
                shiny::req(this_sample)
                as.numeric(unique(this_sample()$tc))
            })

            ###Assign original ploidy in sample
            ploidy<-shiny::reactive({
                shiny::req(this_sample)
                as.numeric(unique(this_sample()$ploidy))
            })

            ###Observe changes in tc and update
            shiny::observeEvent(tc(),{
                    shiny::req(tc)
                    shiny::updateNumericInput(
                    inputId="tc",
                    value=tc()
                    )      
            })
            
             ###Observe changes in ploidy and update
            shiny::observeEvent(ploidy(),{
                    shiny::req(ploidy)
                    shiny::updateNumericInput(
                    inputId="ploidy",
                    value=ploidy()
                    )      
                })  
            
            ###Observe reset button and reset tc and ploidy if truee
            shiny::observeEvent(input$reset,{
                shiny::req(tc,ploidy)
                shiny::updateNumericInput(
                    inputId="tc",
                    value=tc()
                    ) 
                shiny::updateNumericInput(
                    inputId="ploidy",
                    value=ploidy()
                    )           
                
            })



            ### Observe selection in LoH table and assign as current tc and ploidy solution
            shiny::observeEvent(
                input$delSol_rows_selected,
                {
                shiny::req(input$delSol_rows_selected)
                shiny::updateNumericInput(
                    inputId="tc",
                    value=  deletions_solutions()[input$delSol_rows_selected,]$tc
                    ) 
                shiny::updateNumericInput(
                    inputId="ploidy",
                    value=deletions_solutions()[input$delSol_rows_selected,]$pl
                    )           
            })

            ### Get grid for current tc and ploidy
            space<-shiny::reactive({
                if(!is.valid(input_tc(),input_ploidy())){
                    return(NULL)
                }
                sols[[input_tc()*100+1]]%>% 
                dplyr::mutate(X=ifelse(chrom=="X",2**(X2),2*2**(X2))) %>%
                dplyr::mutate(
                    X2=X2-(log2(input_ploidy()/2)),
                    obs_log2=obs_log2-(log2(input_ploidy()/2)))%>%
                dplyr::mutate(
                    X=ifelse(chrom=="X",
                    (2**X2),
                    2*(2**X2)
                    ),
                    obs_total=ifelse(chrom=="X",
                    (2**obs_log2),
                    2*(2**obs_log2)
                    )
                )
            })



            ### Solve all LoH solutions for current grid
            deletions_solutions<-shiny::reactive({
             
                generate_del_solutions(
                    smpl=all_sample_gene(),
                    sols=sols,
                    clonal_thr = input_clonal_thr()
                )
            
            })
            
            ### Solve solution for near tc and ploidy space
            solutions<-reactive({

                    shiny::req(input_tc_range,input_ploidy_range,input_tc_steps,input_ploidy_steps,input_tc,input_ploidy,this_sample,input_clonal_thr)

                    min_tc=input_tc()-input_tc_range()
                    if(min_tc<0){
                    min_tc=0
                    }
                    max_tc=input_tc()+input_tc_range()
                    if(max_tc>1){
                    max_tc=1
                    }
                    min_ploidy=input_ploidy()-input_ploidy_range()
                    if(min_ploidy<1.5){
                    min_ploidy=1.5
                    }
                    max_ploidy=input_ploidy()+input_ploidy_range()
                    tc_steps=input_tc_steps()
                    ploidy_steps=input_ploidy_steps()
                    range_tc=seq(min_tc,max_tc,tc_steps)
                    range_ploidy=seq(min_ploidy,max_ploidy,ploidy_steps)

                    generate_near_solution(
                        sample=this_sample(),
                        range_tc=range_tc,
                        range_ploidy=range_ploidy,
                        threads=threads,
                        clonal_thr=input_clonal_thr()
                    )

            })

            ### Render grid with current tc and ploidy solution
            output$newPlot <-   plotly::renderPlotly({      
                
                    if(!is.valid(input$input_files$datapath)&!is.valid(idata)){
                            p1=ggplot2::ggplot()+ggplot2::geom_text(ggplot2::aes(x=1,y=1,label="Please import data",size=12))+theme_void()
                    }else{
                            if(!is.valid(current_sample())){
                                p1=ggplot2::ggplot()+ggplot2::geom_text(ggplot2::aes(x=1,y=1,label="Please select sample/s",size=12))+theme_void()
                            }else{
                                if(!is.valid(this_sample(),input_ploidy(),input_tc(),input_clonal_thr(),space(),states())){
                                    p1=ggplot2::ggplot()+ggplot2::geom_text(ggplot2::aes(x=1,y=1,label="Loading data please wait...",size=12))+theme_void()
                                }else{
                                    p1=plot_space(
                                    tc=input_tc(),
                                    ploidy=input_ploidy(),
                                    samples=all_sample_gene(),
                                    loc=space(),
                                    plot_type=input$plot_type,
                                    states=states(),
                                    clonal_thr = input_clonal_thr()
                                    )
                                }  
                            }
                    }
                    p1 
                })


            ### Render grid with original tc and ploidy solution for referencing

            output$oldPlot <-   plotly::renderPlotly({ 
                
                if(!is.valid(input$input_files$datapath)){
                    p1=ggplot2::ggplot()+ggplot2::geom_text(ggplot2::aes(x=1,y=1,label="Please import data",size=12))+theme_void()
                }else{
                        if(!is.valid(current_sample())){
                        p1=ggplot2::ggplot()+ggplot2::geom_text(ggplot2::aes(x=1,y=1,label="Please select sample/s",size=12))+theme_void()
                        }else{

                        if(!is.valid(all_sample_gene(),input_ploidy(),input_tc(),space(),states(),input_clonal_thr())){
                            p1=ggplot2::ggplot()+ggplot2::geom_text(ggplot2::aes(x=1,y=1,label="Loading data please wait...",size=12))+theme_void()
                        }else{
                            old_sols=read_solution(tc=tc(),loc=normalizePath(system.file("data/sols", package="SCRAPT")))%>% 
                                    dplyr::mutate(
                                    X2=X2-(log2(ploidy()/2)),
                                    obs_log2=obs_log2-(log2(ploidy()/2))
                                )%>%
                                dplyr::mutate(
                                X=ifelse(chrom=="X",
                                    (2**X2),
                                    2*(2**X2)
                                ),
                                obs_total=ifelse(chrom=="X",
                                    (2**obs_log2),
                                    2*(2**obs_log2)
                                )
                                )

                                p1=plot_space(
                                    tc=tc(),
                                    ploidy=ploidy(),
                                    samples=this_sample(),
                                    loc=old_sols,
                                    plot_type=input$plot_type,
                                    states=states(),
                                    clonal_thr = input_clonal_thr()
                                    )

                        }}}
                    p1
            })



            ### Render solver solutions
            output$solutions <-  plotly::renderPlotly({
                shiny::req(solutions,this_sample,input_tc,input_ploidy,input_clonal_thr) 
                plot_distances(
                    sample=this_sample(),
                    solutions=solutions(),
                    tc=input_tc(),
                    ploidy=input_ploidy(),
                    clonal_thr=input_clonal_thr()
                )
            })

            build_download_button=function(){
                    shiny::downloadHandler(
                                filename = function() {
                                paste0(tools::file_path_sans_ext(paste0(unique(this_sample()$sample),".T",input_tc()*100,"_P",input_ploidy()*100)), ".tsv")
                                },
                                content = function(file) {
                                vroom::vroom_write(
                                this_sample() %>%
                                    process_sample(new_tc=input_tc(),new_ploidy=input_ploidy(),cutoff_table=cutoff_table()) %>%
                                    update_table(replace=input$replace,notes=input$notes), file)
                                }
                            )
            }


            ### Render download button
            output$download <- build_download_button()

            ### Render datatable with all LoH solution
            output$delSol<-DT::renderDT(
                deletions_solutions(),
                selection = 'single'
            )

        
  
        
        }

        ## Build UI structure
        ui <-function(input){
         

                build_main_panel=function(){
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
                    }

                build_sidebar_panel=function(){
                        
                                    shiny::sidebarPanel(
                                    shiny::fluidRow(
                                        shiny::wellPanel(
                                             style = "background-color: white;",
                                             shiny::fluidRow(h4(strong("Data Import/Export"))),
                                                shiny::fluidRow(
                                                    shiny::fileInput("input_files", 
                                                    "Choose a File:", 
                                                    multiple=TRUE,
                                                    accept = c(".csv",".txt",".tsv",".xlsx",".xlsm")
                                                )),
                                                shiny::fluidRow(shiny::selectInput("sample", "Select Sample:", choices = c("No choice"),multiple=TRUE)),
                                                shiny::fluidRow(
                                                    shiny::column(6, shiny::downloadButton(outputId="download", class = "btn-block")),
                                                    shiny::column(6, shinyWidgets::materialSwitch(inputId="replace","Override Old Calls?",value=TRUE,status="info"))
                                            ),
                                            shiny::fluidRow(shiny::textInput(
                                                inputId = "notes",
                                                label = "Notes:"
                                            ))
                                    ),

                                    shiny::wellPanel(
                                        style = "background-color: white;",
                                        shiny::fluidRow(h4(strong("Main Plot Controls"))),
                                            shiny::fluidRow(shiny::numericInput(
                                                inputId = "tc",
                                                label = "Tumour Content:",
                                                min = 0,
                                                max = 1,
                                                value=0.99,
                                                step = 0.01
                                            )),
                                        shiny::fluidRow(shiny::numericInput(
                                            inputId = "ploidy",
                                            label = "Ploidy:",
                                            min = 0,
                                            max = 8,
                                            step=0.01,
                                            value = 2
                                            )),
                                        shiny::fluidRow(shiny::numericInput(
                                            inputId = "clonal_thr",
                                            label = "Clonal Threshold :",
                                            min = 0,
                                            max = 1,
                                            step=0.05,
                                            value=0.2
                                            )),
                                        shiny::fluidRow(shiny::actionButton(
                                            inputId = "reset",
                                            label = "Reset to default"
                                            )),
                                        shiny::fluidRow(shiny::selectizeInput("gene_list", "Select Genes to Highlight:", 
                                            choices = c("No choice"),
                                            multiple=TRUE)),
                                        shiny::fluidRow(shiny::selectizeInput("gene_list2", "Select Genes to Remove:", 
                                            choices = c("No choice"),
                                            multiple=TRUE)),
                                        shiny::fluidRow(shinyWidgets::radioGroupButtons(
                                            inputId="gene_cat",
                                            label="Gene Categories:",
                                            choices=c("Chromosome","Type"),
                                            selected=c("Chromosome"))),
                                        shiny::fluidRow(shinyWidgets::radioGroupButtons(
                                            inputId="plot_type",
                                            label="Plot Type:",
                                            choices=c("CN","Log2"),
                                            selected=c("CN"))),
                                        shiny::fluidRow(shinyWidgets::radioGroupButtons(
                                            inputId="cn_mode",
                                            label="CN Mode:",
                                            choices=c("all_log2","all_log2_right","all_log2_left","focal_log2"),
                                            selected=c("all_log2"))),
                                        shiny::fluidRow(shinyWidgets::prettyToggle(inputId = "no_snps",label_on="Show All Genes",label_off="Show Only with SNPs", value = TRUE))),
                                shiny::wellPanel(
                                        style = "background-color: white;",
                                        shiny::fluidRow(h4(strong("Solver Controls"))),
                                        shiny::fluidRow(shiny::numericInput(
                                                    inputId = "tc_range",
                                                    label = "Tumour Content Range to Solve:",
                                                    min = 0,
                                                    max = 1,
                                                    step=0.01,
                                                    value=0.1
                                            )),
                                        shiny::fluidRow(shiny::numericInput(
                                                    inputId = "tc_steps",
                                                    label = "Tumour Content Steps Size:",
                                                    min = 0,
                                                    max = 1,
                                                    step=0.01,
                                                    value=0.05
                                            )),
                                        shiny::fluidRow(shiny::numericInput(
                                                    inputId = "ploidy_range",
                                                    label = "Ploidy Range to Solve:",
                                                    min = 0,
                                                    max = 8,
                                                    step=0.01,
                                                    value=1
                                            )),
                                        shiny::fluidRow(shiny::numericInput(
                                                    inputId = "ploidy_steps",
                                                    label = "Ploidy Steps Size:",
                                                    min = 0,
                                                    max = 1,
                                                    step=0.01,
                                                    value=0.1
                                            ))
                                        )
                                            
                                )
                            )
                        

                }

                shiny::fluidPage(
                    # App title ----
                    shiny::titlePanel("Systematic CLONET ReAssignment for Ploidy and Tumour content (SCRAPT)"),
                    # Sidebar panel for inputs ----
                    shiny::sidebarLayout(
                        sidebarPanel = build_sidebar_panel(),
                        # Main panel for inputs ----
                        mainPanel= build_main_panel()
                    )
                )
            }
    ### Launch Shiny App on start
    shiny::shinyApp(ui=ui, server=server)
}






