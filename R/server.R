options(shiny.maxRequestSize=500*1024^2) 

# Define server logic required to draw a histogram ----

server <- function(input, output,session,threads=5){
    sols=read_solutions(tc=seq(0,1,0.01),threads=threads)

    cutoff_table<-reactive({
        cutoff_input <- system.file("extdata","cutoff_table.txt", package="SCRAPT")
        cutoffs=read.table(cutoff_input ,header=TRUE)
    })


    input_tc<-reactive({
        input$tc 
    })%>% debounce(2000)

    input_ploidy<-reactive({
        input$ploidy
    })%>% debounce(2000)

    input_tc_range<-reactive({
        input$tc_range
    })%>% debounce(2000)

     input_tc_steps<-reactive({
        input$tc_steps
    })%>% debounce(2000)


    input_ploidy_range<-reactive({
        input$ploidy_range
    })%>% debounce(2000)

    
    input_ploidy_steps<-reactive({
        input$ploidy_steps
    })%>% debounce(2000)

    
    input_clonal_thr<-reactive({
        input$clonal_thr
    })%>% debounce(2000)


    input_gene_list<-reactive({
       input$gene_list
    })%>% debounce(2000)

    states<-reactive({
        read_states()
    })

      in_data <- reactive({
   
        if(!is.valid(input$input_files$datapath)){
          return(NULL)
        }
        
        ai=data.table::fread(input$input_files$datapath)

        names(ai)<-tolower(names(ai))

         ai=suppressWarnings(ai%>% distinct() %>%
            mutate(id=sample)%>%
            mutate(beta=ifelse(is.na(beta),1,beta))
        )

        if(!is.null(ai$patient_id)){
          ai=ai %>% mutate(id=paste0(patient_id,"_",sample))
        }

        ai=ai%>% arrange(id)
    
        return(ai)
    })

    ids<- reactive({
        req(in_data)
        unique(in_data()$id)
    })
    
  

    gene_list<-reactive({
        if(is.null(all_sample())){
          return(NULL)
        }
        genes=all_sample() %>% 
        distinct(gene_type,gene)
        categories=data.frame(
          gene_type=c("All","All","All","All"),
          gene=c("All","All_target","All_control","All_other")
        )
        genes=dplyr::bind_rows(categories,genes) %>% 
          arrange(desc(gene_type),gene)
        return(genes)
    })

    observeEvent(c(gene_list()), {
      req(gene_list())
      updateSelectizeInput(session, "gene_list", 
        choices=unstack(gene_list(), gene ~ gene_type))
    })



  
    observeEvent(ids(), {
      req(ids)
      updateSelectInput(session, "sample", choices=ids())
    })

    all_sample<-reactive({

          if(!is.valid(input$sample,input$no_snps,in_data())){
            return(NULL)
          }
         
         
          rtd=in_data() %>% filter(id %in% unlist(input$sample))
          if(
            input$cn_mode=="focal_log2"
          ){
            rtd=rtd %>% mutate(all_log2=focal_log2)  
          }else if(
             input$cn_mode=="all_log2_right"
          ){
            rtd=rtd %>% mutate(all_log2=all_log2_right) 
          }else if( input$cn_mode=="all_log2_left"){
            rtd=rtd %>% mutate(all_log2=all_log2_left) 
          }
          rtd=rtd %>% mutate(obs_total=2*(2**all_log2))%>%
          mutate(
                tc=ifelse(is.na(tc),0,tc),
                obs_log2=as.numeric(all_log2),
                obs_total=ifelse(chr=="X",obs_total/2,obs_total),
                obs_beta=ifelse(chr=="X",-0.1,as.numeric(beta))
          ) %>%
          distinct()  %>%
          mutate(
                  state=ifelse(chr=="X",
                          paste0("[",cnb.int,"/",cna.int,"]"),
                          paste0("(",cnb.int,"/",cna.int,")")         
                  )
          )

          if(!input$no_snps){
            rtd=rtd %>% filter(!is.na(beta_n))
          }
          return(rtd)
          
    })

    all_sample_gene<-reactive({
       if(!is.valid(all_sample())){ 
          return(NULL)
       }
      this_dat=all_sample()
      this_dat$highlight=0
      this_dat$highlight2=1
      this_dat$highlight=unlist(lapply(this_dat$gene,FUN=function(x){
        if(any(input_gene_list() %in% x)){
          return(1)
        }
        return(0)
      }))


      this_dat$highlight2=unlist(lapply(this_dat$gene,FUN=function(x){
        if(any(input_gene_list() %in% x)){
          return(1)
        }
        return(0.25)
      }))

      return(this_dat)
    })

  


    current_sample<-reactive({
      unlist(input$sample)[1]
    })

    this_sample<-reactive({
      if(!is.valid(all_sample_gene(),current_sample())){
        return(NULL)
      }
      all_sample_gene() %>% 
      dplyr::filter(id==current_sample())
    })

    tc<-reactive({
      req(this_sample)
      as.numeric(unique(this_sample()$tc))
    })

    ploidy<-reactive({
      req(this_sample)
      as.numeric(unique(this_sample()$ploidy))
    })


  observeEvent(tc(),{
        req(tc)
        updateNumericInput(
          inputId="tc",
          value=tc()
        )      
  })

  observeEvent(ploidy(),{
        req(ploidy)
        updateNumericInput(
          inputId="ploidy",
          value=ploidy()
        )      
    })  

  observeEvent(input$reset,{
      req(tc,ploidy)
      updateNumericInput(
          inputId="tc",
          value=tc()
        ) 
      updateNumericInput(
          inputId="ploidy",
          value=ploidy()
        )           
     
  })

  observeEvent(
    input$delSol_rows_selected,
    {
      req(input$delSol_rows_selected)
      updateNumericInput(
          inputId="tc",
          value=  deletions_solutions()[input$delSol_rows_selected,]$tc
        ) 
      updateNumericInput(
          inputId="ploidy",
          value=deletions_solutions()[input$delSol_rows_selected,]$pl
        )           
  })


    space<-reactive({

          if(!is.valid(input_tc(),input_ploidy())){
            return(NULL)
          }
          sols[[input_tc()*100+1]]%>% 
          mutate(X=ifelse(chrom=="X",2**(X2),2*2**(X2))) %>%
          mutate(
            X2=X2-(log2(input_ploidy()/2)),
            obs_log2=obs_log2-(log2(input_ploidy()/2)))%>%
          mutate(
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



  deletions_solutions<-reactive({
    req(current_sample,all_sample_gene,input_clonal_thr)
      generate_del_solutions(
        smpl=all_sample_gene(),
        sols=sols,
        clonal_thr = input_clonal_thr()
      )
  })
    

  solutions<-reactive({

        req(input_tc_range,input_ploidy_range,input_tc_steps,input_ploidy_steps,input_tc,input_ploidy,this_sample)

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
              threads=threads
        )

    })
  

output$newPlot <- renderPlotly({  
        print(input$gene_list)
      if(!is.valid(input$input_files$datapath)){
            p1=ggplot()+geom_text(aes(x=1,y=1,label="Please import data",size=12))+theme_void()
      }else{
            if(!is.valid(current_sample())){
              p1=ggplot()+geom_text(aes(x=1,y=1,label="Please select sample/s",size=12))+theme_void()
            }else{

              if(!is.valid(all_sample_gene(),input_ploidy(),input_tc(),space(),states(),input_clonal_thr())){
                p1=ggplot()+geom_text(aes(x=1,y=1,label="Loading data please wait...",size=12))+theme_void()
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



# output$allPlot <- renderPlotly({
#       req(
#           input_tc,input_ploidy,
#           current_sample,
#           input$plot_type,
#           all_sample_gene,space,
#           states
#       )

#       plot_space(
#         tc=input_tc(),
#         ploidy=input_ploidy(),
#         samples=all_sample(),
#         loc=space(),
#         plot_type=input$plot_type,
#         states=states()
#       )
# })

output$oldPlot <- renderPlotly({



      if(!is.valid(input$input_files$datapath)){
            p1=ggplot()+geom_text(aes(x=1,y=1,label="Please import data",size=12))+theme_void()
      }else{
            if(!is.valid(current_sample())){
              p1=ggplot()+geom_text(aes(x=1,y=1,label="Please select sample/s",size=12))+theme_void()
            }else{

              if(!is.valid(all_sample_gene(),input_ploidy(),input_tc(),space(),states(),input_clonal_thr())){
                p1=ggplot()+geom_text(aes(x=1,y=1,label="Loading data please wait...",size=12))+theme_void()
              }else{
                  old_sols=read_solution(tc=tc())%>% 
                        mutate(
                          X2=X2-(log2(ploidy()/2)),
                          obs_log2=obs_log2-(log2(ploidy()/2))
                    )%>%
                    mutate(
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

  output$solutions<-renderPlotly({
    req(solutions)
    plot_distances(solutions=solutions())
  })

  output$download <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(paste0(unique(this_sample()$sample),".T",input_tc()*100,"_P",input_ploidy()*100)), ".tsv")
    },
    content = function(file) {
      vroom::vroom_write(
       this_sample() %>%
         process_sample(tc=input_tc(),ploidy=input_ploidy(),cutoff_table=cutoff_table()) %>%
         update_table(replace=input$replace,notes=input$notes), file)
    }
  )



output$delSol<-DT::renderDT(
    deletions_solutions(),
    selection = 'single'
)

   
  
}



   
