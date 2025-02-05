
generate_states=function(max_cn=100){
        strd=generate_comb(max_cn=max_cn)
        strd=strd %>% distinct(state,col)
        strd=strd %>% filter(!grepl("\\.",state))
        tst=data.frame(state=c("(NA/NA)","[NA/NA]"),col=c("black","black"))
        all_states=rbind(strd,tst)
        write.csv("states/states.txt",x=all_states)
}

read_states=function(){
        state_input <- system.file("extdata", "states.txt", package="SCRAPT")
        states=read.csv(state_input)
}


generate_comb=function(min_cn=0,max_cn=20,tc=1,ploidy=0){
        sols=lapply(seq(min_cn,max_cn,0.25),FUN=function(x){
                        library(tidyverse)
                                lapply(seq(min_cn,max_cn,0.25),FUN=function(z){
                                        data.frame(
                                                tc=tc,
                                                ploidy=ploidy+2,
                                                obs_total=((z+x)*tc+(1-tc)*2-ploidy),
                                                real_total=x+z,
                                                cna=as.numeric(x),cnb=as.numeric(z),
                                                chrom="autosomal",
                                                obs_beta=((x*tc+(1-tc))/((tc*(x+z)+(1-tc)*2)))*2,
                                                real_beta=x/(z+x)) %>%
                                                mutate(
                                                   obs_log2=log2(obs_total/2) 
                                                )
                                })
                        }) %>% dplyr::bind_rows() %>% filter(cna<=cnb) %>% 
                        mutate(
                            obs_log2=ifelse(is.infinite(obs_log2),-2,obs_log2),
                            obs_beta=ifelse(is.na(obs_beta)|is.infinite(obs_beta),1,obs_beta),
                            real_beta=ifelse(is.na(real_beta)|is.infinite(real_beta),1,real_beta)) %>%
                        mutate(
                            obs_log2=ifelse(is.na(obs_log2),0,obs_log2)) %>%   
                        mutate(col=ifelse(
                                real_total==0,
                                "#0000FF",
                                ifelse(real_total==1,
                                "#ADD8E6",
                                ifelse(real_total==2,
                                "#BEBEBE",
                                ifelse(
                                real_total==3,
                                "#FA8072",
                                ifelse(
                                        real_total==4,
                                        "#FF0000",
                                ifelse(
                                        real_total==5,
                                        "#EE0000",

                                ifelse(
                                        real_total==6,
                                        "#CD0000",
                                        ifelse(
                                        real_total==7,
                                        "#8B0000",
                                                "#8B0000"

                                )    
                                )
                                
                                )

                ))))))%>%
                mutate(
                        state=ifelse(chrom=="X",paste0("[",cna,"/",cnb,"]"),paste0("(",cna,"/",cnb,")"))
                ) 
                
        
                
                sols2=lapply(0,FUN=function(x){
                                lapply(seq(min_cn,max_cn,0.25),FUN=function(z){
                                        data.frame(
                                                tc=tc,
                                                ploidy=ploidy+2,
                                                obs_total=((z+x)*tc+(1-tc)-ploidy),
                                                real_total=x+z,
                                                cna=x,cnb=z,
                                                chrom="X",
                                                obs_beta=c(-0.05,-0.15),
                                                real_beta=c(-0.05,-0.15)) %>% 
                                        mutate(
                                            obs_log2=log2(obs_total)
                                        )
                                 
                                })
                        }) %>% dplyr::bind_rows() %>% filter(cna<=cnb) %>% 
                        mutate(
                            obs_log2=ifelse(is.infinite(obs_log2),-2,obs_log2),
                            obs_beta=ifelse(is.na(obs_beta)|is.infinite(obs_beta),1,obs_beta),
                            real_beta=ifelse(is.na(real_beta)|is.infinite(real_beta),1,real_beta)) %>%
                          mutate(
                            obs_log2=ifelse(is.na(obs_log2),0,obs_log2))%>%
                        mutate(col=ifelse(
                                real_total==0,
                                "#0000FF",
                                ifelse(real_total==1,
                                "#BEBEBE",
                                ifelse(real_total==2,
                                "#FA8072",
                                ifelse(
                                real_total==3,
                                "#FF0000",
                                ifelse(
                                        real_total==4,
                                        "#EE0000",
                                ifelse(
                                        real_total==5,
                                        "#CD0000",

                                ifelse(
                                        real_total==6,
                                        "#8B0000",
                                        ifelse(
                                        real_total==7,
                                        "#8B0000",
                                                "#8B0000"

                                )    
                                )
                                
                                )

                ))))))%>%
                mutate(
                        state=ifelse(chrom=="X",paste0("[",cna,"/",cnb,"]"),paste0("(",cna,"/",cnb,")"))
                ) 
     

        sols=rbind(sols,sols2) %>% 
        filter(real_total<=max_cn) %>% 
        as_tibble()
}



generate_del_solutions<-function(smpl,sols,beta_thr=1,log2_thr=0,clonal_thr=0.1){
        tryCatch({
                sols=sols %>% dplyr::bind_rows() %>% 
                filter(grepl("\\(0/",state)) %>% 
                select(-c(X.1,X,X2,Y,L1,L2)) %>% 
                arrange(tc,state)

                this_sol=sols %>%
                filter(!grepl("\\[",state)) %>% 
                distinct() %>% 
                pivot_wider(id_cols=tc,names_from=state,values_from=obs_beta)
        

                this_sol2=sols %>%filter(!grepl("\\[",state)) %>% distinct() %>% 
                pivot_wider(id_cols=tc,names_from=state,values_from=obs_log2)

                this_sols=smpl %>% filter(obs_beta<beta_thr,obs_log2<log2_thr,chr!="X")
                


                tc_solutions=lapply(this_sols$obs_beta,
                FUN=function(x){(sapply(abs(this_sol[,-1]-as.numeric(x)),which.min)-1)/100}
                ) %>% dplyr::bind_rows()

                expected_beta=tc_solutions
                for (x in 1:nrow(expected_beta)){
                for(i in 1:ncol(expected_beta)){
                        expected_beta[x,i]<-this_sol[this_sol2$tc==tc_solutions[x,i,drop=TRUE],i+1,drop=TRUE]
                }
                }

                expected_log2=tc_solutions

                for (x in 1:nrow(expected_log2)){
                for(i in 1:ncol(expected_log2)){
                        expected_log2[x,i]<-this_sol2[this_sol2$tc==tc_solutions[x,i,drop=TRUE],i+1,drop=TRUE]
                }
                }

                pl_solutions=expected_log2

                for (x in 1:nrow(pl_solutions)){
                for(i in 1:ncol(pl_solutions)){
                        pl_solutions[x,i]<-round(2*(2**(expected_log2[x,i,drop=TRUE]-as.numeric(this_sols$all_log2[x]))),2)
                }
                }

                pl_solutions$gene=this_sols$gene
                tc_solutions$gene=this_sols$gene

                tc_solutions_longer=tc_solutions %>% 
                pivot_longer(cols=!gene,values_to="tc")
                pl_solutions_longer=pl_solutions %>% 
                pivot_longer(cols=!gene,values_to="pl")

                all_solutions=dplyr::left_join(tc_solutions_longer,pl_solutions_longer)

                potential_solutions=all_solutions %>% distinct(tc,pl)

                all_dist=lapply(1:nrow(potential_solutions),FUN=function(x){
                        tc=potential_solutions[x,]$tc
                        pl=potential_solutions[x,]$pl
                        this_dist=generate_distance(smpls=smpl,new_tc=tc,new_pl=pl) %>% 
                        group_by(id)%>%
                        summarise(
                                clonal=sum(dist_o<clonal_thr,na.rm=TRUE),
                                dist=round(sum(dist_o,na.rm=TRUE),2)) %>%
                        mutate(tc=tc,pl=pl)
                }) %>% dplyr::bind_rows()

                all_solutions=dplyr::left_join(all_solutions,all_dist)
                return(all_solutions)
        },error=function(e){
                return()
        })
}


generate_locations=function(
    min_cn=0,max_cn=20,tc=1,ploidy=0
){              
        generate_comb(max_cn=max_cn,min_cn=min_cn)
        generate_area(sols=sols,max_cn=max_cn,min_cn=min_cn)
}


read_solution=function(tc,loc="extdata"){
        solution_input <- system.file(loc,paste0(tc,".sol"), package="SCRAPT")
        dplyr::as_tibble(read.csv(solution_input,header=TRUE))
}

read_solutions=function(tc=seq(0,1,0.01),loc="data",threads=5){
    cl=parallel::makePSOCKcluster(threads)
    sols<-parallel::parLapply(cl,X=tc,fun=read_solution,loc=loc)
    parallel::stopCluster(cl)
    return(sols)
}


generate_area=function(sols,max_cn,min_cn){
    integer_sols=sols %>% filter(cna%%1==0,cnb%%1==0)
    area_solutions=lapply(1:nrow(integer_sols),FUN=function(x){
    
        this_sol=integer_sols[x,]
        margin_sol=this_sol
        names(margin_sol)=paste0(names(margin_sol),"_margin")
        margin_sol$cna=  this_sol$cna
        margin_sol$cnb=  this_sol$cnb
        margin_sol$chrom=  this_sol$chrom
        margin_sol$tc=  this_sol$tc
        margin_sol$ploidy=  this_sol$ploidy
        margin_sol$state=  this_sol$state
        margin_sol$col=  this_sol$col
        margin_sol$obs_beta=  this_sol$obs_beta
        margin_sol$obs_total=  this_sol$obs_total
        margin_sol$obs_log2=  this_sol$obs_log2

        lapply(seq(-0.5,0.5,by=0.25),FUN=function(z){
                 lapply(seq(-0.5,0.5,by=0.25),FUN=function(d){
                        cna_margin=this_sol$cna+z
                        cnb_margin=this_sol$cnb+d
                        if(cna_margin>cnb_margin){
                                tmp=cna_margin
                                cna_margin==cnb_margin
                                cnb_margin=tmp
                        }
                     
                        if(cna_margin>max_cn|cnb_margin>max_cn|cnb_margin<min_cn|cna_margin<min_cn){
                                return(margin_sol)
                        }else{
                                margin_sol=sols %>% 
                                mutate(cna=as.numeric(cna),cnb=as.numeric(cnb))%>%
                                dplyr::filter(chrom==this_sol$chrom,cna==cna_margin,cnb==cnb_margin)
                                names(margin_sol)=paste0(names(margin_sol),"_margin")
                                margin_sol$cna=  this_sol$cna
                                margin_sol$cnb=  this_sol$cnb
                                margin_sol$chrom=this_sol$chrom
                                margin_sol$tc=  this_sol$tc
                                margin_sol$ploidy=  this_sol$ploidy
                                margin_sol$state=  this_sol$state
                                margin_sol$col=  this_sol$col
                                margin_sol$obs_beta=  this_sol$obs_beta
                                margin_sol$obs_total=  this_sol$obs_total
                                margin_sol$obs_log2=  this_sol$obs_log2
                                return(margin_sol)
                        }
                }) %>% dplyr::bind_rows()
        })%>% dplyr::bind_rows()
       
    })%>% dplyr::bind_rows() %>%
    mutate(
        state=ifelse(chrom=="X",paste0("[",cna,"/",cnb,"]"),paste0("(",cna,"/",cnb,")"))
    ) 
    
   
    sf_polygon=area_solutions%>%
    distinct()%>%
    sf::st_as_sf(coords=c("obs_total_margin","obs_beta_margin"))%>%
    group_by(tc,ploidy,col,chrom,state,obs_beta,obs_total,obs_log2) %>% 
    summarise()
    if(tc!=0){
        sf_polygon=sf_polygon%>%
        sf::st_cast("POLYGON") %>%
        sf::st_convex_hull()%>%
        sf::st_coordinates()%>%
        as.data.frame()%>%
        mutate(
                tc=sf_polygon$tc[L2],
                ploidy=sf_polygon$ploidy[L2],
                col=sf_polygon$col[L2],
                chrom=sf_polygon$chrom[L2],
                state=sf_polygon$state[L2],
                obs_beta=sf_polygon$obs_beta[L2],
                obs_total=sf_polygon$obs_total[L2],
                obs_log2=sf_polygon$obs_log2[L2]
                )
    }else{
         sf_polygon=sf_polygon %>% 
         mutate(Y=obs_beta,X=obs_total) %>% 
        as.data.frame() %>% select(-c(geometry)) 
    }
    
    polygon=sf_polygon %>%
    mutate(state2=gsub("\\[|\\]|\\(|\\)","",state))%>%
    mutate(X2=ifelse(chrom=="X",log2(X),log2(X/2))) %>%
    mutate(X2=ifelse(is.infinite(X2),-2,X2))
    return(polygon)
}


 build_buttons=function(
        ){
                 plot_types=as_tibble(
                    list(
                        types=list(
                            "beta_log2","beta_cn"),
                        values=list(
                            list(
                                'visible',list(TRUE,TRUE,FALSE,FALSE)
                            ),
                            list(
                                'visible',list(FALSE,FALSE,TRUE,TRUE)
                            )
                        )
                    ))



                buttons <- list()
                for (i in 1:length(plot_types)) { 
                                buttons[[i]] <- list(
                                    method = "restyle",
                                    args = list(plot_types$values[[i]]),
                                    label =list(plot_types$types[[i]]))
                return(buttons)
                }
        }


plot_space=function(
        tc,
        ploidy,
        samples,
        loc,
        states,
        plot_type="CN",
        clonal_thr
){
        n_samples=length(unique(samples$id))
        


        samples=generate_distance(
                        smpls=samples,
                        new_tc=tc,
                        new_ploidy=ploidy
                )
 
        scatter_type="none"
        if(n_samples>1){
                scatter_type="multi"
        }else if(n_samples==1){
                scatter_type="single"
        }
        
        if(plot_type=="Log2"){
                loc=loc %>% mutate(X=X2,obs_total=obs_log2)
                samples=samples %>% mutate(obs_total=obs_log2)
        }

        
        samples=samples %>% mutate(state=state_o)     
    
     
      

        state_col=unlist(states$col)
        names(state_col)=unlist(states$state)

        build_space=function(loc,state_col,samples,scatter_type,plot_type){
                fig=plot_ly()
                        fig<-fig %>% add_polygons(
                                        data=loc,
                                        x=~X,
                                        y=~Y,
                                        color=~state,
                                        colors=state_col,
                                        alpha=0.25,
                                        line=list(
                                                width=1
                                        )
                                )%>%
                                add_trace(
                                        data=loc %>%
                                        group_by(state)%>%
                                        summarise(
                                        obs_log2=mean(obs_log2),
                                        obs_total=mean(obs_total),
                                        obs_beta=mean(obs_beta)   
                                        ),
                                        name="Integer",
                                        x=~obs_total,
                                        y=~obs_beta,
                                        text=~state,
                                        marker=list(
                                                symbol='x'
                                        ),
                                        inherit=FALSE
                                )

                                if(scatter_type=="multi"){
                                        fig=fig %>%add_trace(
                                                data=samples,     
                                                x=~obs_total,
                                                y=~obs_beta,
                                                text=~paste0(
                                                        "Gene:",gene,
                                                        "\nChromosome:",chr,
                                                        "\nSNPs:",snps,
                                                        "\ncna/cnb:",round(cna_o,2),"/",round(cnb_o,2),
                                                        "\nState:",state_o,
                                                        "\nDist:",round(dist_o,5)
                                                        ),
                                                symbol=~id,
                                                inherit=FALSE
                                        )
                                }else if(scatter_type=="single"){
                                        
      print(samples %>% distinct(gene,highlight,highlight2))
                                        fig=fig %>%add_trace(
                                                data=samples,
                                                type="scatter",
                                                mode="markers",     
                                                x=~obs_total,
                                                y=~obs_beta,
                                                color=~state,
                                                colors=state_col,
                                                text=~paste0(
                                                        "Gene:",gene,
                                                        "\nChromosome:",chr,
                                                        "\nSNPs:",snps,
                                                        "\ncna/cnb:",round(cna_o,2),"/",round(cnb_o,2),
                                                        "\nState:",state_o,
                                                        "\nDist:",round(dist_o,5)
                                                        ),
                                                marker=list(opacity=~highlight2),
                                                inherit=FALSE
                                        )
                                }
                        
                        fig=fig %>% layout(
                                showlegend=FALSE,
                                title=paste0("Analysis of Samples  ",paste0(unique(samples$id),collapse="|")),
                                yaxis=list(title="Beta",range=list(-0.2,1.1)),
                                annotations=list(x = 0.2 , y = 1, text = paste(
                                        "Tumour Content: ",tc,
                                        "| Ploidy: ",ploidy,
                                        "| Subclonal Distance:",round(sum(samples$dist_o,na.rm=TRUE),5),
                                        "| Subclonal Index:",round(sum(samples$dist_o,na.rm=TRUE)/sum(!is.na(samples$dist_o)*0.5),5),
                                        "| Clonal Genes: ",sum(samples$dist_o<clonal_thr,na.rm=TRUE)
                                ),
                                showarrow = F, xref='paper', yref='paper')
                        )

                        if(plot_type=="CN"){
                                fig %>% layout(xaxis=list(title="Total Copy Number"))
                        }else{
                                fig %>% layout(xaxis=list(title="Log2R"))   
                        }     
        }

        build_space(loc=loc,state_col=state_col,samples=samples,scatter_type = scatter_type,plot_type)

      
       
}




generate_distance=function(smpls,new_tc,new_ploidy){
        
        this_tc=new_tc
        if(new_tc==0){
                this_tc=1
        }
        smpls=smpls %>% 
        rowwise()%>%
        dplyr::mutate(
                cna_o=ifelse(
                        chr=="X",
                        (2**(obs_log2+log2(new_ploidy/2))-(1-this_tc))/this_tc,
                        ((2-obs_beta)*(obs_beta*2**(obs_log2+log2(new_ploidy/2))-(1-new_tc))+2*(1-this_tc)*(1-obs_beta))/(this_tc*obs_beta)),
                cnb_o=ifelse(
                        chr=="X",0,(obs_beta*2**(obs_log2+log2(new_ploidy/2))-(1-new_tc))/this_tc)
                )

        smpls=smpls%>%
        dplyr::mutate(
            cna.int_o=ifelse(new_tc==0,1,round(cna_o)),
            cnb.int_o=ifelse(new_tc==0,1,round(cnb_o))
        )%>%
        dplyr::mutate(
            cna.int_o=ifelse(cna.int_o<0,0,cna.int_o),
            cnb.int_o=ifelse(cnb.int_o<0,0,cnb.int_o)
        )%>%
        dplyr::mutate(
            ploidy_o=new_ploidy,
            tc_o=new_tc,
            dist_o=min(0.5,((cna.int_o-cna_o)**2+(cnb.int_o-cnb_o)**2)**0.5)
        )%>%
        mutate(
                state_o=ifelse(chr=="X",
                        paste0("[",cnb.int_o,"/",cna.int_o,"]"),
                        paste0("(",cnb.int_o,"/",cna.int_o,")")         
                )
        )
        return(smpls)

}       
               
        
process_sample<-function(smpls,tc,ploidy,cutoff_table){
        smpls=smpls%>% generate_distance(new_tc=tc,new_ploidy=ploidy)%>%
        ungroup()%>%
        mutate(all_log2_corr_o=ifelse(log2_cutoff_pass,correct_log2(all_log2,ploidy,tc),all_log2))%>%
        rowwise()%>%
        mutate(cn_call_corr_o=get_lesion_type(
                evidence_n = evidence_n,
                evidence=evidence,
                gene=gene,
                all_log2 = all_log2,
                all_log2_corr = all_log2_corr_o,
                tc=tc,
                chr=chr,
                snps=snps,
                cutoff_table=cutoff_table
        )) 
        return(smpls)
}


  generate_near_solution=function(
                                sample,
                                range_tc,
                                range_ploidy,
                                threads=5
        ){
        lapply(range_ploidy,
                FUN=function(p,range_tc){
                lapply(range_tc,FUN=function(t){
                        sample %>% generate_distance(new_tc=t,new_ploidy=p) %>% 
                        group_by(tc_o,ploidy_o) %>%
                        summarise(dist=sum(dist_o,na.rm=TRUE))
                        }) %>% dplyr::bind_rows()
                },range_tc=range_tc) %>% dplyr::bind_rows()
    }


plot_distances=function(solutions){
        best_sol=solutions[which.min(solutions$dist),]
        plot_ly(solutions) %>%plotly::add_heatmap(x=~tc_o,y=~ploidy_o,z=~dist) %>%
        add_trace(
                name="Best Solution",
                data=best_sol,
                x=~tc_o,y=~ploidy_o,
                marker=list(
                        symbol='x'
                ),
                inherit=FALSE
        )%>%layout(
                title=paste0("Distances"),
                yaxis=list(title="Ploidy"),
                xaxis=list(title="Tumour Content")
        )
}


correct_log2 <- function(log2R, ploidy, tc) {
        log2_shift <- -log2(ploidy/2)
        log2_pl_corr <- log2R - log2_shift
        return(log2(pmax((2^(log2_pl_corr) - (1-tc)), 0)/tc))
}


get_lesion_type <- function( 
        evidence_n,
        evidence,
        gene,
        all_log2,
        all_log2_corr, 
        tc,
        chr, 
        thr_evidence_n = 0.5, 
        thr_evidence = 0.2,
        thr_tc = 0.15, 
        snps,
        cutoff_table,
        cutoff_pval=0.005
) {
        this_gene=gene
        # Define dynamic thresholds
        
        this_cutoff_loss=cutoff_table %>% filter(p_val==cutoff_pval,type=="Loss",gene==this_gene)
        this_cutoff_gain=cutoff_table %>% filter(p_val==cutoff_pval,type=="Gain",gene==this_gene)
        log2_cutoff_pass=!((this_cutoff_gain$value>=all_log2)&(all_log2>=this_cutoff_loss$value))
        log2_corr_cutoff_pass=!((this_cutoff_gain$value>=all_log2_corr)&(all_log2_corr>=this_cutoff_loss$value))
        ai_cutoff_pass= evidence>=thr_evidence
        ai_n_cutoff_pass= evidence_n>=thr_evidence_n

        # Gains
        thr.uncbalgain <- thr.likelygain <- 0.3
        thr.balgain <- thr.unbgain <- 0.5
    

        # Deletions
        thr.hemi <- thr.unchomo <- thr.likelyloss <- -0.5
        thr.homo <- -1

        if(is.na(snps)) {
                evidence_n <- 0
                evidence <- 0
        }

      


        # Apply cutoff on uncorrected log2R
        if (chr== "X") {
                if(!log2_cutoff_pass){ return("wt") }
        } else {
                if(ai_n_cutoff_pass){ return("Germline") }
                if(!log2_cutoff_pass&!ai_cutoff_pass){ return("wt") }                                                                                   ### in v3
        }

      

        # Calls for chrX
        if (chr == "X") {
                if (is.na(all_log2_corr)){ return(NA) }
                else{
                        if (tc <= thr_tc | is.na(tc)) {
                                if (all_log2_corr >= 0.5) { return("Likely Gain") }
                                else if (all_log2_corr<= -1){ return("Likely Loss") }
                        } else {
                                if (all_log2_corr >= 0.5) { return("Bal.Gain") }
                                else if (all_log2_corr <= -1){ return("Deletion on chrX") }
                        }
                }
        }

        # Calls on other chromosomes
        if(is.na(all_log2_corr)){ return(NA) }
        else {
                if (ai_cutoff_pass) {
                        if (all_log2_corr<= thr.hemi){ return("HemiDel") }
                        else if(all_log2_corr > thr.unbgain){ return("Unb.Gain") }
                        else if(!log2_corr_cutoff_pass){ return("CNNL") }
                        else { return("Imb_WTGrayZoneCN") }
                } else {
                        if (tc <= thr_tc | is.na(tc)) {
                                if (all_log2_corr <= thr.likelyloss) { return("Likely Loss") }
                                else if(all_log2_corr > thr.likelygain) { return("Likely Gain") }
                                else { return("NoImb_WTGrayZoneCN") }
                        } else {
                                if(all_log2_corr<= thr.unchomo & all_log2_corr> thr.homo){ return("Uncertain HomoDel") }
                                else if(all_log2_corr<= thr.homo){ return("HomoDel") }
                                else if(all_log2_corr > thr.uncbalgain & all_log2_corr <= thr.balgain) { return("Uncertain Bal.Gain") }
                                else if(all_log2_corr > thr.balgain) { return("Bal.Gain") }
                                else { return("NoImb_WTGrayZoneCN") }
                        }
                }
        }
}

update_table<-function(table,replace=TRUE,notes=""){
    table=table %>% 
    select(
        -starts_with("obs_"),-id
    ) %>%
    mutate(notes=notes)

    if(replace){
        table=table%>%
        mutate(
                tc=tc_o,
                ploidy=ploidy_o,
                cna=cna_o,
                cnb=cnb_o,
                cna.int=cna.int_o,
                cnb.int=cnb.int_o,
                state=state_o,
                all_log2_corr=all_log2_corr_o,
                cn_call_corr=cn_call_corr_o
        ) %>% 
        select(-ends_with("_o"))
    }

    

    return(table)

  }


is.valid <- function(...){
        vars=list(...)
        is.valid=Vectorize(
                function(x) {
                        return(is.null(shiny::need(x, message = FALSE)))        
        })
        return(all(is.valid(vars)))
}






