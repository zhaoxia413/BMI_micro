#' return the frequency and a bar plot of the query microbes and the plotdata
#' 
#' @param Disease: the sudty coditions including: 
#' 1: acute_diarrhoea             2: ACVD                     
#' 3: adenoma                     4: AS                       
#' 5: BD                          6: bronchitis               
#' 7: CDI                         8: cephalosporins           
#' 9: cirrhosis                  10: control                  
#' 11: cough                      12: CRC                      
#' 13: cystitis                   14: fatty_liver              
#' 15: fever                      16: FMT                      
#' 17: hypertension               18: IBD                      
#' 19: IGT                        20: infectiousgastroenteritis
#' 21: melanoma                   22: metabolic_syndrome       
#' 23: NK                         24: otitis                   
#' 25: pneumonia                  26: pre-hypertension         
#' 27: premature_born             28: pyelonefritis            
#' 29: pyelonephritis             30: respiratoryinf           
#' 31: salmonellosis              32: sepsis                   
#' 33: skininf                    34: STEC                     
#' 35: stomatitis                 36: suspinf                  
#'37: T1D                        38: T2D                      
#' 39: tonsillitis                
#' @param TaxonClass: Kindom, Phylum, Class, Order, Family, Genus, Species, Strain
#' @param MicroName: the microbes name: "Proteobacteria","Firmicutes","Bacteroides"
#' @param meta: body_sites including: stool, oral, skin
#' 
#' DiseaseQuery<-c("CRC_VitC","HC_PD1","CRC_PD1","HC")
#' TaxonClass<-"Species"
#' MicroName<-c("Fusobacterium","Neisseria","Akk")
#' meta<-c("stool","oral")
#' School-Age (6 to 12 Years)
microbiomeAnalysor(TaxonClass ="Species",DiseaseQuery = c("control","CRC"),MicroName=c("Fusobacterium","Neisseria"), meta=c("stool","oral"))
microbiomeAnalysor<-function(TaxonClass,DiseaseQuery,
                             MicroName,
                             meta = sitesName|NULL|allsites){
  source("../R/requirements.R")
  expr<-fread("../../database/ExMicroSy/TuMicroSy/allOTUtable_split.csv",data.table = F)
  data<-fread("../../database/ExMicroSy/TuMicroSy/combined_metadata.csv")[,1:10]
  data<-data[!duplicated(data$sampleID),]
  meta<-subset(data,body_site%in%meta)
  meta_data<- subset(meta,Disease%in%DiseaseQuery)
  ID<-meta_data$sampleID
  Index<-c(TaxonClass,ID)
  if(sum(Index%in%colnames(expr)==F)==0){
    message("Query ID check passed !!!")
    newData<-expr[Index]
  } else {
    message(paste0("Missed sampleID: ",setdiff(Index,colnames(expr))))
    message(paste0("dropped ",length(setdiff(ID,colnames(expr)))," unmatched samples"))
    ID1<-intersect(ID,colnames(expr))
    Index<-c(TaxonClass,ID1)
    newData<-expr[Index]
  }
  newData_Melt<-melt(newData,id.vars=TaxonClass,
                     value.name = "Relative_abundance", 
                     variable.name = "sampleID")%>%
    filter(Relative_abundance!=0)%>%
    filter(.[,1]!="Super_Class")
  colnames(newData_Melt)[1]<-"className"
  meta_data[1:3,1:3] 
  newData_Melt[1:3,1:3]
  newData_Melt_merge<-merge(meta_data,newData_Melt,by="sampleID")%>%
    group_by(dataset_name,sampleID,body_site,Disease,
             className,dataset_name,age_category,country,
             Age,BMI,antibiotics_current_use)%>%
    summarise(Relative_abundance=mean(Relative_abundance))
  MicroNameHit<-newData_Melt_merge$className[grep(paste(MicroName,collapse = "|"),newData_Melt_merge$className)]
  message("There are ",length(levels(factor(MicroNameHit)))," ",paste(MicroName,collapse = "|")," !")
  plotdata<-newData_Melt_merge %>% subset(.,.[,5]==MicroNameHit)
  message("plotDataMake finished!")
  message("Plotting ... ... ")
  plotdata<-plotdata[!duplicated(plotdata$sampleID),]
  Freq1<-ggplot(plotdata, 
                 aes(Relative_abundance,body_site, fill = ..density..)) + 
      geom_density_ridges_gradient(aes(height = ..density..),scale = 1,size = 0.3)+
      theme_minimal(base_size = 12)+
      scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'Spectral')))(32))+
      labs(x="Relative_abundance (%)",y="Micobes")+
      ggtitle(paste0("Frequence of ",paste(MicroName,collapse = "|")))+
      theme(plot.title = element_text(hjust = 0.5))+
      facet_wrap(~Disease,scales = "free_x")
  plotdata$sampleID<-plotdata$sampleID[order(plotdata$Disease)]
  heatmap_data<-acast(plotdata[,c(2,5,12)],className~sampleID)
  anno_row<-data.frame(row.names = plotdata$sampleID,Disease=plotdata$Disease)
  heatmap_data[is.na(heatmap_data)]=0
  pheatmap(log(heatmap_data+1),annotation_col = anno_row,border_color = "white",
           color = colorRampPalette(c("white","navy","firebrick3"))(50),
             cluster_cols = F,
             show_colnames = F)
  Heatmap(log(heatmap_data+1))
    Freq2<-ggplot(plotdata, 
                  aes(Relative_abundance,className, fill = ..density..)) + 
      geom_density_ridges_gradient(aes(height = ..density..),scale = 1,size = 0.3)+#sacle���÷�???
      theme_minimal(base_size = 14)+
      scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'Spectral')))(32))+
      labs(x="Relative_abundance (%)",y="Micobes")+
      ggtitle(paste0("Frequence of ",paste(MicroName,collapse = "|")))+
      theme(plot.title = element_text(hjust = 0.5))+
      facet_wrap(body_site~Disease,scales = "free_x")
    
    average1<-plotdata%>%group_by(body_site,Disease,className)%>% summarise(Average=mean(Relative_abundance))
    average2<-plotdata%>%group_by(BMI,Disease,className)%>% summarise(Average=mean(Relative_abundance))
    average3<-plotdata%>%group_by(body_site,age_category,className)%>% summarise(Average=mean(Relative_abundance))
    my_tag1 <- paste0("n = ",summary(factor(plotdata$body_site)))
    av1<-ggplot(average1,
               aes(body_site,Average,fill=className))+
      geom_bar(stat = "identity",width = 1)+
      theme_minimal(base_size = 12)+
      scale_fill_manual(values = c("#303841","#D72323","#377F5B","#6B4F36","#00B8A9","#f0027f",
                          "#FAF8DE","#666666","#BDF1F6","#023782","#5e4fa2","#F1C40F",
                          "#ff7f00","#cab2d6","#240041","#ffff99","#0E3BF0","#a65628",
                          "#f781bf","#808FA6","#2EB872","#F0FFE1","#F33535","#011F4E",
                          "#82B269","#D3C13E","#3F9DCD","#014E1F","#AFFFDF","#3D002E",
                          "#3A554A","#fff2ae","#f1e2cc","#cccccc"),
                        guide_legend(title="Microbes"))+
      facet_wrap(~body_site,scales = "free_x")+
      labs(y="Relative_abundance (Average%)")+
      ggtitle(paste("Comparision of ",paste(MicroName,collapse = "|")))+
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x =  element_blank(),
            axis.title.x =  element_blank())
    av1_1<- tag_facet(av1+theme_presentation(base_size = 12),
                      x = -Inf, y = -Inf,
                      vjust = -0.5, hjust = -1,
                      fontface = 1,#bold
                      size =4,
                      open = "", close = "",
                      tag_pool = my_tag1)
    my_tag2 <- paste0("n = ",summary(factor(plotdata$Disease)))
    av2<-ggplot(average2,aes(Disease,Average,fill=className))+
      geom_bar(stat = "identity",width = 1)+
      theme_minimal(base_size = 12)+
      scale_fill_manual(values = c("#303841","#D72323","#377F5B","#6B4F36","#00B8A9","#f0027f",
                                   "#FAF8DE","#666666","#BDF1F6","#023782","#5e4fa2","#F1C40F",
                                   "#ff7f00","#cab2d6","#240041","#ffff99","#0E3BF0","#a65628",
                                   "#f781bf","#808FA6","#2EB872","#F0FFE1","#F33535","#011F4E",
                                   "#82B269","#D3C13E","#3F9DCD","#014E1F","#AFFFDF","#3D002E",
                                   "#3A554A","#fff2ae","#f1e2cc","#cccccc"))+
      facet_wrap(~Disease,scales = "free_x")+
      labs(y="Relative_abundance (Average%)")+
        theme_presentation(base_size = 12)+
        theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none",
            axis.title.x = element_blank())
    av2_1<- tag_facet(av2, x = -Inf, y = -Inf,
                      vjust = -0.5, hjust = -0.5,
                      fontface = 1,#bold
                      size =4,
                      open = "", close = "",
                      tag_pool = my_tag2)
    average3<-subset(average3,age_category!="NA")
    average3$age_category<-factor(average3$age_category,levels = 
                                    c("newborn","child","schoolage","adult","senior"))
    average3$className<-reorder(average3$className,average3$Average)
    av3<-ggplot(average3,aes(className,Average,fill=age_category))+
      geom_bar(stat = "identity") + 
      theme_ridges()+
      coord_polar()+
      theme(axis.text.x = element_text(size = 8,angle = 15,hjust = 1),
            axis.text.y = element_blank(),
            axis.title = element_blank())+
      scale_fill_manual(values = c("#BDF1F6","#061296","#AAAEDB",
                                   "#F9ED69","#69EDF9","#FFE2E1"))+
      facet_wrap(~body_site)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1,1))) 
    vplayout <- function(x,y){
      viewport(layout.pos.row = x, layout.pos.col = y)
    }
    print(Freq2, vp = vplayout(1,1))
    ggsave("./Plots/Frequence.pdf",
           scale = 1, width = 20, height = 30, units =  "cm",
           dpi = 300, limitsize = TRUE)
    ggsave("./Plots/Frequence.png",
           scale = 1, width = 20, height = 30, units =  "cm",
           dpi = 300, limitsize = TRUE)
    message("FrequencePlot finish! (in Plots folder)")
    if(sum(summary(average3$age_category))!=0){
      pdf(file="./Plots/CompostionPlot.pdf",width = 12,height = 12)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(5,3))) 
      vplayout <- function(x,y){
        viewport(layout.pos.row = x, layout.pos.col = y)
      }
    print(av2_1, vp = vplayout(1:2,1))
    print(av1_1, vp = vplayout(1:2,2:3))
    print(Freq1, vp = vplayout(3,1:3))
    print(av3, vp = vplayout(4:5,1:3))
    dev.off()
    png(file="./Plots/CompostionPlot.png",width = 1000,height = 1000)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(5,3))) 
    vplayout <- function(x,y){
      viewport(layout.pos.row = x, layout.pos.col = y)
    }
    print(av2_1, vp = vplayout(1:2,1))
    print(av1_1, vp = vplayout(1:2,2:3))
    print(Freq1, vp = vplayout(3,1:3))
    print(av3, vp = vplayout(4:5,1:3))
    dev.off()
  message("CompostionPlot finish! (in Plots folder)")
    }else { 
      pdf(file="./Plots/CompostionPlot.pdf",width = 12,height = 12)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(3,3))) 
      vplayout <- function(x,y){
        viewport(layout.pos.row = x, layout.pos.col = y)
      }
      print(av2_1, vp = vplayout(1:2,1))
      print(av1_1, vp = vplayout(1:2,2:3))
      print(Freq1, vp = vplayout(3,1:3))
      dev.off()
      png(file="./Plots/CompostionPlot.png",width = 1000,height = 1000)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(3,3))) 
      vplayout <- function(x,y){
        viewport(layout.pos.row = x, layout.pos.col = y)
      }
      print(av2_1, vp = vplayout(1:2,1))
      print(av1_1, vp = vplayout(1:2,2:3))
      print(Freq1, vp = vplayout(3,1:3))
      dev.off()
      message("CompostionPlot finish! (in Plots folder)")
    }
  return(list(plotdata,average1,average2,average3))
}
