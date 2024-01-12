# Shiny app test

library(shiny)
library(shinyWidgets)
library(ggplot2)
library(gplots)
library(cowplot)
library(gridExtra)
library(grid)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyverse)
library(stringr)
library(data.table)
library(scales)
library(DT)
library(ggbeeswarm)
library(sp)
library(raster)
library(lattice)

#Load data from separate gradients

table.ctrl1<-read.table("table.ctrl1.csv",sep="\t",header=TRUE)

table.ctrl2<-read.table("table.ctrl2.csv",sep="\t",header=TRUE)

table.ctrl3<-read.table("table.ctrl3.csv",sep="\t",header=TRUE)

table.rnase1<-read.table("table.rnase1.csv",sep="\t",header=TRUE)

table.rnase2<-read.table("table.rnase2.csv",sep="\t",header=TRUE)

table.rnase3<-read.table("table.rnase3.csv",sep="\t",header=TRUE)


#Load the mean and sd

ctrl_norm_mean<-read.table("ctrl_norm_mean.csv",header=TRUE,sep="\t")
rownames(ctrl_norm_mean)<-ctrl_norm_mean$new_Gene

rnase_norm_mean<-read.table("rnase_norm_mean.csv",header=TRUE,sep="\t")
rownames(rnase_norm_mean)<-rnase_norm_mean$new_Gene

ctrl_norm_sd<-read.table("table_sd_ctrl_combined_proportions_annotation.csv",sep="\t",header=TRUE)
rownames(ctrl_norm_sd)<-ctrl_norm_sd$new_Gene

rnase_norm_sd<-read.table("table_sd_rnase_combined_proportions_annotation.csv",sep="\t",header=TRUE)
rownames(rnase_norm_sd)<-rnase_norm_sd$new_Gene



## Load SVM data

table_SVM<-read.table("table_SVM_plus_annotations.csv",sep="\t",header=TRUE)

levels.for.factor<-rev(c("Arthrospira platensis NIES-39",
                         "Crinalium epipsammum PCC 9333",
                         "Geitlerinema sp. PCC 7407",
                         "Leptolyngbya sp. PCC 7376",
                         "Microcoleus sp. PCC 7113",
                         "Oscillatoria acuminata PCC 6304",
                         "Oscillatoria nigro-viridis PCC 7112",
                         "Pseudanabaena sp. PCC 7367",
                         "Spirulina major PCC 6313",
                         "Trichodesmium erythraeum IMS101",
                         "Fischerella sp. NIES-3754",
                         "Fischerella sp. NIES-4106",
                         "Anabaena cylindrica PCC 7122",
                         "Calothrix parietina PCC 6303",
                         "Calothrix sp. PCC 7507",
                         "Cylindrospermopsis raciborskii Cr2010",
                         "Cylindrospermum stagnale PCC 7417",
                         "Microchaete diplosiphon NIES-3275",
                         "Nodularia spumigena CCY9414",
                         "Nostoc punctiforme PCC 73102",
                         "Nostoc sp. PCC 7107",
                         "Nostoc sp. PCC 7120",
                         "Nostoc sp. PCC 7524",
                         "Rivularia sp. PCC 7116",
                         "Tolypothrix sp. PCC 7910",
                         "Trichormus azollae 0708",
                         "Trichormus variabilis ATCC 29413",
                         "Chroococcidiopsis thermalis PCC 7203",
                         "Cyanothece sp. PCC 7425",
                         "Stanieria cyanosphaera PCC 7437",
                         "Pleurocapsa minor PCC 7327",
                         "Acaryochloris marina MBIC11017",
                         "Chamaesiphon minutus PCC 6605",
                         "Cyanobacterium aponinum PCC 10605",
                         "Cyanobacterium stanieri PCC 7202",
                         "Cyanobium gracile PCC 6307",
                         "Dactylococcopsis salina PCC 8305",
                         "Gloeobacter violaceus PCC 7421",
                         "Gloeocapsa sp. PCC 7428",
                         "Microcystis aeruginosa NIES-843",
                         "Microcystis aeruginosa PCC 7806SL",
                         "Prochlorococcus marinus MIT9211",
                         "Prochlorococcus marinus MIT9303",
                         "Prochlorococcus marinus MIT9313",
                         "Prochlorococcus marinus NATL1A",
                         "Prochlorococcus marinus pastoris CCMP 1986",
                         "Synechococcus elongatus PCC 6301",
                         "Synechococcus elongatus PCC 7942",
                         "Synechococcus sp. CC9311",
                         "Synechococcus sp. JA-2-3B'a(2-13)",
                         "Synechococcus sp. JA-3-3Ab",
                         "Synechococcus sp. WH 8103 WH8103",
                         "Synechococcus sp. WH 8109",
                         "Synechocystis sp. PCC 6803",
                         "Thermosynechococcus elongatus BP-1"))


# Load annotation and shift data

table_annotation<-read.table("table_annotation_proteins.csv",header=TRUE,sep="\t") # I have annotated manually the proteins of Refseq and not present in UNIPROT
rownames(table_annotation)<-table_annotation$new_Gene

# Load table clustering data

factor_soft_cluster<-table_annotation$Cluster
names(factor_soft_cluster)<-table_annotation$new_Gene

factor_soft_cluster<-factor_soft_cluster[!is.na(factor_soft_cluster)]
vector_soft_cluster<-paste0("Cluster ",factor_soft_cluster)
names(vector_soft_cluster)<-names(factor_soft_cluster)

table_clustering<-ctrl_norm_mean[names(vector_soft_cluster),2:19]

table_clustering$Cluster<-vector_soft_cluster
table_clustering$Gene<-rownames(table_clustering)

# Load table with individual gradients

table_individual_gradients<-read.table("table_expression_combined_proportions_annotation.csv",header=TRUE,sep="\t")
table_individual_gradients<-table_individual_gradients[,c(1:19,21)]

###### ui

  ui<-fluidPage(title="GradR Nostoc",fluidRow(
                column(width=4),
                column(h1("GradR", em("Nostoc"), style="font-family: 'Arial'; font-weight: bolder; font-size: 50px; color: #006594",align="center"),width=4),
                column(img(src='CyanoRBP.png', height="80%", width="80%", align = "right"),width=4)),
              tabsetPanel(id="tabset",
                tabPanel(h4("GradR",style="color: #006594"),fluidRow(h2("")),
                         fluidRow(
                           column(wellPanel(p(h3(em("Nostoc"), "Protein of Interest"),br()),
                                            radioButtons(inputId = "Search_mode",  choices = c("All", "Shifting proteins"), label="Search for all proteins or only shifting proteins:", selected = c("All"), inline=TRUE),
                                            pickerInput(inputId = "Gene_ID_GradR", label = "Please select a protein: ", choices = rownames(table.ctrl1), selected=c("rplA(alr5301)"), 
                                                        options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple = F),
                                            downloadButton("downloadGradRPlot", "Download GradR Plot"),
                                            downloadButton("downloadGradRTable", "Download table with distribution per gradient")
                                            ),
                                  wellPanel(p(h4(textOutput("RNA_Dependance")),br(),textOutput("High_Decision"),textOutput("Low_Decision")),
                                            conditionalPanel(condition= "output.RNA_Dependance == 'RNA-dependent shift'",
                                                             conditionalPanel(condition="output.High_Decision == ' '",p("Fraction hits with shift (foldchange > 2 and fdr < 0.05): ",textOutput("High_significant"))),
                                                             conditionalPanel(condition="output.Low_Decision == ' '",p("Fraction candidates with shift (foldchange > 1.5 and fdr < 0.02): ",textOutput("Low_significant")))
                                            )
                                            ),
                                  wellPanel(p(h4(textOutput("SVM_prediction")),br()),
                                            conditionalPanel(condition="output.SVM_prediction == 'Predicted as RBP'",p(style="text-align: justify;","This protein is predicted as an RBP by our SVM approach. 
                                                            To check its score and the conservation of its RNA-binding capacity, click here:"),
                                                             actionButton(inputId = "SVMbutton",label="Go to SVM Prediction"))
                                    
                                  )
                                  ,width = 3),
                           column(wellPanel(uiOutput("GradRPlot_UI"),style="background: white"),width = 6, style="border: grey"),
                           column(wellPanel(p(h3(textOutput("Description"),align="center"),br()),
                                            wellPanel(fluidRow(column(h5("Protein Accesion:",align="left"),width=6),
                                                               column(h5(textOutput("P_accesion"),align="right"),width=6)),
                                                      
                                                      fluidRow(column(h5("Molecular Weight:",align="left"),width=6),
                                                               column(h5(textOutput("MW"),align="right"),width=6)),
                                                      
                                                      fluidRow(column(h5("Length:",align="left"),width=6),
                                                               column(h5(textOutput("AA_Length"),align="right"),width=6)),
                                                      
                                                      fluidRow(column(h5("Isoelectric Point:",align="left"),width=6),
                                                               column(h5(textOutput("PI"),align="right"),width=6)),
                                            style="background-color: #ffffff",align="center"),
                                            p(h5("Link to protein database:",style="font-family: 'Arial', color: #006594"),uiOutput("Link_UNIPROT"))
                                          ),
                                  wellPanel(p(h4("Clustering"),textOutput("Clustering_decision"),br()),
                                            conditionalPanel(condition="output.Clustering_decision == ' '",p("This protein co-fractionate with
                                            other proteins in Cluster",textOutput("Clustering",inline=T),". To look for potential interaction partners click here to go to Cluster visualization: "),
                                                             actionButton(inputId = "Clusteringbutton",label="Test")),
                                            conditionalPanel(condition="output.Clustering_decision != ' '",p(style="text-align: justify;","This protein did not show a clear co-fractionation pattern with other potential interaction partners"))
                                            
                                  )
                           ,width = 3)
                           )
                         ),
                tabPanel(h4("Co-fractionation Explorer",style="color: #006594"),fluidRow(h2("")),
                         fluidRow(
                           column(wellPanel(p(h3("Co-fractionation Explorer"),br()),
                                        pickerInput(inputId = "Gene_ID_m", label = "Please select your proteins of interest: ",
                                        choices = rownames(table.ctrl1), selected=c("rplA(alr5301)"), 
                                        options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple = T),
                                        radioButtons(inputId = "Sample",  choices = c("Ctrl", "RNases"), label="Sample", selected = c("Ctrl"), inline=TRUE),
                                        selectInput(inputId = "height_HeatmapMul", label = "Please adjust Plot Height: ",
                                                    choices = c(1:10*10), selected=20),
                                        sliderInput(inputId = "Heatmap_width", label = "Please select Plot width for download:",
                                                    min = 1, max = 100, value = 20,  step = 5 ),
                                        downloadButton("downloadHeatmap", "Download Heatmap")
                                            ),width = 3),
                           column(uiOutput("GradSeqHeatmapMul_summaryUI"),width = 9)
                         )
                         ),
                tabPanel(h4("Clustering",style="color: #006594"),value="Clustering",fluidRow(h2("")),
                         fluidRow(
                           column(wellPanel(p(h3("Clustering analysis (ctrl samples)"),br()),
                                            selectInput(inputId = "Cluster_analysis", label = "Please select Cluster: ",
                                                        choices = paste("Cluster ",1:19,sep=""), selected="Cluster 19"),
                                            downloadButton("downloadClusterPlot", "Download Cluster Plot"),
                                            downloadButton("downloadClusterTable", "Download Cluster Table")
                           ),
                           uiOutput("ClusterlinePlotUI"),h6("Average distribution of proteins in", textOutput("cluster_caption",inline=T),"along the gradient. Grey lines show the distribution of individual proteins. Blue line show the average distribution of the cluster."),
                           width = 4),
                           column(wellPanel(DT::dataTableOutput("table_Cluster"),style="background-color: white")
                             ,width = 8)
                           
                         )
                         ),
                tabPanel(h4("SVM",style="color: #006594"),value="SVM",fluidRow(h2("")),
                         fluidRow(
                           column(wellPanel(p(h3("Sidebar SVM Explorer"),br()),
                                              pickerInput(inputId = "Gene_SVM_ID", label = "Please select a protein: ",choices = unique(table_SVM$new_Gene), 
                                                          selected=c("rplA(alr5301)"), options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple = T),
                                            radioButtons(inputId = "Plot_type",  choices = c("Boxplot", "Heatmap"), label="Please select the Plot type:", selected = c("Boxplot"), inline=TRUE),
                                            selectInput(inputId = "height_SVM", label = "Please adjust Plot Height: ",
                                                          choices = c(1:10*10), selected=30),
                                            sliderInput(inputId = "SVM_width", label = "Please select Plot width for download:",
                                                        min = 1, max = 100, value = 20,  step = 5 ),
                                            downloadButton("downloadSVMplot", "Download SVM plot")
                                            ),width = 3),
                           column(wellPanel(uiOutput("SVM_UI"),
                                            conditionalPanel(condition="input.Plot_type == 'Boxplot'",h6("SVM score of homologs of selected proteins. Red dot: Score for Nostoc proteins. Black dots: Score for homologs. Red dashed line: threshold for considering a RBP. Scores higher than 0.25 are considered putative RBPs.")),
                                            conditionalPanel(condition="input.Plot_type == 'Heatmap'",h6("SVM score of homologs of selected proteins in different cyanobacterial strains. No found homologs are shown in grey. Red color show higher SVM score which indicate a putative higher capacity to bind RNA.")),
                                            style="background: white; border: grey")
                                  ,width = 9)
                         )
                         ),
                tabPanel(h4("Tutorial",style="color: #006594"),value="Tutorial",fluidRow(h2("")),
                         h4("Introduction"),
                         p("This is a comprehensive database for the identification of putative RNA-binding proteins (RBPs) in the multicellular cyanobacterium", em("Nostoc"),"sp. PCC 7120. We have used GradR/R-DeeP to detect the RNA-dependent proteome of this cyanobacterium. These methods rely on the fractionation of the RNA-protein or protein-protein complexes in sucrose gradients. After the treatment of a gradient with RNAses, if the apparent distribution of a protein changes, that would point out to the participation of the protein in a big RNA-protein complex (see Figure 1).",
                           br(),br(),
                           "The apparent shift of a protein in the gradients could not be due to their RNA-binding capacity. To filter out proteins that are bound to other RBPs (RNA-dependent proteins) from direct RNA-binding proteins, we have implemented a support vector machine (SVM) based on TriPepSVM approach. The overlay between the proteins shifting and the proteins with good SVM scores have yielded a number of potential RBPs.",br(),br()),
                          img(src='Fig1.png', height="45%", width="45%", style="display: block; margin-left: auto; margin-right: auto;"),br(),br(),
                           p("Finally, the clustering analysis of the control gradients has allowed us to predict new components of big macromolecular complexes. The quality of the data is shown in the co-fractionation of 50S and 30S ribosomal subunits, PSI and PSII core and components of the RNA polymerase.",br(),br()),
                         h4("GradR tab"),
                         p("This is the main page of the app. From here, the user can access all the information we have for a selected protein. The user can select all the detected proteins by MS or only the ones showing a significant shift. For each selected protein it is displayed if a protein has a significant shift, if a protein was predicted by our SVM approach, some general information and the clustering results.",
                           br(),br(),
                           "We have used ",em("limma")," to test if the distribution of a protein in each fraction differs between three control gradients or three RNAse-treated gradients.  If the user selects a protein, the app will show if there are any fractions with significant changes. The user can check the distribution of the protein in the control and RNAse-treated gradients in the central plot.",
                           br(),br(),
                           "If the protein was predicted to be an RBP by our SVM approach that will be displayed in a different panel. If the user clicks on the button, he will go directly to the 'SVM Score' section for that protein.",
                           br(),br(),
                           "At the same time, if the protein was clustered, the clustering information will be displayed below the general information for the protein. If the user clicks on the button, he will go directly to the particular cluster in which the user can check potential interaction partner that had the same sedimentation profile."),
                         img(src='Fig_S1.png', height="60%", width="60%", style="display: block; margin-left: auto; margin-right: auto;"), br(), br(),
                         h4("Co-fractionation tab"),
                         p("The user can select any number of proteins to test their sedimentation profiles along the gradient. It can be selected which gradients should be displayed, either the control gradients or the RNAse treated gradients. The selection of height and width for downloading the heatmap plot is necessary to get a proper plot.",br(),br()),
                         img(src='Fig_S2.png', height="60%", width="60%", style="display: block; margin-left: auto; margin-right: auto;"), br(), br(),
                         h4("Clustering tab"),
                         p("Our approach has predicted 19 distinct cluster of protein. The user can select any cluster. The sedimentation profile for the cluster will be plotted and a table with the components of the cluster are displayed and ready for downloading.",br(),br()),
                         img(src='Fig_S3.png', height="60%", width="60%", style="display: block; margin-left: auto; margin-right: auto;"), br(), br(),
                         h4("SVM tab"),
                         p("This tab shows the results of our SVM approach. Any number of protein can be selected and the SVM scores of that protein in ", em("Nostoc"), "and their cyanobacteria homologs will be plotted. The user can select either a Boxplot or a Heatmap plot. The boxplot is recommended for a first check of how good the RNA-binding capacity is conserved in cyanobacteria. The heatmap allow the user a better analysis of which particular cyanobacterial homologs shows good SVM score and therefore a putative RNA-binding capacity.",br(),br()),
                         img(src='Fig_S4.png', height="60%", width="60%", style="display: block; margin-left: auto; margin-right: auto;"), br(), br()
                         )
              ),
              div(class="footer", tags$style(type="text/css", "body {padding-bottom: 70px;}")),
              tags$footer(strong("Brenes-Alvarez et al., unpublished"), style = "
                 text-align: center;
                 font-size: 14px;
                 font-family: Arial;
                 color: #006594;
                 position:fixed;
                 bottom: 0px;
                 width:100%;
                 height:50px;
                 border: 1px solid;
                 border-color: #e7e7e7;
                 margin-bottom: 0;
                 background-color: #f8f8f8;
                 z-index: 1030;
                 padding: 15px;
                ")

                            )


server<-function(input,output,session)
{
  GradR<-reactive({
    pos<-input$Gene_ID_GradR
    
    x <- c(1:18)
    ctrl <- as.numeric(ctrl_norm_mean[pos,2:19])
    title <- pos
    df_ctrl <- data.frame(x, ctrl)
    colnames(df_ctrl) <- c("Fractions","Amount")
    rnase <- as.numeric(rnase_norm_mean[pos,2:19])
    df_rnase <- data.frame(x, rnase)
    colnames(df_rnase) <- c("Fractions","Amount")
    
    # standard deviation of mean values
    df_sd_ctrl <- df_ctrl
    df_sd_ctrl$upper <- ctrl + as.numeric(ctrl_norm_sd[pos,2:19])
    df_sd_ctrl$lower <- ctrl - as.numeric(ctrl_norm_sd[pos,2:19])
    
    df_sd_rnase <- df_rnase
    df_sd_rnase$upper <- rnase + as.numeric(rnase_norm_sd[pos,2:19])
    df_sd_rnase$lower <- rnase - as.numeric(rnase_norm_sd[pos,2:19])
    
    scatterPlot <- ggplot() +
      
      # curve CTRL
      geom_line(data=df_ctrl, aes(x=Fractions, y=Amount, color="green"), size=.5) +
      geom_point(data=df_ctrl, mapping=aes(x=Fractions, y=Amount), color="green", shape=1, size=1) +
      
      # shadow CTRL
      geom_area(data=df_ctrl, aes(x=Fractions, y=Amount), fill="green", alpha=.1, show.legend=F) +
      
      # shadded errors
      geom_ribbon(data=df_sd_ctrl,aes(x=Fractions,ymin=lower,ymax=upper),alpha=0.3) +
      
      # curve RNASE
      geom_line(data=df_rnase, aes(x=Fractions, y=Amount, color="red"), size=.5) +
      geom_point(data=df_rnase, mapping=aes(x=Fractions, y=Amount), color="red", shape=2, size=1) +
      
      # shadded errors
      geom_ribbon(data=df_sd_rnase,aes(x=Fractions,ymin=lower,ymax=upper),alpha=0.3) +
      
      # shadow RNASE
      geom_area(data=df_rnase, aes(x=Fractions, y=Amount), fill="red", alpha=.1, show.legend=F) +
      
      # New axis label
      ylab("Normalized protein amount (%)") +
      
      # title
      labs(title = title) +
      
      # Legend
      scale_color_manual(values = c("#09b957","#f42730"), labels=c("Control","RNases")) +
      theme(legend.position=c(0.75,0.98),
            legend.text=element_text(size=12),
            legend.direction="horizontal",
            legend.title.align = 0,
            legend.justification=c(0,1),
            legend.title=element_blank(),
            axis.title.x=element_blank(),
            plot.title = element_text(size = rel(2),
                                      colour = "black",
                                      hjust = 0.5, face = "bold"),
            legend.background = element_rect(colour = "black"),
            legend.key = element_rect(fill = "white")) +
      
      # panel options
      theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(colour = "lightgrey"), panel.grid.minor = element_line(colour = "lightgrey"), axis.title.y = element_text(color="black",size=12,face="bold"), axis.text.x=element_text(size=10),plot.margin = margin(t=1, r=0, b=1, l=0, "cm"))  +
      
      # x-axis scale and ticks
      scale_x_discrete(name ="Fractions", limits=c("1_2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"))
    
    
    # plot a raster with ggplot
    dataf1 <- data.frame(ctrl_norm_mean[pos,2:19],rnase_norm_mean[pos,2:19])
    
    r1 <- raster(xmn=0, xmx = 18, ymn = 0, ymx = 1, nrows = 1, ncols = 18)
    r2 <- raster(xmn=0, xmx = 18, ymn = 0, ymx = 1, nrows = 1, ncols = 18)
    
    r1[] <- as.numeric(dataf1[1,1:18])
    r2[] <- as.numeric(dataf1[1,19:36])
    r.spdf1 <- as(r1, "SpatialPixelsDataFrame")
    r1.df <- as.data.frame(r.spdf1)
    colnames(r1.df) <- c("Amount","Fraction","ctrl")
    r.spdf2 <- as(r2, "SpatialPixelsDataFrame")
    r2.df <- as.data.frame(r.spdf2)
    colnames(r2.df) <- c("Amount","Fractions","RNases")
    
    # now plot the whole
    plot1 <- ggplot(r1.df, aes(x=Fraction+0.5, y=ctrl)) + geom_tile(aes(fill = Amount)) + coord_equal(2) + theme(legend.position="none", axis.title.y = element_text(color="black",size=12,face="bold"),axis.text.x=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_fill_gradient(low="white", high="green",n.breaks=100) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank())  + scale_x_discrete(name ="Fractions", limits=c("1_2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"))
    plot2 <- ggplot(r2.df, aes(x=Fractions+0.5, y=RNases)) + geom_tile(aes(fill = Amount)) + coord_equal(2) + theme(legend.position="none",axis.text.x=element_text(size=10),axis.title.y = element_text(color="black",size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_fill_gradient(low="white", high="red",n.breaks=100) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank()) + scale_x_discrete(name ="Fractions", limits=c("1_2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"))
    
    # combine all graphics in one panel
    grid.arrange(scatterPlot, plot1, plot2, ncol=1, nrow=3, heights=c(1.5,.5,.5), bottom = textGrob("Fractions", gp = gpar(fontsize = 12,fontface="bold")))
    
  })
  
  extract_table_GradR<- reactive({
    pos <- pos<-input$Gene_ID_GradR
    table_GradR<-table_individual_gradients %>% filter( new_Gene %in% pos )
    table_GradR<-table_GradR[,c(1:19)]
    return(table_GradR)
  })
  
  get_Description<-reactive({
    pos<-input$Gene_ID_GradR
    table_annotation[pos,"description"]
  })
  
  get_ID<-reactive({
    pos<-input$Gene_ID_GradR
    if(is.na(table_annotation[pos,"UNIPROT"]))
    {
      table_annotation[pos,"PROTEIN_REFSEQ"] 
    }
    else
    {
      table_annotation[pos,"UNIPROT"]
    }
  })
  
  get_MW<-reactive({
    pos<-input$Gene_ID_GradR
    table_annotation[pos,"MW"]
  })
  get_AA_length<-reactive({
    pos<-input$Gene_ID_GradR
    table_annotation[pos,"LENGTH_AA"]
  })
  get_PI<-reactive({
    pos<-input$Gene_ID_GradR
    table_annotation[pos,"PI"]
  })
  get_link<-reactive({
    pos<-input$Gene_ID_GradR
    if(is.na(table_annotation[pos,"UNIPROT"]))
    {
      table_annotation[pos,"link_NCBI"] 
    }
    else
    {
      table_annotation[pos,"link_UNIPROT"]
    }
    })
  get_database<-reactive({
    pos<-input$Gene_ID_GradR
    if(is.na(table_annotation[pos,"UNIPROT"]))
    {
      "NCBI" 
    }
    else
    {
      "UNIPROT"
    }
  })
  RNA_dependance_function<-reactive({
    pos<-input$Gene_ID_GradR
    if(table_annotation[pos,"Shift"])
    {
      "RNA-dependent shift" 
    }
    else
    {
      "No significant shift"
    }
  })
  
  hits_decision<-reactive({
    pos<-input$Gene_ID_GradR
    if(!is.na(table_annotation[pos,"High_significant_fractions"]))
    {
     " " 
    }
  })
  
  candidates_decision<-reactive({
    pos<-input$Gene_ID_GradR
    if(!is.na(table_annotation[pos,"Low_significant_fractions"]))
    {
      " " 
    }
  })
  
  extract_hits<-reactive({
    pos<-input$Gene_ID_GradR
    if(!is.na(table_annotation[pos,"High_significant_fractions"]))
      {
      paste(strsplit(table_annotation[pos,"High_significant_fractions"],split=" ")[[1]],collapse=",")
      }
    })
 
  extract_candidates<-reactive({
    pos<-input$Gene_ID_GradR
    if(!is.na(table_annotation[pos,"Low_significant_fractions"]))
    {
      paste(strsplit(table_annotation[pos,"Low_significant_fractions"],split=" ")[[1]],collapse=",")
    }
  })
  
  SVM_prediction_decision<-reactive({
    pos<-input$Gene_ID_GradR
    if(table_annotation[pos,"SVM_prediction"])
    {
      "Predicted as RBP" 
    }
    else
    {
      "No predicted as RBP"
    }
  })
  
  cluster_decision<-reactive({
    pos<-input$Gene_ID_GradR
    if(!is.na(table_annotation[pos,"Cluster"]))
    {
      " "
    }
  })
  
  extract_cluster<-reactive({
    pos<-input$Gene_ID_GradR
    if(!is.na(table_annotation[pos,"Cluster"]))
    {
      table_annotation[pos,"Cluster"]
    }
    else
    {
      " "
    }
  })
  
  HeatmapSize_Mul <- reactive({ req(input$height_HeatmapMul)
    as.numeric(input$height_HeatmapMul)
  })
  HeatmapHeight_Mul <- reactive( (20 * HeatmapSize_Mul()) )
  
  
  GradSeqHeatmap_mul<- reactive({
    
    Gene_ID_in_m <- input$Gene_ID_m
    Sample<-input$Sample
    
    if(Sample == "Ctrl")
    {
      data2<-ctrl_norm_mean
    }
    
    if(Sample == "RNases")
    {
      data2<-rnase_norm_mean
    }

    # Grad-Seq Heatmap Data.frame:
    
    data2 <- data2 %>% filter( new_Gene %in% Gene_ID_in_m ) 
    data2 <-  reshape2::melt(data2, id.vars=c("Gene","Sample","Feature.Subtype","ID","new_Gene"))
    data2$Gene<-factor(x=data2$new_Gene,levels=rev(Gene_ID_in_m),ordered=TRUE)
    
    ggplot(data2, aes(x = variable, y = new_Gene, fill = value)) +  geom_tile()+
      scale_fill_gradient("% Protein:  ",low = "white",high="black") + theme_bw() + labs(x="Fractions",y="Protein") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing.y = unit(-0.1, "lines"), 
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1), legend.position = "bottom",
            text=element_text(size=15,colour="black")) 
  })
  
  cluster_lineplot<-reactive({
    cluster<-input$Cluster_analysis
    table_clustering_filtered <- table_clustering %>% filter( Cluster %in% cluster ) 
    m_clust.data<-reshape2::melt(table_clustering_filtered,id.vars=c("Gene","Cluster"))
    m_clust.data
    
    head(m_clust.data)
    colnames(m_clust.data)<-c("Gene","Cluster","sample","mean.value")
    head(m_clust.data)
    
    customPlot <- list(
      theme_bw(base_size = 14), 
      scale_shape_manual(values = c(16, 17, 15, 3, 7, 8)), 
      scale_fill_brewer(palette = "Set1"), 
      scale_colour_brewer(palette = "Set1")
    )
    
    ggplot(data = m_clust.data, aes(sample, mean.value)) +
      geom_hline(yintercept = 0) +
      geom_line(aes(group = paste(Gene, Cluster)), alpha = 0.1) +
      geom_smooth(fun.data = "mean_se", stat = "summary", aes(group = Cluster)) +
      facet_wrap(~ Cluster) +
      customPlot +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  +
      ylab("% protein") + xlab("Fractions")
    

  })
  
  extract_table_clustering<- reactive({
    cluster <- as.numeric(strsplit(input$Cluster_analysis,split=" ")[[1]][2])
    table_cluster<-table_annotation %>% filter( Cluster %in% cluster )
    table_cluster<-table_cluster[,c("new_Gene","UNIPROT","PROTEIN_REFSEQ","description","MW","LENGTH_AA","Cluster_confidence","Shift")]
    colnames(table_cluster)<-c("Entry","UNIPROT","REFSEQ","Description","Molecular Weight (KDa)", "Length","Cluster confidence (%)","Shift")
    return(table_cluster)
  })
  
  SVM_Ind<-reactive({
    Gene_SVM <- input$Gene_SVM_ID
    plot_type<-input$Plot_type
    if(plot_type == "Boxplot")
    {
      data_SVM <- table_SVM %>% filter( new_Gene %in% Gene_SVM )
      data_SVM$highlight<-rep(FALSE,nrow(data_SVM))
      data_SVM$highlight[data_SVM$Query.Genome == "Nostoc sp. PCC 7120"]<-TRUE
      data_SVM<-data_SVM[,c(8,7,9)]
      
      ggplot(subset(data_SVM, !highlight),aes(x=new_Gene, y=Score, fill=new_Gene)) + geom_boxplot() + stat_boxplot(geom="errorbar") + 
        geom_beeswarm(cex=1, size=1) + geom_hline(yintercept = 0.25, col = "red",lty="dashed") +
        geom_point(data=subset(data_SVM, highlight), aes(x=new_Gene, y=Score), color="red", size=2) +
        theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=90,color="black",size=rel(1.5)),legend.position="none",axis.title.y=element_text(size=rel(1.5)),axis.text.y=element_text(size=rel(1.5)))
      
    }
    else
    {
      data_SVM_X <- table_SVM %>% filter( new_Gene %in% Gene_SVM )
      data_SVM_X<-data_SVM_X[,c("new_Gene","Query.Genome","Score")]
      data_SVM_X2 <- expand.grid(new_Gene=unique(data_SVM_X$new_Gene),Query.Genome=unique(table_SVM$Query.Genome)) %>% left_join(data_SVM_X)
      
      data_SVM_X2 <-  reshape2::melt(data_SVM_X2, id.vars=c("new_Gene","Query.Genome"))
      data_SVM_X2$Query.Genome<-factor(data_SVM_X2$Query.Genome,levels=levels.for.factor)
      
      ggplot(data_SVM_X2, aes(x = new_Gene, y = Query.Genome, fill = value)) +  geom_tile(colour="white",size=0.25) + 
        scale_fill_gradient2("Score  ",low = "blue",mid="white",high="red",na.value="grey80",limits=c(-2,2)) + labs(x="",y="Organism") +
        theme(panel.background=element_blank(),panel.border=element_blank(),legend.background = element_rect(linetype = 2, size = 0.2, colour = 1), legend.position = "bottom",
              axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=8),legend.text=element_text(size=9))
      
    }
    
  })
  
  HeatmapHeight_Mul <- reactive( (20 * HeatmapSize_Mul()) )
  
  
  Size_SVM <- reactive({ req(input$height_SVM)
    as.numeric(input$height_SVM)
  })
  Height_SVM <- reactive({
    if(input$Plot_type == "Boxplot")
    {
      20*20
    }
    else
    {
      20 * Size_SVM()
    }
  })
  
  output$Description <- renderText({get_Description()})
  output$P_accesion<-renderText({get_ID()})
  output$MW<-renderText({paste0(get_MW()," KDa")})
  output$AA_Length<-renderText({paste0(get_AA_length()," AA")})
  output$PI<-renderText({get_PI()})
  output$Link_UNIPROT<-renderUI({HTML(paste0("<a href=",get_link(),"> Link to ",get_database(),"</a>"))})
  output$GradRPlot <- renderPlot({GradR()})
  output$GradRPlot_UI<-renderUI({plotOutput("GradRPlot",height=20*30)})
  output$RNA_Dependance<-renderText({RNA_dependance_function()})
  output$High_Decision<-renderText({hits_decision()})
  output$Low_Decision<-renderText({candidates_decision()})
  output$High_significant<-renderText({extract_hits()})
  output$Low_significant<-renderText({extract_candidates()})
  output$downloadGradRPlot<-downloadHandler(filename= function(){paste0("GradR_",input$Gene_ID_GradR, ".pdf")},
                                           content= function(file){ggsave(file,plot=GradR(),device="pdf",width=10,height=10)}
                                            )
  output$SVM_prediction<-renderText({SVM_prediction_decision()})
  output$Clustering_decision<-renderText({cluster_decision()})
  output$Clustering<-renderText({extract_cluster()})
  
  output$GradSeqHeatmapMul_summary <- renderPlot({GradSeqHeatmap_mul()})
  
  output$GradSeqHeatmapMul_summaryUI <- renderUI({
    plotOutput("GradSeqHeatmapMul_summary", height = HeatmapHeight_Mul())})
  
  Heatmap_Width <- reactive({ req(input$Heatmap_width)
    as.numeric(input$Heatmap_width)
  })
  HeatmapWidth <- reactive( (Heatmap_Width()))
  
  output$downloadHeatmap<-downloadHandler(filename= "Co_fractionation_heatmap.pdf",
                                            content= function(file){
                                              pdf(file, width = HeatmapWidth(), height = HeatmapHeight_Mul()/50)
                                              print(GradSeqHeatmap_mul())
                                              dev.off()
                                            }
                                          )
  output$ClusterlinePlot<-renderPlot({cluster_lineplot()})
  output$ClusterlinePlotUI<-renderUI({plotOutput("ClusterlinePlot")})
  output$downloadClusterPlot<-downloadHandler(filename= function(){paste0(input$Cluster_analysis,"_lineplot.pdf")},
                                            content= function(file){ggsave(file,plot=cluster_lineplot(),device="pdf",width=10,height=10)}
  )
  output$table_Cluster<-DT::renderDataTable(extract_table_clustering(),options=(list(scrollX=FALSE,paging=TRUE)),rownames=FALSE)
  output$downloadClusterTable<-downloadHandler(filename= function(){paste0(input$Cluster_analysis,"_table.txt")},
                                              content= function(file){write.table(extract_table_clustering(), file,sep="\t",row.names=FALSE)}
  )
  output$SVM_plot<-renderPlot({SVM_Ind()})
  output$SVM_UI<-renderUI({plotOutput("SVM_plot", height = Height_SVM())})

  observeEvent(input$SVMbutton, {
    updatePickerInput(session = session, inputId="Gene_SVM_ID",selected =input$Gene_ID_GradR)
    updateTabsetPanel(session = session, inputId = "tabset", selected = "SVM")
  })
  
  observeEvent(input$Gene_ID_GradR, {
    updateActionButton(session = session, inputId="Clusteringbutton",label=paste0("Search Cluster ",extract_cluster()) )
  })
  
  observeEvent(input$Clusteringbutton, {
    updatePickerInput(session = session, inputId="Cluster_analysis",selected = paste0("Cluster ",extract_cluster()) )
    updateTabsetPanel(session = session, inputId = "tabset", selected = "Clustering")
  })
  
  output$cluster_caption<-renderText({input$Cluster_analysis})
  
  SVM_Width <- reactive({ req(input$SVM_width)
    as.numeric(input$SVM_width)
  })
  SVMWidth <- reactive( (SVM_Width()))
  
  output$downloadSVMplot<-downloadHandler(filename= paste0("SVM_",input$Plot_type,".pdf"),
                                          content= function(file){
                                            pdf(file, width = SVMWidth(), height = Height_SVM()/50)
                                            print(SVM_Ind())
                                            dev.off()
                                          }
  )
  
  output$downloadGradRTable<-downloadHandler(filename= function(){paste0(input$Gene_ID_m,"_table.txt")},
                                               content= function(file){write.table(extract_table_GradR(), file,sep="\t",row.names=FALSE)}
  )
  observeEvent(input$Search_mode, {
    if(input$Search_mode == "All")
    {
      updatePickerInput(session = session, inputId="Gene_ID_GradR",choices= rownames(table.ctrl1),selected=c("rplA(alr5301)"))
    }
    
    if(input$Search_mode == "Shifting proteins")
    {
      updatePickerInput(session = session, inputId="Gene_ID_GradR",choices= rownames(table_annotation)[table_annotation$Shift],selected=c("rplA(alr5301)"))
    }
  })
  
}

#    

shinyApp(ui=ui,server=server)
