#list of packages required
list.of.packages <- c("devtools","shiny","tidyverse","magrittr","stringr","plotly","shinyTypeahead")

#checking missing packages from list
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#install missing ones
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

if ("shinyTypeahead" %in% new.packages) devtools::install_github("ThomasSiegmund/shinyTypeahead")

library(shiny)
library(shinyTypeahead)
library(tidyverse)
library(magrittr)
library(stringr)
library(plotly)


##### Transcriptional landscape #####
EC_TPMs <- readRDS("EC_TPMs.rds")
EC_TPMs_averages <- readRDS("EC_TPMs_averages.rds")
EC_TPMs_tidy <- readRDS("EC_TPMs_tidy.rds")
fraction_names <- c(
  "GFPneg" = "Non-EC",
  "GFPpos" = "EC",
  "Total" = "Total Tissue"
)

p7_brain_intersect_DEG <- readRDS("BrainEC_intersect_DEG.rds")
adult_brain_intersect_DEG <- readRDS("AdultBrainEC_intersect_DEG.rds")
liver_intersect_DEG <- readRDS("LiverEC_intersect_DEG.rds")
lung_intersect_DEG <- readRDS("LungEC_intersect_DEG.rds")
kidney_intersect_DEG <- readRDS("KidneyEC_intersect_DEG.rds")
CulturedEC_DEG <- readRDS("CulturedEC_intersect_DEG.rds")
pan_endothelial <- readRDS("pan_endothelial_genes.rds")
BBB_genes <- readRDS("BBB_genes.rds")
BBB_genes_lost <- readRDS("BBB_genes_lost.rds")


##### Chromatin Accessiblity Landscape #####
P7_Brain_EC_uniqueAccessibleChromatin <- read_tsv("P7_Brain_EC_differentialAccessibleChromatin.bed",col_names = F,col_types = "ccc")
P7_Liver_EC_uniqueAccessibleChromatin <- read_tsv("P7_Liver_EC_differentialAccessibleChromatin.bed",col_names = F,col_types = "ccc")
P7_Lung_EC_uniqueAccessibleChromatin <- read_tsv("P7_Lung_EC_differentialAccessibleChromatin.bed",col_names = F,col_types = "ccc")
P7_Kidney_EC_uniqueAccessibleChromatin <- read_tsv("P7_Kidney_EC_differentialAccessibleChromatin.bed",col_names = F,col_types = "ccc")

ATAC_principalComponents_differentialPeaks <- readRDS("ATAC_principalComponents_differentialPeaks.rds")
highVariance_motifs_differentialPeaks <- readRDS("highVariance_motifs_differentialPeaks.rds")


##### Methylation Landscape #####
P7_Brain_EC_uniqueDMR <- read_tsv("P7_Brain_EC_uniqueDMR.bed",col_names = F,col_types = "ccc")
P7_Liver_EC_uniqueDMR <- read_tsv("P7_Liver_EC_uniqueDMR.bed",col_names = F,col_types = "ccc")
P7_Lung_EC_uniqueDMR <- read_tsv("P7_Lung_EC_uniqueDMR.bed",col_names = F,col_types = "ccc")
P7_Kidney_EC_uniqueDMR <- read_tsv("P7_Kidney_EC_uniqueDMR.bed",col_names = F,col_types = "ccc")

DMR_methyl_principalComponents <- readRDS("DMR_methyl_principalComponents.rds")
DMR_methyl_highVariance_motifs <- readRDS("DMR_methyl_highVariance_motifs.rds")

##### Single Cell ####
endothelial_scaled_summary <- readRDS("endothelial_scaled_summary.rds")

##### Beta-catenin stabilization ####
CVO_TPMs <- readRDS("CVO_TPMs.rds")
CVO_TPMs_tidy <- readRDS("CVO_TPMs_tidy.rds")
CVO_TPMs_averages <- readRDS("CVO_TPMs_averages.rds")
CVO_fraction_names <- c(
  "GFPneg" = "Non-EC",
  "GFPpos" = "EC"
)
cerebellar_BBB_transcripts <- readRDS("cerebellar_BBB_transcripts.rds")
pituitary_EC_transcripts <- readRDS("pituitary_EC_transcripts.rds")
anteriorpituitary_EC_transcripts <- readRDS("anteriorpituitary_EC_transcripts.rds")
posteriorpituitary_EC_transcripts <- readRDS("posteriorpituitary_EC_transcripts.rds")
bcat_anteriorpituitary_EC_transcripts <- readRDS("bcat_anteriorpituitary_EC_transcripts.rds")
bcat_posteriorpituitary_EC_transcripts <- readRDS("bcat_posteriorpituitary_EC_transcripts.rds")

ui <- navbarPage ( "Vascular Endothelial Cell Trans-omics Resource Database",
                   tabPanel("Home",
                            fluidPage(
                              fluidRow(
                                h1("VECTRDB"),
                                p("We present the Vascular Endothelial Cell Trans-omics Resource Database (VECTRDB), a database that integrates and compares the transcriptional, chromatin accessibility, and DNA methylation landscapes of vascular endothelial cells (ECs). We encourage you to explore the genomic landscapes presented here using the tabs at the top of the page."),
                                br(),
                                #p("To run this Shiny app on your local computer, please visit ",a("https://github.com/mfsabbagh/EC_Genomics",href="https://github.com/mfsabbagh/EC_Genomics")),
                                #br(),
                                p(strong("Transcriptional Landscape:")),
                                p("Bulk RNA-seq of postnatal day 7 (P7) vascular ECs isolated from brain, liver, lung, and kidney of Tie2-GFP transgenic mice as well as similarily isolated adult brain ECs and cultured adult brain ECs. In this tab, one can explore organ-specific EC gene expression."),
                                br(),
                                p(strong("Chromatin Accessiblity Landscape:")),
                                p("Bulk ATAC-seq of P7 vascular ECs isolated from brain, liver, lung, and kidney of Tie2-GFP transgenic mice as well as similarily isolated adult brain ECs and cultured adult brain ECs. In this tab, one can explore transcription factor binding motifs that exhibit high variability between organ-specific EC accessible chromatin."),
                                br(),
                                p(strong("CG Methylation Landscape:")),
                                p("Bulk MethylC-seq of P7 vascular ECs isolated from brain, liver, lung, and kidney of Tie2-GFP transgenic mice. In this tab, one can explore transcription factor binding motifs that exhibit high variability between organ-specific EC hypomethylated features."),
                                br(),
                                p(strong("Single-Cell RNA-seq:")),
                                p("We assessed EC gene expression in the developing CNS at single cell resolution by performing single-cell RNA-seq on 3,946 FACS-purified GFP-positive ECs from a P7 Tie2-GFP mouse brain. In this tab, one can explore gene expression in arterial, venous, capillary, mitotic, and tip subtypes of brain ECs."),
                                br(),
                                p(strong("Beta-catenin Stabilization:")),
                                p("Bulk RNA-seq of adult vascular ECs isolated from cerebellum, anterior pituitary, and posterior pituitary of Tie2-GFP transgenic mice with and without beta-catenin stabilization. In this tab, one can explore the effects of beta-catenin stabilization on pituitary EC gene expression."),
                                br(),
                                p("You can also view the P7 vascular EC data on our ",a("AnnoJ Genome Browser.",href = "http://neomorph.salk.edu/Endothelial_cell_methylome.php")),
                                br(),
                                p("This project was the result of a colloboration between the ",a("Nathans Lab,",href="http://nathanslab.mbg.jhmi.edu/"),"the ",a("Ecker Lab,",href="http://ecker.salk.edu/"),"and the ",a("Goff Lab",href="http://www.gofflab.org/")),
                                br(),
                                h3("Please feel free to reproduce images generated through this online visualization tool. We only ask that you provide proper citations as follows:"),
                                h4("For P7 RNA, ATAC, Methyl, and Single Cell analyses:"),
                                h5("Mark F. Sabbagh, Jacob S. Heng, Chongyuan Luo, Rosa G. Castanon, Joseph R. Nery, Amir Rattner, Loyal A. Goff, Joseph R. Ecker, and Jeremy Nathans. \"Transcriptional and epigenomic landscapes of CNS and non-CNS vascular endothelial cells.\" eLife 2018;7:e36187 DOI: ",a("10.7554/eLife.36187",href="https://doi.org/10.7554/eLife.36187")),
                                h4("For Beta-catenin Stabilization analyses:"),
                                h5("Yanshu Wang, Mark F. Sabbagh, Xiaowu Gu, Amir Rattner, John Williams, Jeremy Nathans. \"Beta-catenin signaling regulates barrier-specific gene expression in circumventricular organ and ocular vasculatures.\" eLife 2019;8:e43257 DOI: ",a("10.7554/eLife.43257",href="https://doi.org/10.7554/eLife.43257"))
                              )
                            )),
                   tabPanel("Transcriptional Landscape",
                            fluidPage(
                              h4("Please cite as:"),
                              h4("Mark F. Sabbagh, Jacob S. Heng, Chongyuan Luo, Rosa G. Castanon, Joseph R. Nery, Amir Rattner, Loyal A. Goff, Joseph R. Ecker, and Jeremy Nathans. \"Transcriptional and epigenomic landscapes of CNS and non-CNS vascular endothelial cells.\" eLife 2018;7:e36187 DOI: ",a("10.7554/eLife.36187",href="https://doi.org/10.7554/eLife.36187")),
                              sidebarLayout(position="left",
                                            sidebarPanel(
                                              typeaheadInput("MS_gene_name",h3("Enter gene name and click submit to view endothelial expression levels"),value="Pecam1",choices = EC_TPMs$Gene),
                                              actionButton("MS_submit_genes",label=h5("Submit")),
                                              br(),
                                              radioButtons("ECTSG",label=h3("Display TPM expression values from which category of genes:"),choices=c("Pan-endothelial","Enriched only in P7 liver ECs","Enriched only in P7 lung ECs","Enriched only in P7 kidney ECs","Enriched only in P7 brain ECs","Enriched only in adult brain ECs","Enriched only in cultured ECs","Blood-brain barrier (BBB) genes","BBB genes lost upon culture")),
                                              downloadButton("downloadData", "Download expression data for selected option")
                                            ),
                                            mainPanel(
                                              h3("Expression of chosen gene in endothelial cells, non-endothelial cells, and total tissue:"),
                                              plotOutput("MS_expression_plot"),
                                              tableOutput("chosen_gene"),
                                              br(),
                                              h3("Table of TPMs for gene list chosen from side panel:"),
                                              dataTableOutput("differential_genes")
                                            )
                              )
                              
                            )
                   ),
                   tabPanel("Chromatin Accessiblity Landscape",
                            fluidPage(
                              h4("Please cite as:"),
                              h4("Mark F. Sabbagh, Jacob S. Heng, Chongyuan Luo, Rosa G. Castanon, Joseph R. Nery, Amir Rattner, Loyal A. Goff, Joseph R. Ecker, and Jeremy Nathans. \"Transcriptional and epigenomic landscapes of CNS and non-CNS vascular endothelial cells.\" eLife 2018;7:e36187 DOI: ",a("10.7554/eLife.36187",href="https://doi.org/10.7554/eLife.36187")),
                              sidebarLayout(position="left",
                                            sidebarPanel(
                                              selectInput("unique_chromatin", "Choose list of accessible chromatin unique to which endothelium:",
                                                          choices = c("P7 Brain", "P7 Liver", "P7 Lung","P7 Kidney")),
                                              
                                              downloadButton("downloadDataChromatin", "Download bed file")
                                              
                                            ),
                                            mainPanel(
                                              h2("Explore high-variance transcription factor motifs across differentially accessible chromatin"),
                                              br(),
                                              selectizeInput("differentialATAC_tf_motif",label="Choose motif to display accessibility deviation Z-score",choices = highVariance_motifs_differentialPeaks$name,selected = "Zic(Zf)"),
                                              plotlyOutput("differentialATAC_PCA_motif"),
                                              dataTableOutput("table_chromatin")
                                              
                                            )
                              )
                            )
                   ),
                   tabPanel("CG Methylation Landscape",
                            fluidPage(
                              h4("Please cite as:"),
                              h4("Mark F. Sabbagh, Jacob S. Heng, Chongyuan Luo, Rosa G. Castanon, Joseph R. Nery, Amir Rattner, Loyal A. Goff, Joseph R. Ecker, and Jeremy Nathans. \"Transcriptional and epigenomic landscapes of CNS and non-CNS vascular endothelial cells.\" eLife 2018;7:e36187 DOI: ",a("10.7554/eLife.36187",href="https://doi.org/10.7554/eLife.36187")),
                              fluidRow(
                                sidebarLayout(position="left",
                                              sidebarPanel(
                                                selectInput("unique_DMR", "Choose list of differentially methylated regions unique to which endothelium:",
                                                            choices = c("P7 Brain", "P7 Liver", "P7 Lung","P7 Kidney")),
                                                
                                                downloadButton("downloadDataMethyl", "Download bed file")
                                                
                                              ),
                                              mainPanel(
                                                h2("Explore high-variance transcription factor motifs across differentially hypomethylated regions (DMRs)"),
                                                br(),
                                                selectizeInput("DMR_Methyl_tf_motif",label="Choose motif to display accessibility deviation Z-score",choices = DMR_methyl_highVariance_motifs$name,selected = "Zic(Zf)"),
                                                plotlyOutput("DMR_Methyl_PCA_motif"),
                                                dataTableOutput("table_DMR")
                                                
                                              )
                                )
                              )
                            )
                   ),
                   tabPanel("Single-Cell RNA-seq",
                            fluidPage(
                              h4("Please cite as:"),
                              h4("Mark F. Sabbagh, Jacob S. Heng, Chongyuan Luo, Rosa G. Castanon, Joseph R. Nery, Amir Rattner, Loyal A. Goff, Joseph R. Ecker, and Jeremy Nathans. \"Transcriptional and epigenomic landscapes of CNS and non-CNS vascular endothelial cells.\" eLife 2018;7:e36187 DOI: ",a("10.7554/eLife.36187",href="https://doi.org/10.7554/eLife.36187")),
                              sidebarLayout(position="left",
                                            sidebarPanel(
                                              h4("Explore the developing brain EC single-cell expression dataset"),
                                              br(),
                                              typeaheadInput("sc_gene_name","Enter gene name",value="Bmx",choices = unique(endothelial_scaled_summary$gene_short_name)),
                                              actionButton("sc_submit_genes",label=h5("Submit"))
                                              
                                            ),
                                            mainPanel(
                                              plotOutput("sc_plot")
                                            )
                                            
                              )
                            )
                   ),
                   tabPanel("Beta-catenin Stabilization",
                            fluidPage(
                              h4("Please cite as:"),
                              h4("Yanshu Wang, Mark F. Sabbagh, Xiaowu Gu, Amir Rattner, John Williams, Jeremy Nathans. \"Beta-catenin signaling regulates barrier-specific gene expression in circumventricular organ and ocular vasculatures.\" eLife 2019;8:e43257 DOI: ",a("10.7554/eLife.43257",href="https://doi.org/10.7554/eLife.43257")),
                              sidebarLayout(position="left",
                                            sidebarPanel(
                                              typeaheadInput("CVO_gene_name",h3("Enter gene name and click submit to view endothelial expression levels"),value="Pecam1",choices = CVO_TPMs$Gene),
                                              actionButton("CVO_submit_genes",label=h5("Submit")),
                                              br(),
                                              radioButtons("CVO_ECTSG",label=h3("Display TPM expression values from which category of genes:"),choices=c("Cerebellar blood-brain barrier","Pituitary endothelial cells","Enriched only in anterior pituitary ECs","Enriched only in posterior pituitary ECs","Beta-catenin stabilization responsive in anterior pituitary ECs","Beta-catenin stabilization responsive in posterior pituitary ECs")),
                                              downloadButton("downloadDataCVO", "Download expression data for selected option")
                                            ),
                                            mainPanel(
                                              h3("Expression of chosen gene in endothelial cells and non-endothelial cells with and without beta-catenin stabilization:"),
                                              plotOutput("CVO_expression_plot"),
                                              tableOutput("CVO_chosen_gene"),
                                              br(),
                                              h3("Table of TPMs for gene list chosen from side panel:"),
                                              dataTableOutput("CVO_differential_genes")
                                            )
                              )
                              
                            )
                   )
                   
)

server <- function(input, output, session) {
  MS_genes_to_plot <- eventReactive(input$MS_submit_genes,input$MS_gene_name,ignoreNULL = F)
  CVO_genes_to_plot <- eventReactive(input$CVO_submit_genes,input$CVO_gene_name,ignoreNULL = F)
  
  ECTSG_toShow <- reactive (input$ECTSG)
  CVO_ECTSG_toShow <- reactive(input$CVO_ECTSG)
  
  differentialATAC_motifs_to_plot <- reactive(input$differentialATAC_tf_motif)
  DMR_Methyl_motifs_to_plot <- reactive(input$DMR_Methyl_tf_motif)
  
  
  
  # Reactive value for selected dataset ----
  datasetInput <- reactive({
    if (ECTSG_toShow()=="Pan-endothelial") {
      semi_join(EC_TPMs_averages,pan_endothelial,by=c("Gene"))
    } else if (ECTSG_toShow() =="Enriched only in P7 brain ECs") {
      semi_join(EC_TPMs_averages,p7_brain_intersect_DEG,by=c("Gene")) %>% 
        mutate("Fold Change"=`P7_BrainEC_TPMs`/(((`P7_LiverEC_TPMs` + `P7_LungEC_TPMs` + `P7_KidneyEC_TPMs`)/3)+1)) %>%
        arrange(desc(`Fold Change`)) 
    } else if (ECTSG_toShow() =="Enriched only in P7 liver ECs") {
      semi_join(EC_TPMs_averages,liver_intersect_DEG,by=c("Gene")) %>% 
        mutate("Fold Change"=`P7_LiverEC_TPMs`/(((`P7_BrainEC_TPMs` + `P7_LungEC_TPMs` + `P7_KidneyEC_TPMs`)/3)+1)) %>% 
        arrange(desc(`Fold Change`)) 
    } else if (ECTSG_toShow() =="Enriched only in P7 lung ECs") {
      semi_join(EC_TPMs_averages,lung_intersect_DEG,by=c("Gene")) %>% 
        mutate("Fold Change"=`P7_LungEC_TPMs`/(((`P7_LiverEC_TPMs` + `P7_BrainEC_TPMs` + `P7_KidneyEC_TPMs`)/3)+1)) %>% 
        arrange(desc(`Fold Change`)) 
    } else if (ECTSG_toShow() =="Enriched only in P7 kidney ECs") {
      semi_join(EC_TPMs_averages,kidney_intersect_DEG,by=c("Gene")) %>% 
        mutate("Fold Change"=`P7_KidneyEC_TPMs`/(((`P7_LiverEC_TPMs` + `P7_LungEC_TPMs` + `P7_BrainEC_TPMs`)/3)+1)) %>% 
        arrange(desc(`Fold Change`)) 
    } else if (ECTSG_toShow() =="Enriched only in cultured ECs") {
      semi_join(EC_TPMs_averages,CulturedEC_DEG,by=c("Gene")) %>% 
        mutate("Fold Change"=`Cultured_BrainEC_TPMs`/(`P7_BrainEC_TPMs` + 1)) %>% 
        arrange(desc(`Fold Change`)) 
    } else if (ECTSG_toShow() =="Blood-brain barrier (BBB) genes") {
      semi_join(EC_TPMs_averages,BBB_genes,by=c("Gene")) %>% 
        mutate("Fold Change"=`P7_BrainEC_TPMs`/(((`P7_LiverEC_TPMs` + `P7_LungEC_TPMs` + `P7_KidneyEC_TPMs`)/3)+1)) %>% 
        arrange(desc(`Fold Change`)) 
    } else if (ECTSG_toShow() =="BBB genes lost upon culture") {
      semi_join(EC_TPMs_averages,BBB_genes_lost,by=c("Gene")) %>% 
        mutate("Fold Change"=`P7_BrainEC_TPMs`/(((`P7_LiverEC_TPMs` + `P7_LungEC_TPMs` + `P7_KidneyEC_TPMs`)/3)+1)) %>% 
        arrange(desc(`Fold Change`)) 
    } else if (ECTSG_toShow() =="Enriched only in adult brain ECs") {
      semi_join(EC_TPMs_averages,adult_brain_intersect_DEG,by=c("Gene")) %>% 
        mutate("Fold Change"=`Adult_BrainEC_TPMs`/(((`P7_LiverEC_TPMs` + `P7_LungEC_TPMs` + `P7_KidneyEC_TPMs`)/3)+1)) %>% 
        arrange(desc(`Fold Change`)) 
    } 
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$ECTSG, ".tsv", sep = "")
    },
    content = function(file) {
      write_tsv(datasetInput(), file)
    }
  )
  
  datasetInputCVO <- reactive({
    if (CVO_ECTSG_toShow()=="Cerebellar blood-brain barrier") {
      cerebellar_BBB_transcripts %>% 
        mutate("Fold Change"=`WT_Cerebellum`/(((`WT_AnteriorPituitary` + `WT_PosteriorPituitary`)/2)+1)) %>%
        arrange(desc(`Fold Change`)) 
    } else if (CVO_ECTSG_toShow() =="Pituitary endothelial cells") {
      pituitary_EC_transcripts %>% 
        mutate("Fold Change"=((`WT_AnteriorPituitary` + `WT_PosteriorPituitary`)/2)/(`WT_Cerebellum`+1)) %>%
        arrange(desc(`Fold Change`)) 
    } else if (CVO_ECTSG_toShow() =="Enriched only in anterior pituitary ECs") {
      anteriorpituitary_EC_transcripts %>% 
        mutate("Fold Change"=`WT_AnteriorPituitary`/(((`WT_Cerebellum` + `WT_PosteriorPituitary`)/2)+1)) %>%
        arrange(desc(`Fold Change`)) 
    } else if (CVO_ECTSG_toShow() =="Enriched only in posterior pituitary ECs") {
      posteriorpituitary_EC_transcripts %>% 
        mutate("Fold Change"=`WT_PosteriorPituitary`/(((`WT_Cerebellum` + `WT_AnteriorPituitary`)/2)+1)) %>%
        arrange(desc(`Fold Change`))  
    } else if (CVO_ECTSG_toShow() =="Beta-catenin stabilization responsive in anterior pituitary ECs") {
      bcat_anteriorpituitary_EC_transcripts %>% 
        mutate("Fold Change"=`Mutant_AnteriorPituitary`/(`WT_AnteriorPituitary`+1)) %>%
        arrange(desc(`Fold Change`)) 
    } else if (CVO_ECTSG_toShow() =="Beta-catenin stabilization responsive in posterior pituitary ECs") {
      bcat_posteriorpituitary_EC_transcripts %>% 
        mutate("Fold Change"=`Mutant_PosteriorPituitary`/(`WT_PosteriorPituitary`+1)) %>%
        arrange(desc(`Fold Change`)) 
    }
  })
  
  output$downloadDataCVO <- downloadHandler(
    filename = function() {
      paste(input$CVO_ECTSG, ".tsv", sep = "")
    },
    content = function(file) {
      write_tsv(datasetInputCVO(), file)
    }
  )
  
  
  output$chosen_gene <- renderTable(
    {
      dplyr::filter(EC_TPMs_averages,Gene == MS_genes_to_plot())
    }
  )
  
  output$CVO_chosen_gene <- renderTable(
    {
      dplyr::filter(CVO_TPMs_averages,Gene == CVO_genes_to_plot())
    }
  )
  
  output$differential_genes <- renderDataTable(
    {
      if (ECTSG_toShow()=="Pan-endothelial") {
        semi_join(EC_TPMs_averages,pan_endothelial,by=c("Gene"))
      } else if (ECTSG_toShow() =="Enriched only in P7 brain ECs") {
        semi_join(EC_TPMs_averages,p7_brain_intersect_DEG,by=c("Gene")) %>% 
          mutate("Fold Change"=`P7_BrainEC_TPMs`/(((`P7_LiverEC_TPMs` + `P7_LungEC_TPMs` + `P7_KidneyEC_TPMs`)/3)+1)) %>%
          arrange(desc(`Fold Change`)) 
      } else if (ECTSG_toShow() =="Enriched only in P7 liver ECs") {
        semi_join(EC_TPMs_averages,liver_intersect_DEG,by=c("Gene")) %>% 
          mutate("Fold Change"=`P7_LiverEC_TPMs`/(((`P7_BrainEC_TPMs` + `P7_LungEC_TPMs` + `P7_KidneyEC_TPMs`)/3)+1)) %>% 
          arrange(desc(`Fold Change`)) 
      } else if (ECTSG_toShow() =="Enriched only in P7 lung ECs") {
        semi_join(EC_TPMs_averages,lung_intersect_DEG,by=c("Gene")) %>% 
          mutate("Fold Change"=`P7_LungEC_TPMs`/(((`P7_LiverEC_TPMs` + `P7_BrainEC_TPMs` + `P7_KidneyEC_TPMs`)/3)+1)) %>% 
          arrange(desc(`Fold Change`)) 
      } else if (ECTSG_toShow() =="Enriched only in P7 kidney ECs") {
        semi_join(EC_TPMs_averages,kidney_intersect_DEG,by=c("Gene")) %>% 
          mutate("Fold Change"=`P7_KidneyEC_TPMs`/(((`P7_LiverEC_TPMs` + `P7_LungEC_TPMs` + `P7_BrainEC_TPMs`)/3)+1)) %>% 
          arrange(desc(`Fold Change`)) 
      } else if (ECTSG_toShow() =="Enriched only in cultured ECs") {
        semi_join(EC_TPMs_averages,CulturedEC_DEG,by=c("Gene")) %>% 
          mutate("Fold Change"=`Cultured_BrainEC_TPMs`/(`P7_BrainEC_TPMs` + 1)) %>% 
          arrange(desc(`Fold Change`)) 
      } else if (ECTSG_toShow() =="Blood-brain barrier (BBB) genes") {
        semi_join(EC_TPMs_averages,BBB_genes,by=c("Gene")) %>% 
          mutate("Fold Change"=`P7_BrainEC_TPMs`/(((`P7_LiverEC_TPMs` + `P7_LungEC_TPMs` + `P7_KidneyEC_TPMs`)/3)+1)) %>% 
          arrange(desc(`Fold Change`)) 
      } else if (ECTSG_toShow() =="BBB genes lost upon culture") {
        semi_join(EC_TPMs_averages,BBB_genes_lost,by=c("Gene")) %>% 
          mutate("Fold Change"=`P7_BrainEC_TPMs`/(((`P7_LiverEC_TPMs` + `P7_LungEC_TPMs` + `P7_KidneyEC_TPMs`)/3)+1)) %>% 
          arrange(desc(`Fold Change`)) 
      } else if (ECTSG_toShow() =="Enriched only in adult brain ECs") {
        semi_join(EC_TPMs_averages,adult_brain_intersect_DEG,by=c("Gene")) %>% 
          mutate("Fold Change"=`Adult_BrainEC_TPMs`/(((`P7_LiverEC_TPMs` + `P7_LungEC_TPMs` + `P7_KidneyEC_TPMs`)/3)+1)) %>% 
          arrange(desc(`Fold Change`)) 
      } 
      
      
    },
    options=list(pageLength=5)
  )
  
  output$CVO_differential_genes <- renderDataTable(
    {
      if (CVO_ECTSG_toShow()=="Cerebellar blood-brain barrier") {
        cerebellar_BBB_transcripts %>% 
          mutate("Fold Change"=`WT_Cerebellum`/(((`WT_AnteriorPituitary` + `WT_PosteriorPituitary`)/2)+1)) %>%
          arrange(desc(`Fold Change`)) 
      } else if (CVO_ECTSG_toShow() =="Pituitary endothelial cells") {
        pituitary_EC_transcripts %>% 
          mutate("Fold Change"=((`WT_AnteriorPituitary` + `WT_PosteriorPituitary`)/2)/(`WT_Cerebellum`+1)) %>%
          arrange(desc(`Fold Change`)) 
      } else if (CVO_ECTSG_toShow() =="Enriched only in anterior pituitary ECs") {
        anteriorpituitary_EC_transcripts %>% 
          mutate("Fold Change"=`WT_AnteriorPituitary`/(((`WT_Cerebellum` + `WT_PosteriorPituitary`)/2)+1)) %>%
          arrange(desc(`Fold Change`)) 
      } else if (CVO_ECTSG_toShow() =="Enriched only in posterior pituitary ECs") {
        posteriorpituitary_EC_transcripts %>% 
          mutate("Fold Change"=`WT_PosteriorPituitary`/(((`WT_Cerebellum` + `WT_AnteriorPituitary`)/2)+1)) %>%
          arrange(desc(`Fold Change`))  
      } else if (CVO_ECTSG_toShow() =="Beta-catenin stabilization responsive in anterior pituitary ECs") {
        bcat_anteriorpituitary_EC_transcripts %>% 
          mutate("Fold Change"=`Mutant_AnteriorPituitary`/(`WT_AnteriorPituitary`+1)) %>%
          arrange(desc(`Fold Change`)) 
      } else if (CVO_ECTSG_toShow() =="Beta-catenin stabilization responsive in posterior pituitary ECs") {
        bcat_posteriorpituitary_EC_transcripts %>% 
          mutate("Fold Change"=`Mutant_PosteriorPituitary`/(`WT_PosteriorPituitary`+1)) %>%
          arrange(desc(`Fold Change`)) 
      }
      
    },
    options=list(pageLength=5)
  )
  
  
  
  output$MS_expression_plot <- renderPlot(
    {
      ggplot(filter(EC_TPMs_tidy,Gene %in% MS_genes_to_plot()),aes(x=Tissue,y=Expression,color=Tissue,shape=Fraction,alpha=Replicate)) +
        geom_point(size=4,position = position_dodge(width=0.9))+ 
        theme_gray(base_size = 20) + ylab("Expression (TPM)")+ 
        scale_color_manual(values=c("Brain"="blue","Liver"="Cyan","Lung"="magenta","Kidney"="green3","Cultured"="orange","AdultBrain"="darkblue"),guide=F)+ 
        expand_limits(y=0) + scale_alpha_manual(values=c("R1"=1,"R2"=1),guide=F) + scale_shape_discrete(name="Cell or tissue fraction",labels=c("Non-EC","EC","Total tissue"))+ 
        facet_grid(Gene ~ Fraction,scales="free_x", labeller = labeller(Fraction = as_labeller(fraction_names))) + theme(axis.text.x = element_text(angle = 90),strip.background = element_blank(), strip.text = element_text(face = "italic"))
      
    })
  
  output$CVO_expression_plot <- renderPlot(
    {
      ggplot(filter(CVO_TPMs_tidy,Gene %in% CVO_genes_to_plot()),aes(x=Tissue,y=log2(Expression+1),color=Tissue,shape=Fraction,alpha=Replicate)) +
        geom_point(size=4,position = position_dodge(width=0.9))+ 
        theme_gray(base_size = 20) + ylab("Expression (log2(TPM+1))")+ 
        scale_color_manual(values=c("Cerebellum"="blue","AnteriorPituitary"="yellow3","PosteriorPituitary"="firebrick1"),guide=F)+ 
        expand_limits(y=0) + scale_alpha_manual(values=c("R1"=1,"R2"=1),guide=F) + scale_shape_discrete(name="Cell fraction",labels=c("Non-EC","EC"))+ 
        facet_grid(Gene ~ Genotype,scales="free_x", labeller = labeller(Fraction = as_labeller(CVO_fraction_names))) + theme(axis.text.x = element_text(angle = 90),strip.background = element_blank(), strip.text = element_text(face = "italic"))
      
    })
  
  # Reactive value for selected dataset ----
  datasetInputChromatin <- reactive({
    switch(input$unique_chromatin,
           "P7 Brain" = P7_Brain_EC_uniqueAccessibleChromatin,
           "P7 Liver" = P7_Liver_EC_uniqueAccessibleChromatin,
           "P7 Lung" = P7_Lung_EC_uniqueAccessibleChromatin,
           "P7 Kidney" = P7_Kidney_EC_uniqueAccessibleChromatin)
  })
  
  # Table of selected dataset ----
  output$table_chromatin <- renderDataTable({
    datasetInputChromatin() %>% dplyr::rename(Chromosome=X1,Start=X2,End=X3)
  })
  
  # Downloadable csv of selected dataset ----
  output$downloadDataChromatin <- downloadHandler(
    filename = function() {
      paste(input$unique_chromatin," endothelial cell unique accessible chromatin.bed", sep = "")
    },
    content = function(file) {
      write_tsv(datasetInputChromatin(), file,col_names = F)
    }
  )
  
  
  # Reactive value for selected dataset ----
  datasetInputDMR <- reactive({
    switch(input$unique_DMR,
           "P7 Brain" = P7_Brain_EC_uniqueDMR,
           "P7 Liver" = P7_Liver_EC_uniqueDMR,
           "P7 Lung" = P7_Lung_EC_uniqueDMR,
           "P7 Kidney" = P7_Kidney_EC_uniqueDMR)
  })
  
  # Table of selected dataset ----
  output$table_DMR <- renderDataTable({
    datasetInputDMR() %>% dplyr::rename(Chromosome=X1,Start=X2,End=X3)
  })
  
  # Downloadable csv of selected dataset ----
  output$downloadDataDMR <- downloadHandler(
    filename = function() {
      paste(input$unique_DMR," endothelial cell unique differentially methylated region.bed", sep = "")
    },
    content = function(file) {
      write_tsv(datasetInputDMR(), file,col_names = F)
    }
  )
  
  output$differentialATAC_PCA_motif <- renderPlotly( {
    
    plot_ly(ATAC_principalComponents_differentialPeaks,hoverinfo="text",text=~tissue, x = ~PC1, y = ~PC2, z = ~PC3, marker = list(color = ~get(differentialATAC_motifs_to_plot()), colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
      add_markers() %>% add_text() %>%
      layout(scene = list(xaxis = list(title = 'PC1 (30% variance)'),
                          yaxis = list(title = 'PC2 (25% variance)'),
                          zaxis = list(title = 'PC3 (18% variance)')),
             annotations = list(
               x = 1.05,
               y = 1.05,
               text = 'Z-score',
               xref = 'paper',
               yref = 'paper',
               showarrow = FALSE
             ))
    
    
  }
  )
  
  output$DMR_Methyl_PCA_motif <- renderPlotly( {
    
    plot_ly(DMR_methyl_principalComponents,hoverinfo="text",text=~tissue, x = ~PC1, y = ~PC2, z = ~PC3, marker = list(color = ~get(DMR_Methyl_motifs_to_plot()), colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
      add_markers() %>% add_text() %>%
      layout(scene = list(xaxis = list(title = 'PC1 (40% variance)'),
                          yaxis = list(title = 'PC2 (31% variance)'),
                          zaxis = list(title = 'PC3 (25% variance)')),
             annotations = list(
               x = 1.05,
               y = 1.05,
               text = 'Z-score',
               xref = 'paper',
               yref = 'paper',
               showarrow = FALSE
             ))
    
    
  }
  )
  
  
  sc_genes_to_plot <- eventReactive (input$sc_submit_genes, input$sc_gene_name, ignoreNULL = F)
  
  output$sc_plot <- renderPlot ( {
    ggplot(filter(endothelial_scaled_summary,gene_short_name %in% sc_genes_to_plot()),aes(x=CellType,y=mean,fill=CellType)) +
      geom_bar(stat="identity") + 
      geom_errorbar(aes(x=CellType,ymin=mean-std_error,ymax=mean+std_error),size=0.9,width=0.5) +
      theme_classic(base_size = 20) + ylab("Mean UMIs per 10000") + ggtitle(sc_genes_to_plot())
    
    
  } )
  
  
  
}


# Run the application 
shinyApp(ui = ui, server = server)

