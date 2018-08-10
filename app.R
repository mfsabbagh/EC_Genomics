#list of packages required
list.of.packages <- c("devtools","shiny","tidyverse","magrittr","reshape2","stringr","plotly","shinyTypeahead","monocle","cellrangerRkit")

#checking missing packages from list
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#install missing ones
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

if ("shinyTypeahead" %in% new.packages) devtools::install_github("ThomasSiegmund/shinyTypeahead")

if ("monocle" %in% new.packages) {
  source("https://bioconductor.org/biocLite.R") 
  biocLite("monocle")
}

#devtools::install_github("buenrostrolab/cellrangerRkit")

if("cellrangerRkit" %in% new.packages) source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")

library(shiny)
library(shinyTypeahead)
library(tidyverse)
library(magrittr)
library(stringr)
library(plotly)
library(monocle)
library(cellrangerRkit)
library(reshape2)

visualize_gene_markers2 <-function (gbm, gene_probes, projection, limits = c(0, 0.5), low_col = "lightblue", high_col = "darkblue", marker_size = 0.1, 
                                    title = NULL, axis_line_size = 1.0, axis_tick_size = 1.0, panel_border_size = 1.0) 
{
  gbm_trunc <- trunc_gbm_by_genes(gbm, gene_probes)
  gene_values <- t(as.matrix(Biobase::exprs(gbm_trunc)))
  gene_values[gene_values < limits[1]] <- limits[1]
  gene_values[gene_values > limits[2]] <- limits[2]
  colnames(gene_values) <- gene_probes
  projection_names <- colnames(projection)
  colnames(projection) <- c("Component.1", "Component.2")
  proj_gene <- data.frame(cbind(projection, gene_values))
  proj_gene_melt <- melt(proj_gene, id.vars = c("Component.1", "Component.2"))
  colnames(proj_gene_melt) <- c("Component.1", "Component.2", "gene_probe", "value")
  p <- ggplot(proj_gene_melt, aes(Component.1, Component.2)) + 
    geom_point(aes(colour = value), size = marker_size) + 
    facet_wrap(~gene_probe,ncol=7, scales = "free") + scale_colour_gradient(low = low_col, 
                                                                            high = high_col, name = "value") + theme(strip.background = element_rect(colour = "black", fill = "white"))
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  p <- p + theme_bw() + theme(panel.border = element_rect(linetype = "solid", size = panel_border_size, fill = NA),
                              plot.title = element_text(hjust = 0.5),
                              axis.line.x = element_line(colour = 'black', size = axis_line_size),
                              axis.line.y = element_line(colour = 'black', size = axis_line_size),
                              axis.ticks = element_line(colour = "black", size = axis_tick_size),
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                              strip.background = element_rect(colour = "black", fill = "white"),
                              legend.key.height = unit(0.85, "cm"))
  
  return(p)
}

##### Transcriptional landscape #####
RNA_principalComponents <- readRDS("RNA_principalComponents.rds")
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
ATAC_principalComponents <- readRDS("ATAC_principalComponents.rds")
highVariance_motifs <- readRDS("highVariance_motifs.rds")

ATAC_principalComponents_differentialPeaks <- readRDS("ATAC_principalComponents_differentialPeaks.rds")
highVariance_motifs_differentialPeaks <- readRDS("highVariance_motifs_differentialPeaks.rds")

##### Methylation Landscape #####
methyl_principalComponents <- readRDS("methyl_principalComponents.rds")
methyl_highVariance_motifs <- readRDS("methyl_highVariance_motifs.rds")

DMR_methyl_principalComponents <- readRDS("DMR_methyl_principalComponents.rds")
DMR_methyl_highVariance_motifs <- readRDS("DMR_methyl_highVariance_motifs.rds")

##### Single Cell ####
gbm_log <- readRDS("gbm_log.rds")
tsne_proj3 <- readRDS("tsne_proj3.rds")
gene_list <- readRDS("gene_list.rds")

ui <- navbarPage ( "Vascular Endothelial Cell Trans-omics Resource Database",
                   tabPanel("Home",
                            fluidPage(
                              fluidRow(
                                h1("VECTRDB"),
                                p("We present the Vascular Endothelial Cell Trans-omics Resource Database (VECTRDB), a database that integrates and compares the transcriptional, chromatin accessibility, and DNA methylation landscapes of vascular endothelial cells (ECs). We encourage you to explore the genomic landscapes presented here using the tabs at the top of the page."),
                                br(),
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
                                p("You can also view the data on our ",a("AnnoJ Genome Browser.",href = "http://neomorph.salk.edu/Endothelial_cell_methylome.php")),
                                br(),
                                #p("To explore our developing brain endothelial cell single-cell RNA-seq database, please click ",a("here.",href="https://markfsabbagh.shinyapps.io/developingbrainecscrnaseq/")),
                                #br(),
                                p("This project was the result of a colloboration between the ",a("Nathans Lab,",href="http://nathanslab.mbg.jhmi.edu/"),"the ",a("Ecker Lab,",href="http://ecker.salk.edu/"),"and the ",a("Goff Lab",href="http://www.gofflab.org/")),
                                br(),
                                h4("Please cite:"),
                                h4("Mark F. Sabbagh, Jacob S. Heng, Chongyuan Luo, Rosa G. Castanon, Joseph R. Nery, Amir Rattner, Loyal A. Goff, Joseph R. Ecker, and Jeremy Nathans. Vascular Endothelial Cell Trans-omics Resource Database (VECTRDB). https://markfsabbagh.shinyapps.io/vectrdb/")
                              )
                            )),
                   tabPanel("Transcriptional Landscape",
                            fluidPage(
                              sidebarLayout(position="left",
                                            sidebarPanel(
                                              radioButtons("ECTSG",label=h3("Display TPM expression values from which category of genes:"),choices=c("Pan-endothelial","Enriched only in P7 liver ECs","Enriched only in P7 lung ECs","Enriched only in P7 kidney ECs","Enriched only in P7 brain ECs","Enriched only in adult brain ECs","Enriched only in cultured ECs","Blood-brain barrier (BBB) genes","BBB genes lost upon culture")),
                                              br(),
                                              typeaheadInput("MS_gene_name",h3("Enter gene name and click submit to view endothelial expression levels"),value="Pecam1",choices = EC_TPMs$Gene),
                                              actionButton("MS_submit_genes",label=h5("Submit"))
                                            ),
                                            mainPanel(
                                              h3("Table of TPMs for gene list chosen from side panel:"),
                                              dataTableOutput("differential_genes"),
                                              br(),
                                              h3("Expression of chosen gene in endothelial cells, non-endothelial cells, and total tissue:"),
                                              plotOutput("MS_expression_plot"),
                                              tableOutput("chosen_gene"),
                                              br(),
                                              h3("Principal component analysis using gene expression of EC-enriched genes with log2( TPM + 1) heatmap overlaid for chosen gene:"),
                                              fluidRow(
                                                column(6,
                                                       plot_ly(RNA_principalComponents, x= ~PC1, y= ~PC2, z = ~PC3, color= ~tissue,colors = c("Brain"="blue","Liver"="Cyan","Lung"="magenta","Kidney"="green3","Cultured"="orange","AdultBrain"="darkblue")) %>% 
                                                         add_markers() %>% 
                                                         layout(title = "PCA using EC-enriched gene expression", scene = list(xaxis = list(title="PC1 (26% variance)"),
                                                                                                                              yaxis = list(title="PC2 (20% variance)"),
                                                                                                                              zaxis = list(title="PC3 (15% variance)")))
                                                ),
                                                column(5, offset = 1,plotlyOutput("RNA_heatmap")
                                                )
                                              )
                                            )
                              )
                              
                            )
                   ),
                   tabPanel("Chromatin Accessiblity Landscape",
                            fluidPage(
                              fluidRow(
                                column(6,
                                       plot_ly(ATAC_principalComponents,hoverinfo="text",text=~age, x = ~PC1, y = ~PC2, z = ~PC3, color = ~tissue,colors = c("Brain"="blue","Liver"="Cyan","Lung"="magenta","Kidney"="green3","Cultured"="orange","AdultBrain"="darkblue")) %>%
                                         add_markers() %>%
                                         layout(title = "PCA using accessibility at motif-containing peaks", scene = list(xaxis = list(title = 'PC1 (30% variance)'),
                                                                                                                          yaxis = list(title = 'PC2 (25% variance)'),
                                                                                                                          zaxis = list(title = 'PC3 (18% variance)'))),
                                       br(),
                                       br(),
                                       plot_ly(ATAC_principalComponents_differentialPeaks,hoverinfo = "text",text=~age, x = ~PC1, y = ~PC2, z = ~PC3, color = ~tissue,colors=c("Brain"="blue","Liver"="Cyan","Lung"="magenta","Kidney"="green3","Cultured"="orange","AdultBrain"="darkblue")) %>%
                                         add_markers() %>%
                                         layout(title = "PCA using accessibility at motif-containing differential peaks",scene = list(xaxis = list(title = 'PC1 (30% variance)'),
                                                                                                                                      yaxis = list(title = 'PC2 (25% variance)'),
                                                                                                                                      zaxis = list(title = 'PC3 (18% variance)')))
                                ),
                                column(6,
                                       h2("Explore high-variance transcription factor motifs across accessible chromatin"),
                                       br(),
                                       selectizeInput("ATAC_tf_motif",label="Choose motif to display accessibility deviation Z-score",choices = highVariance_motifs$name,selected = "Zic(Zf)"),
                                       plotlyOutput("ATAC_PCA_motif"),
                                       br(),
                                       br(),
                                       h2("Explore high-variance transcription factor motifs across differentially accessible chromatin"),
                                       br(),
                                       selectizeInput("differentialATAC_tf_motif",label="Choose motif to display accessibility deviation Z-score",choices = highVariance_motifs_differentialPeaks$name,selected = "Zic(Zf)"),
                                       plotlyOutput("differentialATAC_PCA_motif")
                                )
                              )
                            )
                   ),
                   tabPanel("CG Methylation Landscape",
                            fluidPage(
                              fluidRow(
                                column(6,
                                       plot_ly(methyl_principalComponents, x = ~PC1, y = ~PC2, z = ~PC3, color = ~tissue,colors=c("Brain"="blue","Liver"="Cyan","Lung"="magenta","Kidney"="green3","Cultured"="orange","AdultBrain"="darkblue")) %>%
                                         add_markers() %>%
                                         layout(title = "PCA using methylation at motif-containing low-methylated regions (LMRs)", scene = list(xaxis = list(title = 'PC1 (43% variance)'),
                                                                                                                                                yaxis = list(title = 'PC2 (27% variance)'),
                                                                                                                                                zaxis = list(title = 'PC3 (22% variance)'))),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       plot_ly(DMR_methyl_principalComponents, x = ~PC1, y = ~PC2, z = ~PC3, color = ~tissue,colors=c("Brain"="blue","Liver"="Cyan","Lung"="magenta","Kidney"="green3","Cultured"="orange","AdultBrain"="darkblue")) %>%
                                         add_markers() %>%
                                         layout(title = "PCA using methylation at motif-containing differentially hypomethylated regions (DMRs)", scene = list(xaxis = list(title = 'PC1 (40% variance)'),
                                                                                                                                                               yaxis = list(title = 'PC2 (31% variance)'),
                                                                                                                                                               zaxis = list(title = 'PC3 (25% variance)')))
                                ),
                                column(6,
                                       h2("Explore high-variance transcription factor motifs across hypomethylated regions (LMRs)"),
                                       br(),
                                       selectizeInput("Methyl_tf_motif",label="Choose motif to display accessibility deviation Z-score",choices = methyl_highVariance_motifs$name,selected = "FoxL2(Forkhead)"),
                                       plotlyOutput("Methyl_PCA_motif"),
                                       br(),
                                       br(),
                                       h2("Explore high-variance transcription factor motifs across differentially hypomethylated regions (DMRs)"),
                                       br(),
                                       selectizeInput("DMR_Methyl_tf_motif",label="Choose motif to display accessibility deviation Z-score",choices = DMR_methyl_highVariance_motifs$name,selected = "Zic(Zf)"),
                                       plotlyOutput("DMR_Methyl_PCA_motif")
                                )
                              )
                            )
                   ),
                   tabPanel("Single-Cell RNA-seq",
                            fluidPage(
                              sidebarLayout(position="left",
                                            sidebarPanel(
                                              h4("Explore the developing brain EC single-cell expression dataset"),
                                              br(),
                                              typeaheadInput("sc_gene_name","Enter gene name",value="Bmx",choices = gene_list),
                                              actionButton("sc_submit_genes",label=h5("Submit"))
                                              
                                            ),
                                            mainPanel(
                                              plotOutput("tsne"),
                                              br(),
                                              img(src="tsne.png")
                                            )
                                            
                              )
                            )
                   )
)

server <- function(input, output, session) {
  MS_genes_to_plot <- eventReactive(input$MS_submit_genes,input$MS_gene_name,ignoreNULL = F)
  
  ECTSG_toShow <- reactive (input$ECTSG)
  
  ATAC_motifs_to_plot <- reactive(input$ATAC_tf_motif)
  differentialATAC_motifs_to_plot <- reactive(input$differentialATAC_tf_motif)
  
  Methyl_motifs_to_plot <- reactive(input$Methyl_tf_motif)
  DMR_Methyl_motifs_to_plot <- reactive(input$DMR_Methyl_tf_motif)
  
  output$chosen_gene <- renderTable(
    {
      dplyr::filter(EC_TPMs_averages,Gene == MS_genes_to_plot())
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
  
  output$MS_expression_plot <- renderPlot(
    {
      ggplot(filter(EC_TPMs_tidy,Gene %in% MS_genes_to_plot()),aes(x=Tissue,y=Expression,color=Tissue,shape=Fraction,alpha=Replicate)) +
        geom_point(size=4,position = position_dodge(width=0.9))+ 
        theme_gray(base_size = 20) + ylab("Expression (TPM)")+ 
        scale_color_manual(values=c("Brain"="blue","Liver"="Cyan","Lung"="magenta","Kidney"="green3","Cultured"="orange","AdultBrain"="darkblue"),guide=F)+
        expand_limits(y=0) + scale_alpha_manual(values=c("R1"=1,"R2"=1),guide=F) + scale_shape_discrete(name="Cell or tissue fraction",labels=c("Non-EC","EC","Total tissue"))+ 
        facet_grid(Gene ~ Fraction,scales="free_x", labeller = labeller(Fraction = as_labeller(fraction_names))) + theme(axis.text.x = element_text(angle = 90),strip.background = element_blank(), strip.text = element_text(face = "italic"))
      
    })
  
  
  
  output$RNA_heatmap <- renderPlotly( {
    
    
    
    plot_ly(RNA_principalComponents,hoverinfo = "text", x= ~PC1, y= ~PC2, z = ~PC3, text=~tissue, marker=list(color= ~log2(get(MS_genes_to_plot()) + 1), colorscale = "YlOrRd",cmin=0,cmax= ~max(log2(get(MS_genes_to_plot()) + 1)),showscale=TRUE,reversescale=T)) %>% 
      add_markers() %>% 
      layout(title = MS_genes_to_plot(), scene = list(xaxis = list(title="PC1 (26% variance)"),
                                                      yaxis = list(title="PC2 (20% variance)"),
                                                      zaxis = list(title="PC3 (15% variance)")),
             annotations = list(
               x = 1.13,
               y = 1.05,
               text = 'Log2 (TPM + 1)',
               xref = 'paper',
               yref = 'paper',
               showarrow = FALSE
             ))
  })
  
  output$ATAC_PCA_motif <- renderPlotly( {
    
    plot_ly(ATAC_principalComponents,hoverinfo = "text",text= ~tissue, x = ~PC1, y = ~PC2, z = ~PC3, marker = list(color = ~get(ATAC_motifs_to_plot()), colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'PC1 (30% variance)'),
                          yaxis = list(title = 'PC2 (25% variance)'),
                          zaxis = list(title = 'PC3 (18% variance)')),
             annotations = list(
               x = 1.13,
               y = 1.05,
               text = 'Z-score',
               xref = 'paper',
               yref = 'paper',
               showarrow = FALSE
             ))
    
    
    
  }
  )
  
  output$Methyl_PCA_motif <- renderPlotly( {
    
    plot_ly(methyl_principalComponents,hoverinfo="text",text=~tissue, x = ~PC1, y = ~PC2, z = ~PC3, marker = list(color = ~get(Methyl_motifs_to_plot()), colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'PC1 (43% variance)'),
                          yaxis = list(title = 'PC2 (27% variance)'),
                          zaxis = list(title = 'PC3 (22% variance)')),
             annotations = list(
               x = 1.13,
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
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'PC1 (40% variance)'),
                          yaxis = list(title = 'PC2 (31% variance)'),
                          zaxis = list(title = 'PC3 (25% variance)')),
             annotations = list(
               x = 1.13,
               y = 1.05,
               text = 'Z-score',
               xref = 'paper',
               yref = 'paper',
               showarrow = FALSE
             ))
    
    
  }
  )
  
  output$differentialATAC_PCA_motif <- renderPlotly( {
    
    plot_ly(ATAC_principalComponents_differentialPeaks,hoverinfo="text",text=~tissue, x = ~PC1, y = ~PC2, z = ~PC3, marker = list(color = ~get(differentialATAC_motifs_to_plot()), colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'PC1 (30% variance)'),
                          yaxis = list(title = 'PC2 (25% variance)'),
                          zaxis = list(title = 'PC3 (18% variance)')),
             annotations = list(
               x = 1.13,
               y = 1.05,
               text = 'Z-score',
               xref = 'paper',
               yref = 'paper',
               showarrow = FALSE
             ))
    
    
  }
  )
  
  sc_genes_to_plot <- eventReactive (input$sc_submit_genes, input$sc_gene_name, ignoreNULL = F)
  
  output$tsne <- renderPlot( {
    visualize_gene_markers2(gbm_log[, colnames(gbm_log) %in% rownames(tsne_proj3)], sc_genes_to_plot() ,tsne_proj3[,c("tsne3.1","tsne3.2")],limits=c(0,10), low_col = "lightblue", high_col = "darkblue", marker_size=0.6, axis_line_size = 0.5, axis_tick_size = 0.5) +
      xlab("TSNE 1") + ylab("TSNE 2") + theme(aspect.ratio = 1,text = element_text(size=20), legend.title = element_text(colour="black", size=18, face="bold"))
  }
  )
  
}


# Run the application 
shinyApp(ui = ui, server = server)

