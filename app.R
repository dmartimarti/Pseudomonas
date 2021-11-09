#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(openxlsx)
library(cowplot)
library(plotly)
library(glue)
library(ggtext)

data_long = read_csv('data/gene_counts_batch_normalized.csv')
stats_shiny = read_csv('data/complete_stats.csv') %>%
    select(gene_id, gene_name, Contrast_description, Contrast, 
           Media:Reference, baseMean:Direction)


theme_set(theme_cowplot(15))



# Define UI for application that draws a histogram
ui <- fluidPage(
    title = "Normalised counts RNA seq",
    
        sidebarPanel(
            selectizeInput('gene', choices = NULL, label = 'Gene',
                           multiple = TRUE, 
                           options = list(maxItems = 4, create = TRUE)),
            sliderInput("height", label = "Plot height", 
                        min = 100, max = 900, value = 500),
            sliderInput("width", label = "Plot width", 
                        min = 100, max = 1500, value = 500),
            br(),
            h3('For volcano'),
            selectInput('contrast', label = 'Contrast', 
                        choices = unique(stats_shiny$Contrast)),
            h4('Use these two sliders to control for extreme values and plot density.'), 
            sliderInput("padj_thres", label = "Max -(log Padj)", 
                        min = 0.1, max = 80, value = 15),
            sliderInput("log2_thres", label = "Max log2FC", 
                        min = 0.1, max = 50, value = 3.5)
    ),
    
        mainPanel(
            tabsetPanel(
                tabPanel("Boxplot",
                    h3('Normalised and corrected counts for batch effect'),
                    h4('Up to 4 genes to be plotted together'),
                    plotOutput("plot"),
                    br(),
                    br(),
                    br(),
                    br(),
                    br(),
                    br(),
                    br(),
                    br(),
                    br(),
                    br(),
                    br(),
                    dataTableOutput("stats")
                ),
                tabPanel("Volcano plot",
                    h3('Volcano plot'),
                    br(),
                    plotlyOutput(outputId  = 'VolcanoPlot')), 
                    br(),
            )
        )
            

)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    # update select size input
    updateSelectizeInput(session, 'gene', label = 'Gene',
                         choices = unique(data_long$gene_name),
                         selected = 'argk-1',
                         server = TRUE)
    
    output$plot <- renderPlot(
        
        data_long %>%
            mutate(Replicate = factor(Replicate,
                                      levels=c(1,2,3,4))) %>% 
            filter(gene_name %in% input$gene) %>%
            ggplot(aes(y = counts, x = Sample)) +
            geom_boxplot(aes(fill = Sample),
                         outlier.shape = NA) +
            geom_point(position = position_jitter(width = 0.2),
                       aes(shape = Replicate)) +
            facet_wrap(~gene_name, scales = 'free_y') +
            labs(x = 'Sample',
                 y = 'Normalised counts') +
            theme_cowplot(15) +
            panel_border() +
            guides(fill = 'none') +
            theme(axis.text.x = element_text(angle=45, hjust = 1)),
        
        width = function() input$width,
        height = function() input$height
    ) 
    
    output$VolcanoPlot = renderPlotly(

        {
            # prep data for volcano plot
            
            stats_shiny_volcano = stats_shiny %>% 
                dplyr::filter(Contrast %in% input$contrast) %>% 
                mutate(padj_log = -log10(padj),
                       significant = ifelse(padj < 0.05, 'Significant', ''),
                       log2FC = ifelse(abs(log2FoldChange) > input$log2_thres, input$log2_thres, log2FoldChange),
                       logPval = ifelse(padj_log > input$padj_thres, input$padj_thres, padj_log)) %>% 
                drop_na(significant) 
            
            # specify plot title and axes names
            plot_title = unique(stats_shiny_volcano$Contrast_description)
            target = unique(stats_shiny_volcano$Target)
            ref = unique(stats_shiny_volcano$Reference)

            p = stats_shiny_volcano %>%
                ggplot(aes(x = log2FC, y = logPval, 
                           text = paste('Gene: ',gene_name))) +
                geom_point(aes(color = significant), alpha = 0.7) +
                scale_colour_manual(values=c('#BABABA', '#0000FF')) +
                labs(x = glue::glue('log<sub>2</sub> Fold Change<br>Contrast: {target} vs {ref}'),
                     y = '-log<sub>10</sub>(P-value adj)',
                     title = plot_title,
                     color = NULL) +
                theme_cowplot(15) +
                theme(axis.title.x = element_markdown(),
                      axis.title.y = element_markdown()) +
                guides(color = 'none')
            
            
            ggplotly(p, height=800)
            
            }
        )
    
    
    output$stats = renderDataTable(stats_shiny %>% filter(gene_name %in% 
                                                             input$gene),
                                    options = list(pageLength = 20))
}

# Run the application 
shinyApp(ui = ui, server = server)
