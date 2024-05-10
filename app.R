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

stats_int = read_csv('data/interaction_stats.csv')


theme_set(theme_cowplot(15))

data_groups = c('WT_0', 'WT_50', 'WTN_0', 'WTN_50',
                'B_0', 'B_50', 'G_0', 'OP50')

# Define UI for application that draws a histogram
ui <- fluidPage(
    title = "Normalised counts RNA seq",
    
        sidebarPanel(
            # general
            selectizeInput('gene', choices = NULL, label = 'Gene',
                           multiple = TRUE, 
                           options = list(maxItems = 4, create = TRUE)),
            selectizeInput('data_groups', label= 'Groups to plot',
                        choices = data_groups, multiple = TRUE),
            sliderInput("height", label = "Plot height", 
                        min = 100, max = 900, value = 500),
            sliderInput("width", label = "Plot width", 
                        min = 100, max = 1500, value = 500),
            br(),
            # volcano plot
            h3('For volcano'),
            selectInput('contrast', label = 'Contrast', 
                        choices = unique(stats_shiny$Contrast)),
            h4('Use these two sliders to control for extreme values and plot density.'), 
            sliderInput("padj_thres", label = "Max -(log Padj)", 
                        min = 0.1, max = 80, value = 15),
            sliderInput("log2_thres", label = "Max log2FC", 
                        min = 0.1, max = 50, value = 3.5),
            br(),
            # interaction
            h3('For interaction terms'),
            h4('Select here the genes you want to plot'),
            selectizeInput('gene_int', choices = NULL, label = 'Gene int',
                           multiple = F, 
                           options = list(create = TRUE))
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
                    dataTableOutput("stats"),
                    downloadButton('downloadData', 'Download data'),
                    # download button for plot in pdf
                    downloadButton('downloadStats', 'Download stats')
                ),
                
                tabPanel("Volcano plot",
                    h3('Volcano plot'),
                    br(),
                    plotlyOutput(outputId  = 'VolcanoPlot')
                    
                ), 
                    
                
                tabPanel('Interaction',
                         h3('Interaction plots'),
                         plotOutput('plot_int'),
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
                         dataTableOutput("stats_int")
                )        
            
        ),
            
                
        )
            

)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    # update select size input for gene
    updateSelectizeInput(session, 'gene', label = 'Gene',
                         choices = unique(data_long$gene_name),
                         selected = 'argk-1',
                         server = TRUE)
  
    # update select size input for gene
    updateSelectizeInput(session, 'data_groups', label = 'Groups to plot',
                         choices = data_groups,
                         selected = data_groups,
                         server = TRUE)
    
    # update select size input for gene_int
    updateSelectizeInput(session, 'gene_int', label = 'Gene int',
                         choices = unique(data_long$gene_name),
                         selected = 'argk-1',
                         server = TRUE)
    
    output$plot = renderPlot(
        
        data_long %>%
            mutate(Replicate = factor(Replicate,
                                      levels=c(1,2,3,4))) %>% 
            filter(gene_name %in% input$gene) %>%
            filter(Sample %in% input$data_groups) %>% 
            ggplot(aes(y = counts, x = Sample)) +
            geom_boxplot(aes(fill = Sample),
                         outlier.shape = NA) +
            geom_point(position = position_jitter(width = 0.2),
                       aes(shape = Replicate)) +
            facet_wrap(~gene_name*gene_id, scales = 'free_y') +
            labs(x = 'Sample',
                 y = 'Normalised counts') +
            theme_cowplot(15) +
            panel_border() +
            guides(fill = 'none') +
            theme(axis.text.x = element_text(angle=45, hjust = 1)),
        
        width = function() input$width,
        height = function() input$height
    ) 
    
    #### Volcano plot ####
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
    
    
    #### interaction plot ####
    
    # output$plot_int = renderPlot(
    #     
    #     gns = stats_int %>% filter(gene_id %in% input$gene_int),
    #     
    #     gns_stats = stats_int %>% 
    #         filter(gene_id %in% gns) %>% 
    #         mutate(p_adj_stars = ifelse(p_adj_stars==' ','no sig',p_adj_stars)),
    #     
    #     gns_counts = data_long %>%
    #         dplyr::filter(gene_id %in% gns) %>%
    #         filter(Sample %in% c('WT_0','WT_50','WTN_0','WTN_50')) %>% 
    #         mutate(max_val = max(counts)),
    #     
    #     max_val = gns_counts$max_val[1],
    #     
    #     fc.wt = gns_stats %>% filter(contrast=='WT') %>% pull(log2FoldChange) %>% round(2),
    #     fc.wtn = gns_stats %>% filter(contrast=='WTN') %>% pull(log2FoldChange)%>% round(2),
    #     fc.int = gns_stats %>% filter(contrast=='interaction') %>% pull(log2FoldChange)%>% round(2),
    #     
    #     pval.wt = gns_stats %>% filter(contrast=='WT') %>% pull(p_adj_stars),
    #     pval.wtn = gns_stats %>% filter(contrast=='WTN') %>% pull(p_adj_stars),
    #     pval.int = gns_stats %>% filter(contrast=='interaction') %>% pull(p_adj_stars),
    #     
    #     
    #     gns_counts %>% 
    #         ggplot(aes(y = counts, x = Sample)) +
    #         geom_boxplot(aes(fill = Sample),
    #                      outlier.colour = NULL,
    #                      outlier.shape = NA) +
    #         geom_point(position = position_jitter(width = 0.2)) +
    #         stat_summary(
    #             geom = "point",
    #             fun = "mean",
    #             col = "black",
    #             size = 6,
    #             shape = 21,
    #             fill = "blue"
    #         ) +
    #         # WT segment
    #         geom_segment(aes(x = 1, y = max_val*1.05, 
    #                          xend = 2, yend = max_val*1.05),
    #                      colour = 'black', size = 1) +
    #         annotate('text', x = 1.5, y = max_val*1.07,
    #                  label = glue('log2FC:{fc.wt}, pval:{pval.wt}')) +
    #         # WTN segment
    #         geom_segment(aes(x = 3, y = max_val*1.05, 
    #                          xend = 4, yend = max_val*1.05),
    #                      colour = 'black', size = 1) +
    #         annotate('text', x = 3.5, y = max_val*1.07,
    #                  label = glue('log2FC:{fc.wtn}, pval:{pval.wtn}')) +
    #         # Int segment
    #         geom_segment(aes(x = 1.5, y = max_val*1.15, 
    #                          xend = 3.5, yend = max_val*1.15),
    #                      colour = 'red', size = 1) +
    #         annotate('text', x = 2.5, y = max_val*1.17,
    #                  label = glue('log2FC:{fc.int}, pval:{pval.int}')) +
    #         facet_wrap(~gene_name, scales = 'free_y') +
    #         labs(x = 'Sample',
    #              y = 'Normalised counts (log scale)') +
    #         theme_cowplot(15) +
    #         panel_border() +
    #         theme(axis.text.x = element_text(angle=45, hjust = 1))
    # 
    # )
    
    
    
    output$stats = renderDataTable(stats_shiny %>% filter(gene_name %in% 
                                                             input$gene),
                                    options = list(pageLength = 20))
    
    # create a button to download the data from the selected genes from data_long
    output$downloadData = downloadHandler(
        filename = function() {
            paste('data-', Sys.Date(), '.csv', sep='')
        },
        content = function(file) {
            write.csv(data_long %>% filter(gene_name %in% input$gene), file)
        }
    )
    
  # download the stats from the selected genes as a csv
    output$downloadStats = downloadHandler(
      filename = function() {
        paste('stats-', Sys.Date(), '.csv', sep='')
      },
      content = function(file) {
        write.csv(stats_shiny %>% filter(gene_name %in% input$gene), file)
    })

  }

# Run the application 
shinyApp(ui = ui, server = server)
