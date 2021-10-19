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

data_long = read_csv('data/gene_counts_batch_normalized.csv')
stats_shiny = read_csv('data/complete_stats.csv')

theme_set(theme_cowplot(15))

# Define UI for application that draws a histogram
ui <- fluidPage(
    title = "Normalised counts RNA seq",
    sidebarPanel(
        selectizeInput('gene', label = 'Gene', multiple = TRUE,
                       selected = 'aap-1',
                       choices = unique(data_long$gene_name),
                       options = list(maxItems = 4, create = TRUE)),
        # selectInput("gene", label = "Gene", 
        #             choices = unique(data_long$gene_name)),
        sliderInput("height", label = "Plot height", 
                    min = 100, max = 900, value = 500),
        sliderInput("width", label = "Plot width", 
                    min = 100, max = 1500, value = 500)
    ),
    mainPanel(
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
        tableOutput("mytable")
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
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
    
    output$mytable <- renderTable(stats_shiny %>% filter(gene_name %in% 
                                                             input$gene))
}

# Run the application 
shinyApp(ui = ui, server = server)
