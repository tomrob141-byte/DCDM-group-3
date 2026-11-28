library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)

#----- Importing data
#remove all NA p values
#assign NA values in group and parameter_stable_id into Unclassified and Not recorded
#do log transformation for p values, set cap for p_value = 0 into to avoid infinity  
dat <- read.csv("Cleaned_analysis.csv", na.strings = c("", "NA")) %>%
  filter(!is.na(p_value)) %>%
  mutate(group = coalesce(group, "Unclassified"),
         parameter_stable_id = coalesce(parameter_stable_id, "Not recorded"),
         log_p = -log10(pmax(p_value, 1e-85)))

#significant threshold, based on IMPC guidelines
threshold <- 0.0001   

#short names so the labels fit under the plot
short_names <- c(
  "Blood, clinical chemistry and metabolism"   = "Blood",
  "Brain, behaviour and neurological function" = "Brain",
  "Cardiovascular and respiratory function"    = "Heart",
  "Histopathology and Tissue Analysis"         = "Histology",
  "Imaging and morphology"                     = "Imaging",
  "Immune system and immunophenotyping"        = "Immune",
  "Reproductive and developmental biology"     = "Reproduction",
  "Weight and body composition"                = "Weight",
  "Unclassified"                               = "Unclassified"
)

## Assigning variables
groups <- names(short_names)
n_groups <- length(groups)
genes <- sort(unique(dat$gene_symbol))   

#------ Page Layout
ui <- fluidPage(
titlePanel("IMPC Mouse Phenotype Dashboard"),
  
h3("1. Phenotype scores for a knockout mouse"),
  selectizeInput("gene", "Select gene",
                 choices = genes, selected = "Marchf6",
                 options = list(maxOptions = length(genes))),   #allow all avail genes to appear
  plotlyOutput("plot1", height = "500px"),
  
  h4("Results ranked from most to least significant"),
  div(style = "height: 400px; overflow: auto;",   # makes the table scroll
      tableOutput("table1")),
  hr(),
  
h3("2. All knockout mice for a selected phenotype"),
  p("Coming soon"),
  hr(),
  
h3("3. Genes with similar phenotype scores"),
  p("Coming soon")
)

#----- Server
server <- function(input, output, session) {
  
  # Data for the chosen gene
  selected_data <- reactive({                #allow immediate respond when switching genes
    req(input$gene)                          #checks whether required input is selectec
    filter(dat, gene_symbol == input$gene)  #filter data of only the selected genes
  })
  
  #extract data from selected gene into variable called gene_data
  #group the data by phenotype group
 
  output$plot1 <- renderPlotly({
    #arrange parameter within each group alphabetically
    #create variable called position for each parameter using the formula : 
    # match(group,groups) - 0.5 + row_number / n+1. 
    # ex: Blood is the first group, RBC is 2nd parameter in blood group
    # there are 3 parameter in the blood group.  
    # 1 - 0.5 + 2/4 = 1.00
    # this allows 0.25 intervals between parameters
    #create variable called significant to store whether p value < 0.0001
    #create variable called info to store information to display during hover
    gene_data <- selected_data() %>%
      group_by(group) %>%
      arrange(parameter_name, .by_group = TRUE) %>%
      mutate(position = match(group, groups) - 0.5 + row_number() / (n() + 1)) %>%
      ungroup() %>%
      mutate(significant = p_value < threshold,
             info = paste0("Parameter: ", parameter_name,
                           "<br>Stable ID: ", parameter_stable_id,
                           "<br>p-value: ", signif(p_value, 3)))
  
    line_y <- -log10(threshold)          # assign line for p = 0.0001
    top <- max(gene_data$log_p, line_y) * 1.1 + 0.5  # calculate upper limit of graph, changes for each gene
    shaded <- seq(1, n_groups, by = 2)   # create sequence 1, 3, 5 for shading purpose
    
    #Generate plot
    g <- ggplot(gene_data, aes(x = position, y = log_p)) +
      #create shaded rectangles from position of xmin to xmax and ymin to ymax
      annotate("rect", xmin = shaded - 0.5, xmax = shaded + 0.5,
               ymin = 0, ymax = top, fill = "grey92") +
      
      #create dashed line for significance threshold, and add label
      geom_hline(yintercept = line_y, linetype = "dashed") +
      annotate("text", x = 1, y = line_y + top * 0.03,
               label = "p = 0.0001") +
      
      #generate dots, assign color based on significancem and add info to each dot for hover
      geom_point(aes(colour = significant, text = info), size = 1.5) +
      scale_colour_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey60")) +
      
      #edit labels in x and y axis
      #create tick marks, put label of group name, add spacing around each
      scale_x_continuous(breaks = 1:n_groups, labels = unname(short_names),
                         limits = c(0.5, n_groups + 0.5), expand = c(0, 0)) +
      labs(x = "Parameters, in blocks by group", y = "-log10(p-value)") +
      theme_bw() +
      #remove legend explaining dots and background grid lines
      theme(legend.position = "none", panel.grid = element_blank())
    
    #converts ggplot to interactive Ploty, and allow display of text on hover
    ggplotly(g, tooltip = "text")
  })
  
  #creates a table that list p value from smallest
  output$table1 <- renderTable({
    selected_data() %>%
      arrange(p_value) %>%                      # smallest p-value first
      transmute("Rank" = row_number(),
                "Significant" = ifelse(p_value < threshold, "Yes", "No"),
                "p-value" = sprintf("%.3g", p_value),
                "Parameter" = parameter_name,
                "Stable ID" = parameter_stable_id,
                "Group" = group,
                "Procedure" = suffix_desc,
                "Allele" = allele_symbol,
                "Strain" = strain_name,
                "Life stage" = mouse_life_stage)
    #striped alternates row background colors
    #hover highlights a row when mouse pointer moves over it 
    #na = "Not recorde" displayes "Not recorded" whenever a value is missing
  }, striped = TRUE, hover = TRUE, na = "Not recorded")
}

shinyApp(ui, server)