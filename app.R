# Interactive look at the results tables. Run from the repo root:
#   Rscript -e 'shiny::runApp()'
# Needs shiny and ggplot2: install.packages(c("shiny", "ggplot2"))

library(shiny)
library(ggplot2)

counts <- read.csv("results/celltype_counts.csv")
markers <- read.csv("results/marker_genes.csv")

ui <- fluidPage(
    titlePanel("PBMC 3k results"),
    p("2,638 cells after QC, 8 clusters mapped to cell types from marker-gene scores.",
      "The full analysis is in the README and the notebook."),
    sidebarLayout(
        sidebarPanel(
            selectInput("cluster",
                        "Cluster",
                        choices = c("all", as.character(sort(unique(markers$cluster)))),
                        selected = "all"),
            helpText("Filters the marker table below.")
        ),
        mainPanel(
            imageOutput("umap"),
            plotOutput("counts", height = "300px"),
            h4("Marker genes"),
            tableOutput("markers")
        )
    )
)

server <- function(input, output) {
    output$umap <- renderImage({
        list(src = "figures/11_umap_celltypes.png",
             width = "100%", contentType = "image/png")
    }, deleteFile = FALSE)

    output$counts <- renderPlot({
        ggplot(counts, aes(x = n_cells, y = reorder(cell_type, n_cells))) +
            geom_col(fill = "#4C72B0") +
            labs(x = "number of cells", y = NULL) +
            theme_minimal(base_size = 13)
    })

    output$markers <- renderTable({
        m <- markers
        if (input$cluster != "all") {
            m <- m[m$cluster == as.numeric(input$cluster), ]
        }
        m <- m[order(m$cluster, m$rank), c("cluster", "rank", "gene",
                                           "log2_fold_change", "pval_adj")]
        m$log2_fold_change <- sprintf("%.2f", m$log2_fold_change)
        m$pval_adj <- format(m$pval_adj, digits = 2)
        names(m) <- c("cluster", "rank", "gene", "log2 FC", "p adj")
        head(m, 100)
    }, rownames = FALSE, digits = 2)
}

shinyApp(ui, server)
