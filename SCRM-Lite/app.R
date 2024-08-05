#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinyjs)
library(bslib)
library(shinyWidgets)
library(ggplot2)

page_01 <- fluidPage(
    useShinyjs(),
    sidebarLayout(
        sidebarPanel(width = 2,
                     # style = "position:fixed;",
                     fluidRow(
                         column(
                             sliderTextInput("sldtxt_variability", "Contamination Event Variability", choices = list("Low", "High"), "High", FALSE, TRUE, hide_min_max = TRUE),
                             h4("Scenario Inputs"),
                             sliderTextInput("sldtxt_additional_testing", "Additional Product Testing", choices = list("None", "Some"), "None", FALSE, TRUE, hide_min_max = TRUE),
                             sliderTextInput("sldtxt_process_controls", "Process Wash", choices = list("Standard", "Improved"), "Standard", FALSE, TRUE, hide_min_max = TRUE),
                             width = 12, align="center"
                         ),
                     ),
        ),
        mainPanel(
            layout_columns(
                plotOutput("plot_overall_risk"),
                # plotOutput("plot_n_highest_risk_lots"),
                col_widths = c(8)
            ),
            layout_columns(
                plotOutput("plot_n_highest_risk_lots"),
                col_widths = c(8)
            ),
            width = 10,
        )
    )
)

# Page 02 ----
page_02 <- fluidPage(
    useShinyjs(),
    tags$style(".uiucimg {
                            margin-left:65px;
                            margin-right:0px;
                            margin-top:65px;
                          }"),
    tags$style(".cornellimg {
                            margin-left:130px;
                            margin-right:65px;
                            margin-top:65px;
                          }"),
    titlePanel("Acknowledgements"),
    h3(""),
    h5("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."),
    h5(HTML("")),
    h3("Contacts"),
    h5("Please feel free to contact us if you have any questions."),
    h5(HTML("<b>Matthew Stasiewicz</b>: Principal Investigator, Associate Professor of Applied Food Safety, PhD | Email: "), a("mstasie@illinois.edu",
                                                                                                                               target = "_blank",
                                                                                                                               href = "mailto:mstasie@illinois.edu"
    )),
    h5(HTML("<b>Martin Wiedmann</b>: Co-Principal Investigator, Gellert Family Professor in Food Safety, PhD | Email: "), a("martin.wiedmann@cornell.edu",
                                                                                                                            target = "_blank",
                                                                                                                            href = "mailto:martin.wiedmann@cornell.edu"
    )),
    h5(HTML("<b>Cecil Barnett-Neefs</b>: App Creator, Model Author | Email: "), a("cecilwb2@illinois.edu",
                                                                                  target = "_blank",
                                                                                  href = "mailto:cecilwb2@illinois.edu"
    )),
    h5(HTML("<b>Gabriella Pinto</b>: Model Author | Email: "), a("gnpinto2@illinois.edu",
                                                                 target = "_blank",
                                                                 href = "mailto:gnpinto2@illinois.edu"
    )),
    h3("Acknowledgements"),
    fluidRow(
        # column(1, div(class="uiucimg", (imageOutput("img_uiuc")))),
        # column(1, div(class="cornellimg", (imageOutput("img_cornell"))))
        column(2, div(class="uiucimg", img(src = "University-Wordmark-Full-Color-RGB-TM.png", width = "250px", align = "left"))),
        column(2, div(class="cornellimg", img(src = "bold_cornell_logo_pms187_red.png", width = "250px", align = "left")))
    )
)

# UI ----
ui <- navbarPage(
    theme = bs_theme(bootswatch = "flatly", version = 5),
    title = "SCRM (in R) Lite",
    tabPanel("Model", page_01),
    tabPanel("Acknowledgements", page_02)
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    rvs <- reactiveValues()
    rvs$overall_risk <- 4500
    rvs$n_highest_risk_lots <- 16
    rvs$variability <- "High"
    rvs$process_controls <- "Standard"
    rvs$additional_testing <- "None"
    lv_df <- as.data.frame(cbind(
        c("Low", "Low", "Low", "Low"),
        c("Standard", "Standard", "Improved", "Improved"),
        c("None", "Some", "None", "Some"),
        c(20000, 21200, 113200, 115000),
        c(0, 0, 0, 0),
        c("Baseline", "Additional Product Testing", "Improved Process Controls", "Both Practices"),
        c(FALSE, FALSE, FALSE, FALSE)
    ))
    colnames(lv_df) <- c("variability", "process_controls", "additional_testing", "overall_risk", "n_highest_risk_lots", "practice", "selected")
    lv_df$overall_risk <- as.integer(lv_df$overall_risk)
    lv_df$n_highest_risk_lots <- as.integer(lv_df$n_highest_risk_lots)
    lv_df$practice <- factor(lv_df$practice, levels = c("Baseline", "Additional Product Testing", "Improved Process Controls", "Both Practices"))
    # lv_df$selected[which(lv_df$variability == rvs$variability & lv_df$process_controls == rvs$process_controls & lv_df$additional_testing == rvs$additional_testing)] <- TRUE
    rvs$results_dataframe_low_variability <- lv_df
    
    hv_df <- as.data.frame(cbind(
        c("High", "High", "High", "High"),
        c("Standard", "Standard", "Improved", "Improved"),
        c("None", "Some", "None", "Some"),
        c(4500, 10800, 25700, 35600),
        c(16, 0, 0, 0),
        c("Baseline", "Additional Product Testing", "Improved Process Controls", "Both Practices"),
        c(TRUE, FALSE, FALSE, FALSE)
    ))
    colnames(hv_df) <- c("variability", "process_controls", "additional_testing", "overall_risk", "n_highest_risk_lots", "practice", "selected")
    hv_df$overall_risk <- as.integer(hv_df$overall_risk)
    hv_df$n_highest_risk_lots <- as.integer(hv_df$n_highest_risk_lots)
    hv_df$practice <- factor(hv_df$practice, levels = c("Baseline", "Additional Product Testing", "Improved Process Controls", "Both Practices"))
    rvs$results_dataframe_high_variability <- hv_df
    
    toListenAllInputs <- reactive({
        list(
            input$sldtxt_variability,
            input$sldtxt_process_controls,
            input$sldtxt_additional_testing
        )
    })
    
    observeEvent(toListenAllInputs(), {
        rvs$variability <- input$sldtxt_variability
        rvs$process_controls <- input$sldtxt_process_controls
        rvs$additional_testing <- input$sldtxt_additional_testing
        
        rvs$results_dataframe_low_variability$selected <- FALSE
        rvs$results_dataframe_low_variability$selected[which(rvs$results_dataframe_low_variability$variability == rvs$variability & rvs$results_dataframe_low_variability$process_controls == rvs$process_controls & rvs$results_dataframe_low_variability$additional_testing == rvs$additional_testing)] <- TRUE
        rvs$results_dataframe_high_variability$selected <- FALSE
        rvs$results_dataframe_high_variability$selected[which(rvs$results_dataframe_high_variability$variability == rvs$variability & rvs$results_dataframe_high_variability$process_controls == rvs$process_controls & rvs$results_dataframe_high_variability$additional_testing == rvs$additional_testing)] <- TRUE
    
        if (rvs$variability == "Low") {
            rvs$plot_dataframe <- rvs$results_dataframe_low_variability
            rvs$overall_risk <- rvs$results_dataframe_low_variability$overall_risk[which(rvs$results_dataframe_low_variability$selected == TRUE)]
            rvs$n_highest_risk_lots <- rvs$results_dataframe_low_variability$n_highest_risk_lots[which(rvs$results_dataframe_low_variability$selected == TRUE)]
        } else if (rvs$variability == "High") {
            rvs$plot_dataframe <- rvs$results_dataframe_high_variability
            rvs$overall_risk <- rvs$results_dataframe_high_variability$overall_risk[which(rvs$results_dataframe_high_variability$selected == TRUE)]
            rvs$n_highest_risk_lots <- rvs$results_dataframe_high_variability$n_highest_risk_lots[which(rvs$results_dataframe_high_variability$selected == TRUE)]
        }
        
        
    })
    
    output$plot_overall_risk <- renderPlot({
        ggplot(rvs$plot_dataframe, aes(x=practice, y=overall_risk)) +
            geom_bar(stat="identity", aes(fill=selected)) +
            geom_text(aes(label=overall_risk), vjust=-0.2, position = position_dodge(0.9), size=5) +
            xlab("") +
            ylab("The Overall Risk is 1 in...") +
            scale_fill_manual(values = c("grey", "green")) +
            theme_bw() +
            theme(
                plot.title = element_text(size = 20),
                axis.title = element_text(size = 15),
                axis.text = element_text(size = 15),
                legend.position = "none"
            ) +
            ggtitle(paste("Overall Risk of a Positive Test at Retail Under", rvs$variability, "Variability System"))
    })
    
    output$plot_n_highest_risk_lots <- renderPlot({
        ggplot(rvs$plot_dataframe, aes(x=practice, y=n_highest_risk_lots)) +
            geom_bar(stat="identity", aes(fill=selected)) +
            geom_text(aes(label=n_highest_risk_lots), vjust=-0.2, position = position_dodge(0.9), size=5) +
            xlab("") +
            ylab("") +
            scale_fill_manual(values = c("grey", "green")) +
            theme_bw() +
            ylim(0,20) +
            theme(
                plot.title = element_text(size = 20),
                axis.title = element_text(size = 15),
                axis.text = element_text(size = 15),
                legend.position = "none"
            ) +
            ggtitle(paste("Number of Highest Risk Lots Under", rvs$variability, "Variability System"))
    })
    
}
    
# Run the application 
shinyApp(ui = ui, server = server)