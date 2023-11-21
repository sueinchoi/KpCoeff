library(shiny)
library(tidyverse)
library(gridExtra)
library(shinydashboard)
library(DT)
source("Kp_calculator/code/CalcKp_P&T.R")
source("Kp_calculator/code/CalcKp_R&R.R")
source("Kp_calculator/code/CalcKpu_R&R.R")
source("Kp_calculator/code/CalcKp_Berez.R")
source("Kp_calculator/code/CalcKp_Schmitt.R")
source("Kp_calculator/code/CalcKp_pksim.R")
source("Kp_calculator/code/CalcVss_P&T.R")
source("Kp_calculator/code/getPartitionCoeff.R")

dat_uni <- read.csv("Kp_calculator/data/unified_tissue_comp.csv") # unified tissue composition
dat_PT <- read.csv("Kp_calculator/data/tissue_comp_P&T.csv") # data reported by Poulin and Theil
dat_Berez <- read.csv("Kp_calculator/data/PKSim_tissue_comp_PT_Berez.csv") # data used by PK-Sim for PT and Berez methods
dat_RR <- read.csv("Kp_calculator/data/tissue_comp_R&R.csv") # data reported by Rodgers and Rowland
dat_Schmitt_rep <- read.csv("Kp_calculator/data/tissue_comp_Schmitt.csv") # data reported by Schmitt
dat_Schmitt <- read.csv("Kp_calculator/data/PKSim_tissue_comp_Schmitt.csv") # data used by PK-Sim for Schmitt method
dat_pksim <- read.csv("Kp_calculator/data/PKSim_tissue_comp_pksim.csv") # data used by PK-Sim for PK-Sim method

dat_PT_RR <- read.csv("Kp_calculator/data/tissue_comp_R&R_for_P&T.csv") # data reported by R&R formatted to use with PT method
dat_pksim_RR <- read.csv("Kp_calculator/data/tissue_comp_R&R_for_pksim.csv") # data reported by R&R formatted to use with PK-Sim method
dat_Schmitt_RR <- read.csv("Kp_calculator/data/tissue_comp_R&R_for_Schmitt.csv")



body <- dashboardBody(
    fluidPage(
        fluidRow(
            column(
                width = 2,
                radioButtons('fit_Kp', "Choose the prediction method",
                             choices = list("Poulin & Theil" = "P&T",
                                            "Berez" = "Berez",
                                            "Rowland & Rodgers" = "R&R",
                                            "Schmitt" = "Schmitt",
                                            "PK-sim" = "pksim")
                                ),
                radioButtons('fit_dat', "Choose data type",
                             choices = list("Original" = 1,
                                            "Unified" = 2)
                                ),               
                numericInput('logP', 'Partition coefficient (logP)', min = 0, max = 30, value = 1),
                numericInput('fup', 'Plasma unbound fraction (fu)', min = 0, max = 1, value = 0.9),
                numericInput('BP', 'Blood:plasma ratio (BP)', min = 0, max = 5, value = 1),
                radioButtons('fit_type', "Acid/Base selection",
                             choices = list("Neutral" = 1,
                                            "Monoprotic acid" = 2,
                                            "Monoprotic base" = 3,
                                            "Diprotic acid" = 4,
                                            "Diprotic base" = 5,
                                            "Zwitterion" = 6)
                                            ),
                conditionalPanel(
                    condition = "input.fit_type != 1",
                    numericInput('pKa1', 'pKa1', min = -20, max = 30, value = 1)
                ),
                conditionalPanel(
                    condition = "input.fit_type == 4 | input.fit_type == 5 | input.fit_type ==6",
                    numericInput('pKa2', 'pKa2', min = -20, max = 30, value = 1)
                )
            ),
            column(
                width = 7,
                # downloadButton("downloadData", "Parameter download"),
                numericInput('Kpscalar', 'Kp scalar', min = 0, max = 100, value = 1),
                dataTableOutput('Kp_result'),
                verbatimTextOutput('pred')
                # uiOutput('Vss')
            )
        )
    )
)

sidebar <- dashboardSidebar(
    disable = TRUE
)
ui <- dashboardPage(
    dashboardHeader(
        title = "Kp prediction"
    ),
    sidebar = sidebar,
    body)

`%nin%` = Negate(`%in%`)
server <- function(input, output){
    
    #output the datatable based on the dataframe (and make it editable)

    data <- reactive({
        if(input$fit_dat == 1){
            switch(input$fit_Kp,
                "P&T" = dat_PT,
                "Berez" = dat_Berez,
                "R&R" = dat_RR,
                "Schmitt" = dat_Schmitt,
                "pksim" = dat_pksim
            )
        } else {
            dat_uni
        }
    })

    output$pred <- renderPrint({
            data()
    })
    
    Kp <- reactive({
        if(input$fit_type %in% c(1, 2, 3))
        {
            pKa_vec <- ifelse(input$fit_type == 1, 0, input$pKa1)
        } else {
            pKa_vec <- sort(c(input$pKa1, input$pKa2))
        }
        
        result <- pcoeffs(input$logP, pKa_vec, input$fup, input$BP, as.numeric(input$fit_type), input$fit_Kp, data()) # nolint
        map_df(result, ~.x, .names = "id") %>%
            gather("tissue", "Kp") %>%
            mutate(Tissue = data() %>% filter(tissue %nin% c('RBC', 'Plasma', 'RBCs', 'plasma')) %>% pull(tissue), Kp = signif(Kp, 4)) %>%
            select(Tissue, Kp)
    })
    output$Kp_result <- renderDataTable({
       Kp()
    })


    # output$Vss <- renderUI({
    #     HTML(paste0("<b>", "Vss (L/kg) = ", as.character(Vss_raw()), "</b>"))
    # })

    # Vss_raw <- reactive({
    #     Vss <- Kp() %>%
    #         left_join(dat_PT, by = "tissue") %>%
    #         mutate(Vss = Kp * FV) %>%
    #         pull(Vss)
    #     sum(Vss)
    # })
    # values <- reactive({
    #     value1 <- data.frame(Parameter = c('Model', 'ka', 'Fa', 'Kpliver', 'Vc_kg', 'Vp_kg', 'Q_kg', 'fu', 'CLint_noncyp', 'CLscalar', 'CLint_cyp'),
    #                          Values = c(ifelse(input$fit_PK==1, '1comp', '2comp'),
    #                                     input$ka, input$fa, input$Kpliver, input$Vc, 
    #                                     ifelse(input$fit_PK == 1, NA, input$Vp), 
    #                                     ifelse(input$fit_PK == 1, NA, input$Q), 
    #                                     input$fu,  input$CLs,  input$CL, CLint_v()))
    #     colnames(v$data) <- c('Parameter', 'Values')
    #     rbind(value1, v$data) %>%
    #         cbind(unit = c(" ", "/h", " ", " ", "L/kg", "L/kg", "L/h/kg", " ", "uL/min/mg ptn", " ", rep("uL/min/mg ptn", "uL/min/mg ptn", 9)))
    # })
    
#     output$downloadData <- downloadHandler(
#         filename = function(){
#             "Parameters.csv"
#         },
#         content = function(file){
#             write.csv(values(), file, row.names = FALSE)
#         },
#         contentType = "text/csv"
#     )
    
 }


shinyApp(ui = ui, server = server)

