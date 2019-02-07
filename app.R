#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Dehydration experiment Data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        img(src="rstudio.png", height=80,width=150), br(),
         selectInput("yvalues",
                     "Y-axis choice:",
                     choices=c("Total weight loss"="total.loss",
                               "Total percent body weight lost"="total.percent",
                               "Average body weight lost per day"="avg.daily",
                               "Average percent body weight lost per day"="avgpercentdaily",
                               "Number of days survived at 0% relative humidity"="days.survived"))
      ),
      
      mainPanel(
         plotOutput("waterplot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$waterplot <- renderPlot({
      # generate bins based on input$bins from ui.R
      if(input$yvalues=="total.loss"){
        y=waterstats$total.loss
        ylabs="Total weight lost"
      }
      else if (input$yvalues=="total.percent"){
        y=waterstats$total.percent
      ylabs="Total percent body weight lost"
      }
      else if (input$yvalues=="avg.daily"){
        y=waterstats$avg.daily
      ylabs="Total percent body weight lost"
      }
     else if (input$yvalues=="avgdailypercent"){
       y=waterstats$avgdailypercent
       ylabs="avg percent body weight lost per day"
     }
     else{
       y=waterstats$days.survived
       ylabs="Days survived dried"
     }
      # make ggplot obj
      ggplot(waterstats,aes(factor(status),y))+
        geom_boxplot()+
        stat_boxplot(geom="errorbar")+
        ylab(ylabs)+
        xlab("Treatment group")
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

