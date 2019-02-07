#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/

library(shiny)
library(DESeq2)
library(DBI)
library(RSQLite)
library(reshape2)
library(plotly)
library(dplyr)

#Create a connection to the SQLite database (change path if needed !)
finalproj = dbConnect(dbDriver("SQLite"), dbname='finalproj.db')

dbExecute( finalproj, "CREATE TABLE  IF NOT EXISTS countstable AS
SELECT NPF1.CDS, NPF1.NP1, NPF2.NP2, NPF3.NP3, PreL1.PL1, PreL2.PL2, PreL3.PL3,	EarL1.EL1,EarL2.EL2,EarL3.EL3,LateL1.LL1,	LateL2.LL2,LateL3.LL3,Male1.M1,Male2.M2,Male3.M3
           from NPF1
           JOIN NPF2
           ON NPF1.CDS=NPF2.CDS
           JOIN NPF3
           ON NPF1.CDS=NPF3.CDS
           JOIN PreL1
           ON NPF1.CDS=PreL1.CDS
           JOIN PreL2
           ON NPF1.CDS=PreL2.CDS
           JOIN PreL3
           ON NPF1.CDS=PreL3.CDS
           JOIN EarL1
           ON NPF1.CDS=EarL1.CDS
           JOIN EarL2
           ON NPF1.CDS=EarL2.CDS
           JOIN EarL3
           ON NPF1.CDS=EarL3.CDS
           JOIN LateL1
           ON NPF1.CDS=LateL1.CDS
           JOIN LateL2
           ON NPF1.CDS=LateL2.CDS
           JOIN LateL3
           ON NPF1.CDS=LateL3.CDS
           JOIN Male1
           ON NPF1.CDS=Male1.CDS
           JOIN Male2
           ON NPF1.CDS=Male2.CDS
           JOIN Male3
           ON NPF1.CDS=Male3.CDS")


countframe=data.frame(dbReadTable(finalproj,"countstable"))
countframe$NP1=as.integer(countframe$NP1)
countframe$NP2=as.integer(countframe$NP2)
countframe$NP3=as.integer(countframe$NP3)
countframe$PL1=as.integer(countframe$PL1)
countframe$PL2=as.integer(countframe$PL2)
countframe$PL3=as.integer(countframe$PL3)
countframe$EL1=as.integer(countframe$EL1)
countframe$EL2=as.integer(countframe$EL2)
countframe$EL3=as.integer(countframe$EL3)
countframe$LL1=as.integer(countframe$LL1)
countframe$LL2=as.integer(countframe$LL2)
countframe$LL3=as.integer(countframe$LL3)
countframe$M1=as.integer(countframe$M1)
countframe$M2=as.integer(countframe$M2)
countframe$M3=as.integer(countframe$M3)
sampleframe=data.frame(dbReadTable(finalproj,"sampledat"))
countmatrix=as.matrix(countframe[,-1])
rownames(countmatrix)=countframe[,1]
samplematrix=as.matrix(sampleframe[,-1])
rownames(samplematrix)=sampleframe[,1]

desdatset <- DESeqDataSetFromMatrix(countData = countmatrix,
                                    colData = samplematrix,
                                    design = ~Condition)

desdatset <- DESeq(desdatset)

NPvPLres=results(desdatset,contrast=c("Condition","NP","PL"))
NotPvFirstT=data.frame(NPvPLres,row.names=NULL)
NotPvFirstT$CDS=rownames(NPvPLres)
NPvELres=results(desdatset,contrast=c("Condition","NP","EL"))
NotPvSecondT=data.frame(NPvELres,row.names=NULL)
NotPvSecondT$CDS=rownames(NPvELres)
NPvLLres=results(desdatset,contrast=c("Condition","NP","LL"))
NotPvThirdT=data.frame(NPvLLres,row.names=NULL)
NotPvThirdT$CDS=rownames(NPvLLres)
NPvMres=results(desdatset,contrast=c("Condition","NP","M"))
NotPvMale=data.frame(NPvMres,row.names=NULL)
NotPvMale$CDS=rownames(NPvMres)
PLvELres=results(desdatset,contrast=c("Condition","PL","EL"))
FirstTvsSecondT=data.frame(PLvELres,row.names=NULL)
FirstTvsSecondT$CDS=rownames(PLvELres)
PLvLLres=results(desdatset,contrast=c("Condition","PL","LL"))
FirstTvsThirdT=data.frame(PLvLLres,row.names=NULL)
FirstTvsThirdT$CDS=rownames(PLvLLres)
PLvMres=results(desdatset,contrast=c("Condition","PL","M"))
FirstTvsMale=data.frame(PLvMres,row.names=NULL)
FirstTvsMale$CDS=rownames(PLvMres)
ELvLLres=results(desdatset,contrast=c("Condition","EL","LL"))
SecondTvsThirdT=data.frame(ELvLLres,row.names=NULL)
SecondTvsThirdT$CDS=rownames(ELvLLres)
ELvMres=results(desdatset,contrast=c("Condition","EL","M"))
SecondTvsMale=data.frame(ELvMres,row.names=NULL)
SecondTvsMale$CDS=rownames(ELvMres)
LLvMres=results(desdatset,contrast=c("Condition","LL","M"))
ThirdTvsMale=data.frame(LLvMres,row.names=NULL)
ThirdTvsMale$CDS=rownames(LLvMres)

countframe$NPvFT.padj=NotPvFirstT$padj
countframe$NPvFT.FC=NotPvFirstT$log2FoldChange

countframe$NPvST.padj=NotPvSecondT$padj
countframe$NPvST.FC=NotPvSecondT$log2FoldChange

countframe$NPvTT.padj=NotPvThirdT$padj
countframe$NPvTT.FC=NotPvThirdT$log2FoldChange

countframe$NPvM.padj=NotPvMale$padj
countframe$NPvM.FC=NotPvMale$log2FoldChange

countframe$FTvST.padj=FirstTvsSecondT$padj
countframe$FTvST.FC=FirstTvsSecondT$log2FoldChange

countframe$FTvTT.padj=FirstTvsThirdT$padj
countframe$FTvTT.FC=FirstTvsThirdT$log2FoldChange

countframe$FTvM.padj=FirstTvsMale$padj
countframe$FTvM.FC=FirstTvsMale$log2FoldChange

countframe$STvTT.padj=SecondTvsThirdT$padj
countframe$STvTT.FC=SecondTvsThirdT$log2FoldChange

countframe$STvM.padj=SecondTvsMale$padj
countframe$STvM.FC=SecondTvsMale$log2FoldChange

countframe$TTvM.padj=ThirdTvsMale$padj
countframe$TTvM.FC=ThirdTvsMale$log2FoldChange

dbWriteTable(finalproj,"NotPvFirstT",NotPvFirstT,overwrite=TRUE)
dbWriteTable(finalproj,"NotPvSecondT",NotPvSecondT,overwrite=TRUE)
dbWriteTable(finalproj,"NotPvThirdT",NotPvThirdT,overwrite=TRUE)
dbWriteTable(finalproj,"NotPvMale",NotPvMale,overwrite=TRUE)
dbWriteTable(finalproj,"FirstTvsSecondT",FirstTvsSecondT,overwrite=TRUE)
dbWriteTable(finalproj,"FirstTvsThirdT",FirstTvsThirdT,overwrite=TRUE)
dbWriteTable(finalproj,"FirstTvsMale",FirstTvsMale,overwrite=TRUE)
dbWriteTable(finalproj,"SecondTvsThirdT",SecondTvsThirdT,overwrite=TRUE)
dbWriteTable(finalproj,"SecondTvsMale",SecondTvsMale,overwrite=TRUE)
dbWriteTable(finalproj,"ThirdTvsMale",ThirdTvsMale,overwrite=TRUE)

# Define UI for application that draws a heatmap
ui <- fluidPage(
  titlePanel("Pregnancy associated genes of a live bearing cockroach"),
  sidebarLayout(
    sidebarPanel(
      img(src = "cockroach.png",height=120.3,width=254.65 ),
      selectInput("comp", label = "Groups compared:",
                  choices = c("Not Pregnant vs First Trimester"="NPvFT", "Not Pregnant vs Second Trimester"="NPvST", "Not Pregnant vs Third Trimester"="NPvTT", "Not Pregnant vs Male"="NPvM", "First Trimester vs Second Trimester"="FTvST","First Trimester vs Third Trimester"="FTvTT","First Trimester vs Male"="FTvM","SecondTrimester vs Third Trimester"="STvTT","Second Trimester vs Male"="STvM","Third Trimester vs Male"="TTvM"), selected = "Not Pregnant vs First Trimester"),
      
      sliderInput("pval", label = "Maximum P-value",
                  min = 0, max = 1, value = .05, step = 0.001),
      sliderInput("minfc",label="Minimum Log2 Fold Change from reference:", min=-1000,max=1000,value=1,step=5),
      sliderInput("maxfc",label="Maximum Log2 Fold Change from reference:", min=-1000,max=1000,value=10,step=5)
      
    ),
    # Show a plot of the generated distribution
    mainPanel(
      
      textOutput("rowcnt1"),
     plotlyOutput("selecteddata")
      
      #dataTableOutput("table"),
      #tableOutput("colnames")
      
    )
    
  )
)




# Define server logic required to draw a histogram
server <- function(input, output, session) {
  raw<-countframe
  output$selecteddata=renderPlotly({
    if(input$comp=="NPvFT") {
      filt=quote(NPvFT.padj<input$pval & NPvFT.FC<=input$maxfc & NPvFT.FC>=input$minfc)
      
    }
    else if (input$comp=="NPvST"){
      filt=quote(NPvST.padj<=input$pval & NPvST.FC<=input$maxfc & NPvST.FC>=input$minfc)
    }
    else if (input$comp=="NPvTT"){
      filt=quote(NPvTT.padj<=input$pval & NPvTT.FC<=input$maxfc & NPvTT.FC>=input$minfc)
    }
    else if (input$comp=="NPvM"){
      filt=quote(NPvM.padj<=input$pval & NPvM.FC<=input$maxfc & NPvM.FC>=input$minfc)
    }
    else if (input$comp=="FTvST"){
      filt=quote(FTvST.padj<=input$pval & FTvST.FC<=input$maxfc & FTvST.FC>=input$minfc)
    }
    else if (input$comp=="FTvTT"){
      filt=quote(FTvTT.padj<=input$pval & FTvTT.FC<=input$maxfc & FTvTT.FC>=input$minfc)
    }
    else if (input$comp=="FTvM"){
      filt=quote(FTvM.padj<=input$pval & FTvM.FC<=input$maxfc & FTvM.FC>=input$minfc)
    }
    else if (input$comp=="STvTT"){
      filt=quote(STvTT.padj<=input$pval & STvTT.FC<=input$maxfc & STvTT.FC>=input$minfc)
    }
    else if (input$comp=="STvM"){
      filt=quote(STvM.padj<=input$pval & STvM.FC<=input$maxfc & STvM.FC>=input$minfc)
    }
    else {
      filt=quote(TTvM.padj<=input$pval & TTvM.FC<=input$maxfc & TTvM.FC>=input$minfc)
    }
    
    dropit=function(x)(x[,2:16])
    dropit1=function(x)(x[,1])
    collabs<-c("Not.Pregnant1","Not.Pregnant2","Not.Pregnant3",
                                        "First.Trimester1","First.Trimester2","First.Trimester3",
                                        "Second.Trimester1","Second.Trimester2","Second.Trimester3",
                                        "Third.Trimester1","Third.Trimester2","Third.Trimester3",
                                        "Male1","Male2","Male3")
  
    m <- list(
      l = 110,
      t=150
    )  
  names1=as.list(dropit1(raw %>% filter_(filt)))
   t=as.matrix(dropit(raw %>% filter_(filt)))
   rownames(t)=names1
   colnames(t)=collabs
   plot_ly(x=colnames(t),y=rownames(t),z=t,type="heatmap",autosize=F,width=700,height=3000)%>%
     layout(margin = m,xaxis=list(side="top",title="Pregnancy Stage"))
      
    
  })
  
  
  rowcnt=reactive({
    if(input$comp=="NPvFT") {
      filt=quote(NPvFT.padj<input$pval & NPvFT.FC<=input$maxfc & NPvFT.FC>=input$minfc)
      
    }
    else if (input$comp=="NPvST"){
      filt=quote(NPvST.padj<=input$pval & NPvST.FC<=input$maxfc & NPvST.FC>=input$minfc)
    }
    else if (input$comp=="NPvTT"){
      filt=quote(NPvTT.padj<=input$pval & NPvTT.FC<=input$maxfc & NPvTT.FC>=input$minfc)
    }
    else if (input$comp=="NPvM"){
      filt=quote(NPvM.padj<=input$pval & NPvM.FC<=input$maxfc & NPvM.FC>=input$minfc)
    }
    else if (input$comp=="FTvST"){
      filt=quote(FTvST.padj<=input$pval & FTvST.FC<=input$maxfc & FTvST.FC>=input$minfc)
    }
    else if (input$comp=="FTvTT"){
      filt=quote(FTvTT.padj<=input$pval & FTvTT.FC<=input$maxfc & FTvTT.FC>=input$minfc)
    }
    else if (input$comp=="FTvM"){
      filt=quote(FTvM.padj<=input$pval & FTvM.FC<=input$maxfc & FTvM.FC>=input$minfc)
    }
    else if (input$comp=="STvTT"){
      filt=quote(STvTT.padj<=input$pval & STvTT.FC<=input$maxfc & STvTT.FC>=input$minfc)
    }
    else if (input$comp=="STvM"){
      filt=quote(STvM.padj<=input$pval & STvM.FC<=input$maxfc & STvM.FC>=input$minfc)
    }
    else {
      filt=quote(TTvM.padj<=input$pval & TTvM.FC<=input$maxfc & TTvM.FC>=input$minfc)
    }
    
    dropit=function(x)(x[,2:16])
    addnames=function(x)(colnames(x)<-c("Not.Pregnant1","Not.Pregnant2","Not.Pregnant3",
                                        "First.Trimester1","First.Trimester2","First.Trimester3",
                                        "Second.Trimester1","Second.Trimester2","Second.Trimester3",
                                        "Third.Trimester1","Third.Trimester2","Third.Trimester3",
                                        "Male1","Male2","Male3"))
    nrow(dropit(raw %>%
                  filter_(filt)))
    
    
  })
  
  contignames=reactive({
    if(input$comp=="NPvFT") {
      filt=quote(NPvFT.padj<input$pval & NPvFT.FC<=input$maxfc & NPvFT.FC>=input$minfc)
      
    }
    else if (input$comp=="NPvST"){
      filt=quote(NPvST.padj<=input$pval & NPvST.FC<=input$maxfc & NPvST.FC>=input$minfc)
    }
    else if (input$comp=="NPvTT"){
      filt=quote(NPvTT.padj<=input$pval & NPvTT.FC<=input$maxfc & NPvTT.FC>=input$minfc)
    }
    else if (input$comp=="NPvM"){
      filt=quote(NPvM.padj<=input$pval & NPvM.FC<=input$maxfc & NPvM.FC>=input$minfc)
    }
    else if (input$comp=="FTvST"){
      filt=quote(FTvST.padj<=input$pval & FTvST.FC<=input$maxfc & FTvST.FC>=input$minfc)
    }
    else if (input$comp=="FTvTT"){
      filt=quote(FTvTT.padj<=input$pval & FTvTT.FC<=input$maxfc & FTvTT.FC>=input$minfc)
    }
    else if (input$comp=="FTvM"){
      filt=quote(FTvM.padj<=input$pval & FTvM.FC<=input$maxfc & FTvM.FC>=input$minfc)
    }
    else if (input$comp=="STvTT"){
      filt=quote(STvTT.padj<=input$pval & STvTT.FC<=input$maxfc & STvTT.FC>=input$minfc)
    }
    else if (input$comp=="STvM"){
      filt=quote(STvM.padj<=input$pval & STvM.FC<=input$maxfc & STvM.FC>=input$minfc)
    }
    else {
      filt=quote(TTvM.padj<=input$pval & TTvM.FC<=input$maxfc & TTvM.FC>=input$minfc)
    }
    
    dropit1=function(x)(x[,1])
    addnames=function(x)(colnames(x)<-c("Not.Pregnant1","Not.Pregnant2","Not.Pregnant3",
                                        "First.Trimester1","First.Trimester2","First.Trimester3",
                                        "Second.Trimester1","Second.Trimester2","Second.Trimester3",
                                        "Third.Trimester1","Third.Trimester2","Third.Trimester3",
                                        "Male1","Male2","Male3"))
    dropit1(raw %>%
                  filter_(filt))
    
    
  })
  
  output$contigs=renderTable(contignames())
  
  output$rowcnt1=renderText(paste(rowcnt(), "genes match your criteria"))


  table1=renderTable(selecteddata(),rownames=TRUE)
  
 
output$heatmap=renderPlot({heatmaply(selecteddata)})
  
}
# Run the application 
shinyApp(ui = ui, server = server)

