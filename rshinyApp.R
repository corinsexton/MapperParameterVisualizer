library(shiny)
library(readr)
library(tidyr)

library(ggplot2)
library(data.table)

getFilePrefix <- function(input) {
  ign <- ''
  if (input$data == "simulated" && input$mapper == 'bowtie2') {
    sim = 'sim_'
  } else {
    sim = ''
  }
  
  if (input$mapper == 'bowtie2') {
    if (input$ign_quals == "yes") {
      ign <- 'ignore_quals'
      filename <- paste(input$seedlen,
                        input$seedmismatch,
                        input$seedattempts,
                        input$reseed,
                        ign, sep = '_')
    } else {
      filename <- paste(input$seedlen,
                        input$seedmismatch,
                        input$seedattempts,
                        input$reseed,
                        sep = '_')
    }
  } else {

    seed_attempts <- format(round(input$bwa_seedattempts, 1), nsmall = 1)
    filename <- paste(input$bwa_seedlen,
                      seed_attempts,
                      sep = '_')
  }
  
  filename <- paste(sim, filename, sep = '')
}

# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Alignment Parameter Visualizer"),
  tags$h1(headerPanel("an E. coli case study")),
  
  # Sidebar layout with a input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
           selectInput(inputId = "mapper",
                        label = "Choose a mapper:",
                        choices = c("bowtie2", "bwamem")),
           selectInput(inputId = "data",
                       label = "Choose a dataset:",
                       choices = c("miseq", "simulated")),
           
          conditionalPanel(
             condition = 'input.mapper == "bowtie2"',
             
             sliderInput(inputId = "seedlen",
                         label = "Seed substring length (-L default = 20):",
                         min = 4, max = 31, value = 20),
             
             radioButtons(inputId = "seedmismatch",
                          label = "Max number of seed alignment mismatches (-N default = 0):",
                          choices = list(0,1), inline = T),
             
             sliderInput(inputId = "seedattempts",
                         label = "Number of seed extension attempts allowed to fail (-D default = 15):",
                         min = 5, max = 25, value = 15, step = 5),
             
             radioButtons(inputId = "reseed",
                          label = "Number of time to reseed reads with repetitive
                          seeds (-R default = 2):", 
                          choices = list(1,2,3), selected = 2, inline = T),
             
             radioButtons(inputId = "ign_quals",
                          label = "--ignore-quality toggle (default:no):", 
                          choices = list("no", "yes"), selected = "no", inline = T)
             
             ),
          conditionalPanel(
            condition = 'input.mapper == "bwamem"',
            sliderInput(inputId = "bwa_seedlen",
                        label = "Minimum seed length (-k default = 19):",
                        min = 4, max = 31, value = 19),
            sliderInput(inputId = "bwa_seedattempts",
                        label = "Trigger re-seeding for a MEM longer than minSeedLen*FLOAT  (-r default = 1.5):",
                        min = 0.5, max = 3.0, value = 1.5, step = 0.5)
            ),
          width = 3
        ),
           
    
    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Verbatim text for data summary ----
      
      fluidRow(
        splitLayout(cellWidths = c("48.75%", "48.75%"), 
                    verbatimTextOutput("command"),
                    verbatimTextOutput("accuracy"))
      ),
      

      fluidRow(
        splitLayout(cellWidths = c("33%", "33%", "33%"), 
                    plotOutput("picard", height="300px"),
                    plotOutput("picard1", height = "300px"),
                    plotOutput("picard2", height = "300px"))
      ),
      
      
      fluidRow(
        plotOutput("plot", height="350px")
      ),
      
      fluidRow(
        splitLayout(cellWidths = c("99%"),
                    plotOutput("bargraph", height="350px"))
        
      ),
      width=9
      
    )
  )
)



# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  
  output$command <- renderPrint({
    
    if (input$mapper == "bowtie2") {
      ign <- ''
      if (input$ign_quals == "yes") {
        ign <- ' --ignore-quals'
      }
      
      cat(paste("bowtie2 -L", input$seedlen, 
                "-N", input$seedmismatch, 
                "-D", input$seedattempts, 
                "-R", input$reseed), ign, sep = '')
    } else {
      seed_attempts <- format(round(input$bwa_seedattempts, 1), nsmall = 1)
      
      cat(paste("bwa mem -k", input$bwa_seedlen, 
                "-r", seed_attempts, 
                sep = ' ')
      )
    }
    

  })
  
  output$accuracy <- renderPrint({
    
    cat("Accuracy:  ")
    out <- 'accuracy not available for miseq data'
    
    if (input$data == 'simulated') {
      if (input$mapper == "bowtie2") {
        acc <- read_tsv("~/COLLEGE/winter2018/bio465/project/bowtie2/simulated/accuracy.tsv")
      }  else {
        acc <- read_tsv("~/COLLEGE/winter2018/bio465/project/bwamem/simulated/accuracy.tsv")
      }
      
      pref <- paste('bamfiles/', getFilePrefix(input), '.bam', sep = '')
      y <- subset(acc, file %like% pref)

      out <- format(round(y$accuracy, 4), nsmall = 4)
    }
    cat(out)
    
  })
  
  output$plot <- renderPlot({
    
    filename <- getFilePrefix(input)
    filename1 <- paste(filename, ".depth.10000window", sep = '')
    
    x <- read_tsv(paste("~/COLLEGE/winter2018/bio465/project/", input$mapper, '/',
                        input$data, "/depth/", filename1, sep =''))
    
    
    if (input$data == 'miseq') {
      begin = 200
      end = 470
    } else {
      begin = 150
      end = 215
    }
      filename42 <- paste(filename, "_42", ".depth.10000window", sep = '')
      
      y <- read_tsv(paste("~/COLLEGE/winter2018/bio465/project/", input$mapper, '/',
                          input$data, "/depth/", filename42, sep =''))
      
      z <- merge(x,y, by="pos")
      
      p <- ggplot(z,aes(pos)) +
        geom_line(aes(y = depth.y, colour = "High Qual")) +
        geom_line(aes(y = depth.x, colour = "All Quals")) + 
        labs(x = "Position", y = "Coverage Depth per 10,000bp") +
        ylim(begin,end) +
        theme(text = element_text(size=18)) +
        scale_colour_manual(name = NULL, values = c("#FF7F50", "#00B6EB"))
      print(p)
  })
  
  output$bargraph <- renderPlot({

    filename <- getFilePrefix(input)

    filename <- paste(filename, ".mapq", sep = '')

    x <- read_tsv(paste("~/COLLEGE/winter2018/bio465/project/", input$mapper, '/',
                       input$data, "/mapq/", filename, sep =''))
    #x <- x[-1]
    #n <- 10000;
    #y <- aggregate(x,list(rep(1:(nrow(x)%/%n+1),each=n,len=nrow(x))),mean)[-1]

    p <-ggplot(x,aes(x=as.factor(mapq), y = log(count))) +
      geom_bar(stat='identity', fill = 'lightblue') +
      labs(x = "Mapping Quality", y = "log(Number of Reads)") +
      theme(text = element_text(size=18), axis.text.x = element_text(angle = 60, hjust = 1, size = 10))

    print(p)

  })
  
  output$picard <- renderPlot({

    if (input$data == 'miseq') {
      begin = 0
      end = 0.016
    } else {
      begin = 0
      end = 0.005
    }
    
     filename <- getFilePrefix(input)
     
     filename <- paste(filename, ".picard_stats", sep = '')
     x <- read_tsv(paste("~/COLLEGE/winter2018/bio465/project/", input$mapper, '/',
                         input$data, "/picard_stats/", filename, sep =''), comment = '#', skip = 1)
     
     tidied <- gather(x, key = type, value = books, PF_MISMATCH_RATE, PF_INDEL_RATE, -CATEGORY)
     
     p <- ggplot(subset(tidied, CATEGORY != "PAIR"), aes(x=CATEGORY, y = books, fill=type)) +
       geom_bar(stat='identity', position = "dodge") +
       labs(x = NULL, y = NULL) + 
       theme(text=element_text(size= 14), legend.position="bottom", legend.direction="vertical") +
       scale_x_discrete(labels = c("first read in pair", "second read in pair")) +
       scale_fill_discrete(name = NULL, 
                           labels = c("Indel rate per 100 aligned bases", 
                                      "Mismatch rate of all bases")) + 
       ylim(begin,end)

     print(p)
  })
  
  output$picard1 <- renderPlot({
    if (input$data == 'miseq') {
      begin = 0
      end = 0.035
    } else {
      begin = 0
      end = 0.005
    }
    
    filename <- getFilePrefix(input)
    
    filename <- paste(filename, ".picard_stats", sep = '')
    x <- read_tsv(paste("~/COLLEGE/winter2018/bio465/project/", input$mapper, '/',
                        input$data, "/picard_stats/", filename, sep =''), comment = '#', skip = 1)
    tidied <- gather(x, key = type, value = books, PCT_PF_READS_IMPROPER_PAIRS, -CATEGORY)
    
    ggplot(subset(tidied, CATEGORY != "PAIR"),aes(x=CATEGORY, y = books)) +
      geom_bar(stat='identity', fill = "slateblue1") +
      labs(x = NULL, y = "Percent of Reads with Imperfect Pair") + 
      theme(text=element_text(size= 14), legend.position="bottom", legend.direction="vertical") +
      scale_x_discrete(labels = c("first read in pair", "second read in pair")) +
      ylim(begin,end)
    
  })
  
  output$picard2 <- renderPlot({
    if (input$data == 'miseq') {
      begin = 4000000
      end = 6000000
    } else {
      begin = 1000000
      end = 3500000
    }
    filename <- getFilePrefix(input)
    
    filename <- paste(filename, ".picard_stats", sep = '')
    x <- read_tsv(paste("~/COLLEGE/winter2018/bio465/project/", input$mapper, '/',
                        input$data, "/picard_stats/", filename, sep =''), comment = '#', skip = 1)
    tidied <- gather(x, key = type, value = books, TOTAL_READS, PF_READS, PF_READS_ALIGNED,
                     PF_HQ_ALIGNED_READS, -CATEGORY)
    tidied$type <- factor(tidied$type, levels=c("TOTAL_READS", 
                                                "PF_READS", 
                                                "PF_READS_ALIGNED", 
                                                "PF_HQ_ALIGNED_READS"))
    
    
    p <- ggplot(subset(tidied, CATEGORY != "PAIR"),aes(x=CATEGORY, y = books, fill=type)) +
      geom_bar(stat='identity', position = "dodge") +
      labs(x = NULL, y = NULL) + 
      theme(text=element_text(size= 14), legend.position="bottom", legend.direction="vertical") +
      scale_x_discrete(labels = c("first read in pair", "second read in pair")) +
      scale_fill_discrete(name = NULL, labels = c("Total reads", 
                                                  "Reads passed Illumina filter (PF)", 
                                                  "PF reads aligned", 
                                                  "PF reads aligned with >= Q20")
      ) + coord_cartesian(ylim=c(begin, end))
    
    print(p)
    
    
  })
  
}


shinyApp(ui, server)
  