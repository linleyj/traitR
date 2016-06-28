# server.R
source("traitR_pcoa.R")



histfun <- function(bsobj, pvalobj, dataobj, coln) {
  hist(bsobj[, coln], xlab = paste(names(bsobj[coln])), main = paste(names(pvalobj[coln])))
  abline(v = dataobj[, coln], col = "red")
}
shinyServer(function(input, output, session) {
  dsnames <- c()
  
  data_set <- reactive({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    data_set <- read.csv(
      inFile$datapath,
      header = input$header,
      sep = input$sep,
      quote = input$quote,
      stringsAsFactors = T,
      na.strings = c("NA", "N/A", "", " ")
    )
  })
  
  observe({
    dsnames <- names(data_set())
    cb_options <- list()
    cb_options[dsnames] <- dsnames
    updateSelectInput(
      session,
      "speciesCol",
      label = "Species Column",
      choices = cb_options,
      selected = ""
    )
    updateSelectInput(
      session,
      "attackedCol",
      label = "Attacked Column",
      choices = cb_options,
      selected = ""
    )
    updateSelectInput(
      session,
      "subsetCol",
      label = "Column for Subsetting",
      choices = cb_options,
      selected = ""
    )
    updateSelectInput(
      session,
      "traitGrp",
      label = "Traits for PCA",
      choices = cb_options,
      selected = ""
    )
    updateNumericInput(session, "niter",
                       label = "Number of permutations:")
  })
  
  
  
  
  observe({
    df <- data_set()
    if (!is.null(input$speciesCol)) {
      anc_names <- list()
      anc_names <- levels(df[[input$speciesCol]])
      anc_options <- list()
      anc_options[anc_names] <- anc_names
      updateSelectInput(
        session,
        "ancestralGrp",
        label = "Ancestral Species",
        choices = anc_options,
        selected = anc_options[1]
      )
    }
  })
  
  observe({
    df <- data_set()
    if (!is.null(input$attackedCol)) {
      att_names <- list()
      att_names <- levels(df[[input$attackedCol]])
      att_options <- list()
      att_options[att_names] <- att_names
      updateSelectInput(
        session,
        "attackedGrp",
        label = "Attacked Identifier",
        choices = att_options,
        selected = att_options[1]
      )
    }
  })
  
  observe({
    df <- data_set()
    if (!is.null(input$subsetCol)) {
      s_names <- list()
      ss <- as.factor(df[[input$subsetCol]])
      s_names <- levels(ss)
      s_options <- list()
      s_options[s_names] <- s_names
      updateSelectInput(
        session,
        "subsetGrp",
        label = "Subset Identifier",
        choices = s_options,
        selected = s_options[2]
      )
    }
  })
  
  observe({
    df <- data_set()
    if (!is.null(input$traitGrp)) {
      traitgrp <- input$traitGrp
      facgrp<-input$traitFac
      ord <-traitgrp[!traitgrp %in% facgrp]
      updateSelectInput(
        session,
        "traitOrdFac",
        label = "Ordered Factor Traits",
        choices = ord,
        selected = ""
      )}})
  
  observe({
    df <- data_set()
    if (!is.null(input$traitGrp)) {
      traitgrp <- input$traitGrp
     updateSelectInput(
        session,
        "traitFac",
        label = "Factor Traits",
        choices = traitgrp,
        selected = ""
      )
  }})
  
  multanalysis <- reactive({
    df <- data_set()
    gp <- NULL
    if (!is.null(df)) {
      xv <- input$speciesCol
      yv <- input$attackedCol
      zv <- input$traitGrp
      sc <- input$subsetCol
      if (sc == "")
        sc <- NULL
      sg <- input$subsetGrp
      tf <- input$traitFac
      tof <- input$traitOrdFac
      if (!is.null(xv) & !is.null(yv) & !is.null(zv)) {
        if (sum(zv %in% names(df)) > 1) {
          df1 <<- traitR.df(df, xv, yv, zv, sc, sg, tf, tof)
          pscs <- data.frame(pcscores(df1, xv, yv, tf, tof))
          return(pscs)
        }
      }
    }
  })
  
  
  output$plot1 = renderPlot({
    df <- data_set()
    gp <- NULL
    if (!is.null(df)) {
      xv <- input$speciesCol
      yv <- input$attackedCol
      zv <- input$traitGrp
      tf <- input$traitFac
      tof <- input$traitOrdFac
      if (!is.null(xv) & !is.null(yv) & !is.null(zv)) {
        if (sum(zv %in% names(df)) > 1) {
          if (!is.null(tf) | !is.null(tof))
          {
            xl <- "PCOA Component 1"
            yl <- "PCOA Component 2"
          }
          else{
            xl <- "PCA Component 1"
            yl <- "PCA Component 2"
            
          }
          
          plotdata <- data.frame(multanalysis())
          gp <-
            ggplot(data = plotdata, aes(x = Comp.1, y = Comp.2)) + geom_point(aes(col =
                                                                                    Species.Attacked)) + xlab(xl) + ylab(yl)
        }
      }
      return(gp)
    }
  })
  
  observeEvent(input$do,
               {
                 traitRvalues <- NULL
                 df <- data_set()
                 sc <- input$speciesCol
                 ac <- input$attackedCol
                 as <- input$ancestralGrp
                 ag <- input$attackedGrp
                 tg <- input$traitGrp
                 ni <- input$niter
                 if (!is.null(sc) &
                     !is.null(ac) &
                     !is.null(as) &
                     !is.null(ag) & !is.null(tg) & !is.null(ni))
                   pcadf <<- data.frame(multanalysis())
                 traitRvalues <-
                   traitR.perm(
                     pcadf = pcadf,
                     speciescol = sc,
                     attackedcol = ac,
                     ancspecies = as,
                     attackedgroup = ag,
                     niter = ni
                   )
                 
                 output$sum <- renderTable({
                   if (!is.null(pcadf))
                     tabSum <-
                       data.frame(table(pcadf$Species.Attacked))
                   return(tabSum)
                 })
                 
                 output$pvalueLoc <- renderTable({
                   if (!is.null(traitRvalues))
                     tabL <- data.frame(traitRvalues[[3]][1, 1:3])
                   colnames(tabL) <- "P value"
                   return(tabL)
                 })
                 output$pvalueSpread <- renderTable({
                   if (!is.null(traitRvalues))
                     tabS <- data.frame(traitRvalues[[3]][1, 4:5])
                   colnames(tabS) <- "P value"
                   return(tabS)
                 })
                 
                 output$speciescentroids <-
                   renderText("Species Centroids")
                 output$centroidtable <- renderTable({
                   if (!is.null(traitRvalues))
                     tab2 <- data.frame(traitRvalues[[4]])
                   return(tab2)
                 })
                 
                 output$attackedcentroids <-
                   renderText("Attacked Centroids")
                 output$centroidtable2 <- renderTable({
                   if (!is.null(traitRvalues))
                     tab3 <- data.frame(traitRvalues[[5]])
                   return(tab3[c(1, 3),])
                 })
                 
                 
                 output$plot2 <- renderPlot({
                   if (!is.null(traitRvalues))
                     bstab <- data.frame(traitRvalues[[2]])
                   pvalobj <- data.frame(traitRvalues[[3]])
                   datatab <- data.frame(traitRvalues[[1]])
                   
                   p <- function()
                   {
                     par(mfrow = c(2, 3))
                     histfun(bstab, pvalobj, datatab, 1)
                     histfun(bstab, pvalobj, datatab, 2)
                     histfun(bstab, pvalobj, datatab, 3)
                     histfun(bstab, pvalobj, datatab, 4)
                     histfun(bstab, pvalobj, datatab, 5)
                   }
                   return(p())
                 })
                 traitR.out <<- traitRvalues
                 dtableInput <- reactive({
                   if (is.null(traitRvalues))
                     return()
                   datas <- data.frame(traitRvalues[[2]])
                   #row.names(datas)<- c("Calculated Statistic", "Pvalue")
                   return(datas)
                 })
                 
                 # downloadHandler() takes two arguments, both functions.
                 # The content function is passed a filename as an argument, and
                 #   it should write out data to that filename.
                 output$downloadData <- downloadHandler(
                   # This function returns a string which tells the client
                   # browser what name to use when saving the file.
                   filename = function() {
                     ns <- levels(df[[sc]])[levels(df[[sc]]) != as]
                     paste("~\\", as, "_" , ns, ".", "csv", sep = "")
                   },
                   
                   # This function should write data to a file given to it by
                   # the argument 'file'.
                   content = function(file) {
                     # Write to a file specified by the 'file' argument
                     write.csv(dtableInput(), file, row.names = F)
                     return(dtableInput)
                   }
                 )
               })
  
  
  output$filetable <- renderTable({
    if (is.null(input$file1)) {
      return(NULL)
    }  else {
      df <- data_set()
      xv <- input$speciesCol
      yv <- input$attackedCol
      zv <- input$traitGrp
      if (!is.null(xv) & !is.null(yv) & !is.null(zv)) {
        pcadf <- data.frame(multanalysis())
        ht <- rbind(head(pcadf), tail(pcadf))
        return(ht)
      }
    }
  })
})