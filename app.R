library(shiny)
library(rintrojs)
library(DT)

hlas <- readRDS('data/allelelist.RDS')
hlas <- hlas[1:87,1]
HLAfiles <- list.files(path="HLA", pattern = "*.RDS")

gtex_annot <- readRDS("data/gtex_annot.RDS")
genekey <- readRDS("data/genekey.RDS")
gtex_p0 <- readRDS("data/gtex_p0.RDS")
gtex_p1 <- readRDS("data/gtex_p1.RDS")
gtex_p2 <- readRDS("data/gtex_p2.RDS")
gtex_p3 <- readRDS("data/gtex_p3.RDS")
gtex_p4 <- readRDS("data/gtex_p4.RDS")
gtex_p5 <- readRDS("data/gtex_p5.RDS")
gtex_p6 <- readRDS("data/gtex_p6.RDS")
gtex_p7 <- readRDS("data/gtex_p7.RDS")
gtex_p8 <- readRDS("data/gtex_p8.RDS")
gtex_p9 <- readRDS("data/gtex_p9.RDS")

gtex <- rbind(gtex_p0,gtex_p1,gtex_p2,gtex_p3,gtex_p4,gtex_p5,gtex_p6,gtex_p7,gtex_p8,gtex_p9)

combined <- readRDS("data/combined_normal_ligands.RDS")



  ui <- shinyUI(fluidPage(

 title = "sCRAP",
 introjsUI(),
 

 
 introBox(titlePanel(img(src='sCRAP.png', width="50%")), 

  data.step = 1,
  data.intro = "sCRAP is a tool for preemptively predicting potential cross-reactivities for peptides presented on HLA. Predictions are based on the peptide sequence alone and do not require the prior development of a receptor."
 ),

actionButton("Shiny_CRAP_Help", "Tutorial", icon = icon("question"),
                  class = "btn-xs", title = "Help"),



  br(),

   
  sidebarLayout(


    
    sidebarPanel(h4("Query Input", width=3, fluid=F),
  	  introBox(
               	introBox(
                  selectInput("HLA", h5("HLA:"), choices = unique(gsub(".p1|.p2|.p3|.p4|.p5|.RDS", "", HLAfiles))),
 		  data.step = 3,
                  data.intro = "HLA: select HLA allele on which target peptide is presented"
	       	),
		introBox(
                  textInput("PepIn", h5("Peptide Sequence:", value="KVAELVHFL", placeholder="KVAELVHFL")),
 		  data.step = 4,
                  data.intro = "Peptide: Enter target peptide sequence"
		),
		introBox(
                  checkboxGroupInput("hotspots", h5("Hot-spots:"),choices= list("3"=3, "4"=4, "5"=5, "6"=6, "7"=7, "8"=8), selected =list(4,5,6), inline=T),
 		  data.step = 5,
                  data.intro = "Hot-spots: if key interaction residues are known in the target peptide (hot-spots), select the amino acid positions of these residues."
		),
		introBox(
                  actionButton("process", "Process"),
 		  data.step = 6,
                  data.intro = "Press the Process button to generate a cross reactive peptide table."
		),
 		  data.step = 2,
                  data.intro = "Enter information about the HLA allele and peptide in the query panel."

		)
 
  ),               
    
    mainPanel(


              fluidRow(introBox(h4("Cross Reactive Peptide Output"),
		  data.step = 7,
                  data.intro = "Output table with possible cross-reacive peptides will be generated here. For screening cross-reactivity, we recommend selecting top ranking peptides based on both peptides score an overall score. Additionally, we recommend prioritizing peptides that have been detected in the normal ligandome."

		)),
		
              	  h5(textOutput("txt")),
               	  h5(textOutput("pep")),
              	  h5(textOutput("status")),
                  introBox(
		    downloadButton('download',"Download Results"),
		    data.step = 8,
		    data.intro = "An output table with possible cross-reacive peptides will be generated below after clicking the 'process' button. Click this download button to download the data in CSV format."
		  ),

		fluidRow(
		  tags$img(
                   src = "CrossReact.png",
                   style = 'position: absolute; padding:0px'
                  ),		 
              	  DT::dataTableOutput('tab1')

		)
  
))))


  
  

server <- function(input, output, session,Output_peptable = c("")){
  
  output$txt <- renderText({
    if (input$PepIn != ""){
    pepinput <- paste(input$PepIn)
    paste("Peptide", pepinput, "on", input$HLA)
  }
  })
  
  
  observeEvent(input$process, {
    withProgress(message = 'Computing: this will take several minutes', value = 0, {

  
  G=list()
  G[[1]] = c("A","G") #Short Chains
  G[[2]] = c("K","R","H") #Basic
  G[[3]] = c("N","Q") #Polar Uncharged
  G[[4]] = c("D","E") #Acidic
  G[[5]] = c("C", "M") #Contains S
  G[[6]] = c("F","Y","W","H") #Aromatic
  G[[7]] = c("S","T","Y") #Alcohol-Hydroxyl Group
  G[[8]] = c("I","L","V","M","A") #Aliphatic
  G[[9]] = c("P") #Weird - I mean proline
  
  Polarity = list()
  Polarity[[1]] = c("D","E","N","Q","R","K","H","Y","C","S","T") #Polar
  Polarity[[2]] = c("G","A","F","W","P","I","L","V","M") #Non Polar
  
  Charge = list()
  Charge[[1]] = c("D","E") #Negative
  Charge[[2]] = c("K","R","H") #Positive
  
 
  
   HLA_RDS_Files <- system(paste("ls HLA/", input$HLA, ".p*.RDS", sep=""), intern=T)

   combinedpeps <- data.frame()
   for(I in 1:length(HLA_RDS_Files))
   {
	combinedpeps <- rbind(combinedpeps,data.frame(readRDS(HLA_RDS_Files[I])))
   }
 
  
  
  peptable <- combinedpeps[,c(3,11,13,15)]
  peptable$Identity <- as.character(gsub("_HUMAN","", as.character(peptable$Identity)))
  
 
 

  #Function to Compare Peptides
  Get_Peptide_Scores = function(QUERY,PEPTIDE)
  {
    QUERY.array = unlist(strsplit(as.character(QUERY),split=""))
    PEPTIDE.array = unlist(strsplit(as.character(PEPTIDE),split=""))
    
    Peptide_Score = -100
    if(
      length(QUERY.array) == length(PEPTIDE.array) #&
    )#if
    {#Start function
      TEST.Pair=rbind(QUERY.array,PEPTIDE.array)
      
      #Equivalent
      EQU_score = apply(TEST.Pair,MARGIN = 2,function(x) {Score=0;Criteria = x[1]==x[2]; if(Criteria){Score=3};return(Score)})
      #GROUP
      GROUP_score = apply(TEST.Pair,MARGIN=2,function(X) {Score=0;Criteria=any(unlist(lapply(G,function(G) all(X %in% G))));if(Criteria){Score=2};return(Score)} )
      #Polarity
      Polarity_score = apply(TEST.Pair,MARGIN=2,function(X) {Score=0;Criteria=any(unlist(lapply(Polarity,function(G) all(X %in% G))));if(!Criteria){Score=-2};return(Score)})
      #Charge
      Charge_score = apply(TEST.Pair,MARGIN=2,function(X) {Score=0;Criteria=all(unlist(lapply(Charge, function(x) any(x %in% G))));if(!Criteria){Score=-1};return(Score)})
      
      #Double Scores for Selected Hotspots.
      
      HOTSPOTS <- as.numeric(unlist(input$hotspots))

      EQU_score[HOTSPOTS] = EQU_score[HOTSPOTS]*2 
      GROUP_score[HOTSPOTS] = GROUP_score[HOTSPOTS]*2 
      Polarity_score[HOTSPOTS] = Polarity_score[HOTSPOTS]*2 
      Charge_score[HOTSPOTS] = Charge_score[HOTSPOTS]*2 
      
      #Sets scores at Positions 2 and last position to 0.
      EQU_score[c(2,length(EQU_score))] = 0
      GROUP_score[c(2,length(GROUP_score))] = 0
      Charge_score[c(2,length(Charge_score))] = 0
      
      #Final Match Score
      Peptide_Score = sum(c(EQU_score,GROUP_score,Polarity_score))
      
    }#{Start function
    return(Peptide_Score)
  }
  
  #normal tissue RPKM
  normal <- function(peptide){
    for (i in 1:nrow(peptide)){
      if (!is.na(g <- match(peptide$Gene[i], genekey$name))){
        g <- match(peptide$Gene[i], genekey$name)
        peptide$Gene[i] <- genekey$gene[g]
      }
      if (!is.na(match(peptide$Gene[i], rownames(gtex)))){
        r <- match(peptide$Gene[i], rownames(gtex))
        peptide$max_norm[i] <- max(gtex[r,])

        tiss <- match(colnames(gtex)[which.max(gtex[r,])[[1]]], gtex_annot$SAMPID)
        peptide$max_tissue[i] <- as.character(gtex_annot$SMTS[as.numeric(tiss)])
      }
      else {
        peptide$max_norm[i] <- NA
        peptide$max_tissue[i] <- NA
      }
      incProgress(0.1/nrow(peptide), message = paste("Scoring normal RPKM ", i, " of ", nrow(peptide),": this will take several minutes"))
      
    }
    return(peptide)
  }
  
  # compare to  normal ligandome
  ligandome <- function(peptide){
    
      
      peptide$Normal_Ligandome <- "F"
      peptide$MHC_class <- ""
      peptide$HLA <- ""
      peptide$Organism <- ""
      for (i in 1:nrow(peptide)){
        if (!is.na(match(peptide$Peptide[i], combined$search_hit))){
          r <- match(peptide$Peptide[i], combined$search_hit)
          peptide$Normal_Ligandome[i] <- "T"
          peptide$MHC_class[i] <- as.character(combined$MHCClass[r])
          peptide$HLA[i] <- as.character(combined$top_allele[r])
          peptide$Organism[i] <- as.character(combined$Organism[r])
        }
        else{
          peptide$Normal_Ligandome[i] <- "F"
        }
        incProgress(0.1/nrow(peptide), message = paste("Scoring normal ligandome ", i, " of ", nrow(peptide),": this will take several minutes"))
        
      }
      
    return(peptide)
    
  }
  
  
  
  incProgress(1/4, message = paste("Scoring peptide: this will take several minutes"))
  
  {

    PEPTIDE_INPUT = input$PepIn
    incProgress(0.1, message = " Filtering table: this will take several minutes")

    peptable.peptide.subset = peptable$peptide[which(sapply(as.character(peptable$peptide),nchar) == nchar(PEPTIDE_INPUT) &  peptable$BindLevel=="<=SB")]

    peptable.peptide.subset = peptable.peptide.subset[which((substr(peptable.peptide.subset, 1,1) == substr(PEPTIDE_INPUT, 1,1)) | (substr(peptable.peptide.subset, 3,3) == substr(PEPTIDE_INPUT, 3,3)) | (substr(peptable.peptide.subset, 4,4) == substr(PEPTIDE_INPUT, 4,4)) | (substr(peptable.peptide.subset, 5,5) == substr(PEPTIDE_INPUT, 5,5)) | (substr(peptable.peptide.subset, 6,6) == substr(PEPTIDE_INPUT, 6,6)) | (substr(peptable.peptide.subset, 7,7) == substr(PEPTIDE_INPUT, 7,7)) | (substr(peptable.peptide.subset, 8,8) == substr(PEPTIDE_INPUT, 8,8)))]

    incProgress(0.15, message = paste("Calculating", length(peptable.peptide.subset)," scores. This may take a few minutes"))

    quarter_count = floor(length(peptable.peptide.subset)/4)
    mod_count = length(peptable.peptide.subset) %% 4

    TEST_Scores.filter_1 = sapply(peptable.peptide.subset[1:quarter_count],function(X) as.numeric(Get_Peptide_Scores(PEPTIDE_INPUT, X)))

    incProgress(0.2, message = paste("Calculating", (length(peptable.peptide.subset)-quarter_count)," scores. This may take a few minutes"))

    TEST_Scores.filter_2 = sapply(peptable.peptide.subset[(quarter_count+1):(quarter_count*2)],function(X) as.numeric(Get_Peptide_Scores(PEPTIDE_INPUT, X)))

    incProgress(0.25, message = paste("Calculating", (length(peptable.peptide.subset)-quarter_count*2)," scores. This may take a few minutes"))
    TEST_Scores.filter_3 = sapply(peptable.peptide.subset[(quarter_count*2 + 1):(quarter_count*3)],function(X) as.numeric(Get_Peptide_Scores(PEPTIDE_INPUT, X)))

    incProgress(0.3, message = paste("Calculating", (length(peptable.peptide.subset)-quarter_count*3)," scores. This may take a few minutes"))
    TEST_Scores.filter_4 = sapply(peptable.peptide.subset[(quarter_count*3 + 1):length(peptable.peptide.subset)],function(X) as.numeric(Get_Peptide_Scores(PEPTIDE_INPUT, X)))
     
    TEST_Scores.filter = c(TEST_Scores.filter_1,TEST_Scores.filter_2,TEST_Scores.filter_3,TEST_Scores.filter_4)

    peptable$score = -100

    incProgress(0.35, message = " Sorting Scores")
    peptable$score[match(peptable.peptide.subset,as.character(peptable$peptide))] = TEST_Scores.filter
    
    peptable$BindLevel[which(peptable$BindLevel == "<=WB")] <- "Weak Binder"
    peptable$BindLevel[which(peptable$BindLevel == "<=SB")] <- "Strong Binder"


    
  }

  Return_TOP = 100

  colnames(peptable) <- c("Peptide", "Gene", "Affinity_(nM)","Bind_Level", "Peptide_Score")
  peptable$Gene <- as.character(peptable$Gene)
  peptable <- peptable[order(peptable$Peptide_Score, decreasing = T)[1: Return_TOP],]
  
   
  peptable <- normal(peptable)
 
 
  peptable <- ligandome(peptable)
  incProgress(.9, message = paste("Almost done"))
  Sys.sleep(0.3)
  
  
  peptable$Overall_Score = round(peptable$Peptide_Score * peptable$max_norm/peptable$'Affinity_(nM)',2)

  Output_peptable = peptable[,c("Peptide","Gene","Peptide_Score","Overall_Score","Affinity_(nM)","max_norm","max_tissue","Normal_Ligandome","HLA")]
 
  Output_peptable$max_tissue <- as.character(Output_peptable$max_tissue) 
  Output_peptable$max_norm <- round(as.numeric(Output_peptable$max_norm),2)

  Output_peptable = Output_peptable[order(Output_peptable$Overall_Score,decreasing=T),]

  rownames(Output_peptable) <- NULL
  Output_peptable <<- Output_peptable
  output$tab1 = DT::renderDataTable({Output_peptable})


  })
  })

  observeEvent(input$Shiny_CRAP_Help,
               introjs(session, options = list("nextLabel"="next",
                                               "prevLabel"="prev",
                                               "skipLabel"="skip",
                                               events = list("oncomplete"=I('alert("Glad that is over")'))))
  )

  output$download <- downloadHandler(
    filename = function(){"Output_peptable.csv"}, 
    content = function(fname){
      write.csv(Output_peptable, fname)
    }
  )


}

#}
shinyApp(ui = ui, server = server)


