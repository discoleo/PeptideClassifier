
fileInput.csv = function(id, label) {
	shiny::fileInput(id, label,
		multiple = FALSE,
		accept = c(
			"text/csv",
			"text/fasta",
			"text/comma-separated-values,text/plain",
			".csv", ".fasta", ".fa", ".txt")
    );
}


# version = 1: Old variant;
getUI = function(version = 2) {
	### Shiny functions
	# - NO need to be visible in the app;
	# Layout:
	sidebarLayout = function(...) shiny::sidebarLayout(...);
	fluidRow = function(...) shiny::fluidRow(...);
	column = function(...) shiny::column(...);
	# Panels:
	tabPanel = function(...) shiny::tabPanel(...);
	mainPanel = function(...) shiny::mainPanel(...);
	titlePanel = function(...) shiny::titlePanel(...);
	sidebarPanel = function(...) shiny::sidebarPanel(...);
	# Input:
	sliderInput = function(...) shiny::sliderInput(...);
	selectInput = function(...) shiny::selectInput(...);
	checkboxInput = function(...) shiny::checkboxInput(...);
	# <----------------->
	### UI:
	shiny::shinyUI(
	shiny::fluidPage(
	shiny::navbarPage("Peptide Classifier", id="menu.top",
		# theme = init.theme(),
		tabPanel("Data", # icon = icon("upload file"),
			sidebarLayout(
				sidebarPanel(
					titlePanel("Load file"),
					fileInput.csv("file", "Select CSV / FASTA file"),
					# Note: value is initialised in code;
					sliderInput(inputId = "fltGlobalLen", label = "Global Length",
						value = c(6, 60), min = 0, max = 100, step = 1),
					selectInput("fltType",
						label = "Filter Activity:",
						choices = list("All" = "All", "Positive" = "Positive",
							"Negative" = "Negative"),
						selected = "All"),
					checkboxInput("chkRegex", "Regex Search: Data", value = TRUE),
					downloadButton("downloadData", "Download"),
				),
				mainPanel(
				fluidRow(DT::DTOutput("tblData")),
				fluidRow(textOutput("txtDataSummary")),
				fluidRow(
					column(4, DT::DTOutput("tblDataSummary")))
				)
		)),
		# DTM
		tabPanel("DTM", # icon = icon("DTM"),
			panelDTM()
		),
		# Topic Models
		tabPanel("Topics", # icon = icon("Topics"),
			panelTopicModel()
		),
		tabPanel("Clustering", # icon = icon("Clustering"),
			mainPanel(DT::DTOutput("tblClusters"))
		),
	) ) )
}

### DTM
panelDTM = function() {
	sidebarLayout(
		sidebarPanel(
			selectInput("fltDTMSeq",
				label = "Filter sequences:", selected = "Positive",
				choices = list("All" = "All", "Positive" = "Positive",
					"Negative" = "Negative") ),
			sliderInput(inputId = "fltDTMLen", label = "Length",
				value = c(6, 40), min = 0, max = 100, step = 1),
			sliderInput(inputId = "fltTF", label = "TF-IDF",
				value = 0.1, min = 0, max = 2, step = 0.025),
			actionButton("btnDTM", "Build DTM"),
			actionButton("btnDTMFilter", "Filter DTM"),
			# n-Grams:
			fluidRow(HTML("&nbsp;")),
			checkboxGroupInput("chkNGrams", "n-Grams",
				inline = TRUE, width = 360,
				choices = c("2 Ord" = "2", "2 UnOrd" = "2u",
					"3 Ord" = "3", "3 UnOrd" = "3u",
					"4 Ord" = "4", "4 UnOrd" = "4u",
					"PP Length" = "Len",
					"3 Gr-Charge" = "ch3.tot", "3 Gr-Charged AA" = "ch3.aa",
					"4 Gr-Charge" = "ch4.tot", "4 Gr-Charged AA" = "ch4.aa",
					"5 Gr-Charge" = "ch5.tot", "5 Gr-Charged AA" = "ch5.aa"),
				selected = c("2", "2u", "3u", "Len", "ch4.tot", "ch4.aa")
			),
		),
		### Main Panel
		mainPanel(
			fluidRow(column(7,
			fluidRow(tag("h1", "Document Term Matrix")),
			fluidRow(DT::DTOutput("tblDTMSummary")),
			),
			column(5, DT::DTOutput("tblDTMFltTerms"))
			),
		)
	)
}

### Topic Models
panelTopicModel = function() {
	sidebarLayout(
		sidebarPanel(
			sliderInput(inputId = "numClusters", label = "Clusters / Topics",
				value = 6, min = 2, max = 42, step = 1),
			selectInput("fltTMType",
				label = "Model type:", selected = "VEM",
				choices = list("VEM" = "VEM", "fixVEM" = "Fixed VEM",
					"Gibbs" = "Gibbs", "CTM" = "CTM", "All" = "All") ),
			actionButton("btnModelTopics", "Run Model"),
			# Explore Topics:
			fluidRow(tag("h3", "Visualise Topics:")),
			# actionButton("btnVisualise", "Visualise"),
			textInput("inTopicID", "Topic ID", "1"),
		),
		### Main Panel
		mainPanel(
			# Topic Models:
			fluidRow(tag("h1", "Topics")),
			fluidRow(
			column(4,
			fluidRow(DT::DTOutput("tblTopics")),
			fluidRow("Note: may take some time to compute the model.")
			),
			column(8, DT::DTOutput("tblTopicInfo")
			) ),
			fluidRow(DT::DTOutput("tblPPTopic"))
		)
	)
}
