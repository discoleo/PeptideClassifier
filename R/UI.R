
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
