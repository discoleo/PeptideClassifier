
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
					sliderInput(inputId = "fltUNK", label = "Not yet",
						value = 0.2, min = 0, max = 5, step = 0.25),
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
		# Topic Models
		tabPanel("Modeling", # icon = icon("Modeling"),
			panelTopicModel()
		),
		tabPanel("Clustering", # icon = icon("Clustering"),
			mainPanel(DT::DTOutput("tblClusters"))
		),
	) ) )
}


panelTopicModel = function() {
	sidebarLayout(
		sidebarPanel(
			selectInput("fltDTMSeq",
				label = "Filter sequences:", selected = "Positive",
				choices = list("All" = "All", "Positive" = "Positive",
					"Negative" = "Negative") ),
			sliderInput(inputId = "fltLen", label = "Length",
				value = c(6, 40), min = 0, max = 100, step = 1),
			sliderInput(inputId = "fltTF", label = "TF-IDF",
				value = 0.1, min = 0, max = 2, step = 0.025),
			sliderInput(inputId = "numClusters", label = "Clusters / Topics",
				value = 6, min = 2, max = 42, step = 1),
			selectInput("fltTMType",
				label = "Model type:", selected = "VEM",
				choices = list("VEM" = "VEM", "fixVEM" = "Fixed VEM",
					"Gibbs" = "Gibbs", "CTM" = "CTM", "All" = "All") ),
			actionButton("btnDTM", "Build DTM"),
			actionButton("btnDTMFilter", "Filter DTM"),
			actionButton("btnModelTopics", "Run Model"),
		),
		mainPanel(
			fluidRow(column(7,
			fluidRow(tag("h1", "Document Term Matrix")),
			fluidRow(DT::DTOutput("tblDTMSummary")),
			)),
			fluidRow(tag("h1", "Topics")),
			column(4,
			fluidRow(DT::DTOutput("tblTopics")),
			fluidRow("Note: may take some time to compute the model.")
			)
		)
	)
}