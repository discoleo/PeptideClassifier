#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## UI: Main
##
## URL: https://github.com/discoleo/PeptideClassifier
##
## draft v.0.2c


fileInput.data = function(id, label, csv = FALSE) {
	accept = c(
			"text/csv",
			"text/comma-separated-values",
			"text/plain", ".csv", ".txt");
	if(! csv) {
		accept = c(accept, "text/fasta", ".fasta", ".fa");
	}
	shiny::fileInput(id, label,
		multiple = FALSE,
		accept = accept
    );
}
fileInput.rds = function(id, label) {
	accept = c(".rds");
	shiny::fileInput(id, label,
		multiple = FALSE,
		accept = accept
    );
}

NBSP = function(as.row = TRUE) {
	html = HTML("&nbsp;");
	if(as.row) {
		html = fluidRow(html);
	}
	return(html);
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
					fileInput.data("file", "Select CSV / FASTA file"),
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
		# Hierarchical Clustering
		tabPanel("Clustering", # icon = icon("Clustering"),
			panelClustering()
		),
		# Clustering: Diagnostics
		tabPanel("Diagnostics", # icon = icon("ClusterDx"),
			panelClusterDiagnostics()
		),
		# Generalised Beta Distribution
		tabPanel("Gen.Beta", # icon = icon("BetaDist"),
			panelBetaDist()
		),
	) ) )
}

