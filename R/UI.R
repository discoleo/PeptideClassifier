
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
		tabPanel("Clustering", # icon = icon("Clustering"),
			mainPanel(DT::DTOutput("tblClusters"))
		),
	) ) )
}