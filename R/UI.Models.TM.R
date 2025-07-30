
### Tab: Topic Models

### Topic Models
panelTopicModel = function() {
	sidebarLayout(
		sidebarPanel(
			sliderInput(inputId = "numClusters", label = "Clusters / Topics",
				value = 6, min = 2, max = 42, step = 1),
			selectModel.TM(),
			actionButton("btnModelTopics", "Run Model"),
			# Explore Topics:
			fluidRow(tag("h3", "Visualise Topics:")),
			# actionButton("btnVisualise", "Visualise"),
			textInput("inTopicID", "Topic ID", "1"),
			sliderInput(inputId = "numTermsTM", label = "Display Terms",
				value = 10, min = 5, max = 45, step = 5),
		),
		### Main Panel
		mainPanel(
			# Topic Models:
			fluidRow(tag("h1", "Topics")),
			fluidRow(
			column(3,
			fluidRow(DT::DTOutput("tblTopics")),
			),
			column(9, DT::DTOutput("tblTopicInfo"),
			fluidRow("Note: may take some time to compute the model.")
			) ),
			fluidRow(DT::DTOutput("tblPPTopic")),
			fluidRow(DT::DTOutput("tblTopicTerms")),
		)
	)
}

selectModel.TM = function(selected = "VEM", id = "fltTMType") {
	selectInput(id,
		label = "Model type:", selected = selected,
		choices = list("VEM" = "VEM", "fixVEM" = "Fixed VEM",
					"Gibbs" = "Gibbs", "CTM" = "CTM", "All" = "All") );
}
