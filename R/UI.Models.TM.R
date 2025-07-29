
### Tab: Topic Models

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
