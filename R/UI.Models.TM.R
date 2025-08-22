
### Tab: Topic Models

### Topic Models
panelTopicModel = function() {
	sidebarLayout(
		sidebarPanel(
			sliderInput(inputId = "numClusters", label = "Clusters / Topics",
				value = 6, min = 2, max = 42, step = 1),
			selectModel.TM(),
			actionButton("btnModelTopics", "Run Model"),
			downloadButton("downloadTMSummary", "Summary"), # Export Summary
			downloadButton("downloadTMTopics", "Topics"),   # Export Topics
			downloadButton("downloadTMFull", "Full"),       # Export Full Model(s)
			# Explore Topics:
			fluidRow(tag("h3", "Visualise Topics:")),
			# actionButton("btnVisualise", "Visualise"),
			fluidRow(
			column(6, textInput("inModelID", "Model ID", "1", width = 150)),
			column(6, textInput("inTopicID", "Topic ID", "1", width = 150)),
			),
			sliderInput(inputId = "numTermsTM", label = "Display Terms",
				value = 10, min = 5, max = 45, step = 5),
			# Existing Models:
			fileInput.rds("loadTMFull", "Load existing Models"),
			# LDA-Control:
			fluidRow(tag("h3", "Model Control:")),
			fluidRow(
			column(6, textInput("inTMSeed", "Seed", "0", width = 150)),
			column(6, textInput("inTMDownloadTop", "Download Topics", "1", width = 150)),
			),
			sliderInput(inputId = "numTM_Iterations", label = "Max Iterations",
				value = 1000, min = 500, max = 5000, step = 100),
		),
		### Main Panel
		mainPanel(
			# Topic Models:
			fluidRow(tag("h1", "Topics")),
			fluidRow(DT::DTOutput("tblTopicInfo")),
			fluidRow("Note: may take some time to compute the model."),
			# PP in Cluster[i]:
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
