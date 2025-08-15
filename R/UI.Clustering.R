#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## UI: Clustering
##
## URL: https://github.com/discoleo/PeptideClassifier
##
## draft v.0.2c


### Tab: Clustering
panelClustering = function() {
	sidebarLayout(
	sidebarPanel(
	h3("Subtree:"),
	fluidRow(
	column(4, textInput("txtSubTree_Node", "Node/Peptide", "")),
	column(4, textInput("txtSubTree_Size", "Size", "20")),
	),
	actionButton("btnSubtree", "SubTree"),
	actionButton("btnSubT_Details", "Peptides"),
	# SubTree Data:
	downloadButton("downloadSubT_Data", "Data"),
	fluidRow(HTML("&nbsp;")),
	fluidRow(
	column(6, selectClusterType()),
	column(6,)
	),
	),
	mainPanel(
		plotOutput("imgTree"),
		DT::DTOutput("tblClusters"),
		plotOutput("imgSubTree"),
		DT::DTOutput("tblSubTree"),
	)
	)
}

selectClusterType = function(selected = "complete", id = "fltTreeType") {
	selectInput(id, label = "Clustering method:",
		choices = list(
			"Single"  = "single",  "Complete" = "complete",
			"Ward D"  = "ward.D",  "Ward D2"  = "ward.D2",
			"Average" = "average", "McQuitty/WPGMA" = "mcquitty",
			"Median/WPGMC" = "median", "Centroid/UPGMC" = "centroid"),
		selected = selected);
}
