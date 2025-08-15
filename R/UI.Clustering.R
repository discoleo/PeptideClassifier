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
	),
	mainPanel(
		plotOutput("imgTree"),
		DT::DTOutput("tblClusters"),
		plotOutput("imgSubTree"),
		DT::DTOutput("tblSubTree"),
	)
	)
}
