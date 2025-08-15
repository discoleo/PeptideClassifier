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
panelClustering = function(img.height = 800) {
	hImg = if(is.character(img.height)) img.height
		else paste0(img.height, "px");
	# UI:
	sidebarLayout(
	sidebarPanel(
	h3("Subtree:"),
	fluidRow(
	column(4, textInput("txtSubTree_Node", "Node/Peptide", "")),
	column(4, textInput("txtSubTree_Size", "Size", "50")),
	),
	actionButton("btnSubtree", "SubTree"),
	actionButton("btnSubT_Details", "Peptides"),
	# SubTree Data:
	downloadButton("downloadSubT_Data", "Data"),
	fluidRow(HTML("&nbsp;")),
	fluidRow(
	column(6, selectClusterType()),
	# Plot Orientation:
	column(6, selectClusterPlotOrientation())
	),
	),
	mainPanel(
		plotOutput("imgTree", height = hImg),
		DT::DTOutput("tblClusters"),
		plotOutput("imgSubTree"),
		DT::DTOutput("tblSubTree"),
	)
	)
}

selectClusterPlotOrientation = function(selected = TRUE, id = "fltTreePlotOrientation") {
	selectInput(id, label = "Plot orientation:",
		choices = list(
			"Horizontal" = TRUE,  "Vertical" = FALSE),
		selected = selected);
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
