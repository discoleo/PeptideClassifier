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
	hImg = as.img.height(img.height);
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

### Tab: Diagnostics
panelClusterDiagnostics = function(img.height = 640) {
	hImg = as.img.height(img.height);
	# UI:
	sidebarLayout(
	sidebarPanel(
		h3("Diagnostics:"),
		fluidRow(
		column(6, selectDxCorTypes()),
		column(6, selectDxCorOrder()),
		),
		actionButton("btnTreeCor", "Correlation"),
	),
	mainPanel(
		fluidRow(textOutput("txtTreeDx_Warn")),
		plotOutput("imgTreeCor", height = hImg),
	)
	)
}

### Helper:

as.img.height = function(x) {
	if(is.character(x)) x else paste0(x, "px");
}

### Controls

# Options: Plot Orientation
selectClusterPlotOrientation = function(selected = TRUE, id = "fltTreePlotOrientation") {
	selectInput(id, label = "Plot orientation:",
		choices = list(
			"Horizontal" = TRUE,  "Vertical" = FALSE),
		selected = selected);
}

# Clustering Method:
selectClusterType = function(selected = "complete", id = "fltTreeType") {
	selectInput(id, label = "Clustering method:",
		choices = list(
			"Single"  = "single",  "Complete" = "complete",
			"Ward D"  = "ward.D",  "Ward D2"  = "ward.D2",
			"Average" = "average", "McQuitty/WPGMA" = "mcquitty",
			"Median/WPGMC" = "median", "Centroid/UPGMC" = "centroid"),
		selected = selected);
}

### Diagnostics:

# Correlation Plot: Visualisation Type
selectDxCorTypes = function(selected = "circle", id = "fltDxCorTypes") {
	selectInput(id, label = "Visualisation Type:",
		choices = list(
			"Circle" = "circle", "Ellipse" = "ellipse",
			"Number" = "number", "Pie" = "pie"),
		selected = selected);
}

# Correlation Plot: Order
selectDxCorOrder = function(selected = "AOE", id = "fltDxCorOrder") {
	selectInput(id, label = "Visualisation Type:",
		choices = list(
			"AOE" = "AOE", "FPC" = "FPC",
			"Clustering" = "hclust", "Original" = "original"),
		selected = selected);
}
