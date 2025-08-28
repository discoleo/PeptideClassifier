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
		fluidRow(
		column(6, selectClusterType()),
		# Plot Orientation:
		column(6, selectClusterPlotOrientation())
		),
		### Generate Tree
		# Existing Tree:
		fileInput.rds("loadTree", "Load existing Tree"),
		# New Tree:
		fluidRow(
		column(8,
			actionButton("btnTreeBuild", "Build Tree"),
			actionButton("btnTreePruneSize", "Prune"),
		)),
		NBSP(),
		fluidRow(
		column(6, selectClusterSource()),
		column(6, textInput("inTreeMinSize",
			"Prune Branches: Min", "20")),
		),
		# SubTrees:
		h3("Subtree:"),
		fluidRow(
		column(8, textInput("txtSubTree_Node", "Node/Peptide", "")),
		column(4, textInput("txtSubTree_Size", "Size", "50")),
		),
		actionButton("btnSubtree", "SubTree"),
		actionButton("btnSubT_Details", "Peptides"),
		# SubTree Data:
		downloadButton("downloadSubT_Data", "Data"),
		downloadButton("downloadTree", "Full Tree"),
		# Messages
		h3(textOutput("txtTreeWarn")),
		textOutput("txtTreeInfo"),
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
panelClusterDiagnostics = function() {
	# UI:
	sidebarLayout(
	sidebarPanel(
		h3("Diagnostics: Branches"),
		### Summary: Branches
		# Summary: Tree computed on the <Cluster>-Tab
		actionButton("btnDx_BranchSummary_1", "Cluster-Tab"),
		# Summary: All Clustering Methods
		actionButton("btnDx_BranchSummary", "All Methods"),
		# Download Trees
		downloadButton("downloadTrees_All", "All Trees"),
		NBSP(),
		### Dx Info
		tag("b", "Branch Ratio:"),
		"Ratio between number of leaves on larger branch vs smaller one.",
		fluidRow(),
		# Branches with only Leaves:
		tag("b", "NBr0 ="),
		"Number of branches with 2 leaves;",
		NBSP(),
		### Load Trees:
		fileInput.rds("loadTrees_Dx", "Load existing Trees"),
	),
	mainPanel(
		DT::DTOutput("tblDx_BranchSummary"),
		plotOutput("imgDx_HistBranchRatios"),
		plotOutput("imgDx_HistBranchLeafs"),
	)
	)
}

### Tab: Correlations
panelClusterCorrelations = function(img.height = 560) {
	hImg = as.img.height(img.height);
	# UI:
	sidebarLayout(
	sidebarPanel(
		h3("Correlations:"),
		fluidRow(
		column(6, selectDxCorTypes()),
		column(6, selectDxCorOrder()),
		),
		actionButton("btnTreeCor", "Correlation"),
		downloadButton("downloadTreeCor", "Correlation"),
		NBSP(),
		fileInput.data("loadTreeCor", "Select csv file with Correlations"),
		# Messages:
		fluidRow(h3(textOutput("txtTreeDx_Warn"))),
		fluidRow(textOutput("txtTreeDx_Info")),
	),
	mainPanel(
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

# Clustering Source:
# - Source of Data-Matrix for DIST;
selectClusterSource = function(selected = "DTM", id = "fltTreeDataSource") {
	selectInput(id, label = "Data Source:",
		choices = list(
			"DTM" = "DTM", "Topic Model" = "TM",
			"Topic n" = "Tn", "DTM for Tn" = "DTM.Tn"),
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
