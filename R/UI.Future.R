#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## UI: Tutorial & Future Plans
##
## URL: https://github.com/discoleo/PeptideClassifier
##
## draft v.0.2e

style = function(other = NULL) {
	s = "color:#F03232;";
	if(! is.null(other)) {
		other = paste(other, collapse = ";");
		s = paste(s, other);
	};
	return(s);
}
styleB = function() {
	"color:#4432D8;weight:heavy;font-size:150%;";
}
styleB2 = function() {
	"color:#6444B2;";
}


### Tutorial

panelTutorial = function() {
	b = function(x) tag("b", x);
	# UI:
	rowUI = fluidRow(
		style = "padding-left:20px;",
		### Analysis:
		h2("Analysis"),
		
		### Cluster Types:
		h3("Cluster Types"),
		p(span("Micro-Clusters:", style = styleB()),
		"A micro-cluster is a cluster formed by exactly 2 leaves."
		),
		p(span("Mini-Clusters:", style = styleB()),
		"A mini-cluster includes more than 2 leaves, but far less than O(n), where n is the number of leaves in the tree.",
		"A good threshold is log(n): in a tree formed from 1 000 leaves, mini-clusters range between 3 and 10 leaves."
		),
		p(span("Macro-Clusters:", style = styleB()),
		"Macro-clusters can be defined as clusters having at least O(n^p) leaves, where p is a fractional parameter.",
		"Possible thresholds for the number of leaves are sqrt(n) or n^(1/3); scaling these values by a constant may be useful as well."
		),
		p(span("Medium-Sized Clusters:", style = styleB()),
		"Every cluster with the number of leaves between O(log(n)) and O(n^p) can be described as a medium sized cluster."
		),
		
		# Trees:
		h4("Balanced Trees", style = styleB2()),
		"A balanced tree will have sufficiently many medium sized clusters, which join to form non-nested macro-clusters.",
		h4("Unbalanced Trees", style = styleB2()),
		"In some types of trees, a single macro-cluster increases steadily in size: such clusters are not useful for the current project.",
	);
	tagList(rowUI,
		panelTutorial_Linkage(),
		panelTutorial_Shortcomings());
}

panelTutorial_Linkage = function() {
	b = function(x) tag("b", x);
	# UI:
	fluidRow(
		style = "padding-left:20px;",
		### Linkage Types:
		h3("Linkage Type:"),
		# Ward D1 & D2:
		h3("Ward Method:", style = style()),
		"The Ward method implemented initially in function hclust",
		"requires the squaring of the dissimilarity matrix,",
		"while the true node height is the square root of the computed height.",
		fluidRow(
			style = "padding-left:40px;",
			"r.clust = hclust(dist(...)^2, method = \"ward.D\");"),
		fluidRow(
			style = "padding-left:40px;",
			"r.clust$height = sqrt(r.clust$height);"),
		"The R documentation lacks this information!",
		"The above steps are performed automatically using the ward.D2 method.",
		"The square-root of the height is now used to plot a tree generated with method ward.D,",
		"in order to maintain a proper scaling for nodes close to the root.",
		"The actual power can be set using the global option: cluster.Ward.hpow (see in Server.R).",
		"Notice that the dissimilarity matrix is not squared!",
		"Users are advised to use method ward.D2 if they want",
		"the correct minimal-variance partitioning.",
		fluidRow(
		style = "padding-left:20px;",
		span("Note:", style = style()),
		"Ward method tends to produce very homogenous clusters,",
		  "which may NOT reflect real-life data sets.",
		"Other linkage methods, like complete linkage,",
		  "may create more realistic partitions (for this project).",
		),
		# McQuitty:
		h3("McQuitty Method:", style = style()),
		"Around half of the leaves join other leaves to form a micro-cluster.",
		"The remaining leaves are solitary leaves which join a micro- or mini-cluster.",
		"Only very few solitary leaves join a macro-cluster.",
		"The micro-clusters merge together to form larger mini-clusters and then macro-clusters.",
		"There is a good partitioning into medium-sized and large-sized clusters.",
	)
}

panelTutorial_Shortcomings = function() {
	b = function(x) tag("b", x);
	### Shortcomings:
	fluidRow(
		style = "padding-left:20px;",
		### Shortcomings:
		h2("Shortcomings"),
		# Median Linkage:
		h3("Median Linkage:", style = style()),
		"There are very few micro-clusters (2 leaves), a few mini-clusters and virtually no macro-clusters.",
		"Almost all leaves join the tree as solitary leaves.",
		"The solitary leaves join sequentially to a single very large \"branch\".",
		
		# Centroid Linkage:
		h3("Centroid Linkage:", style = style()),
		"There are very few macro-clusters, and, with the exception of the main branch, all are tiny.",
		"Most of them are actually at the limit between medium-sized and macro-clusters.",
		"Most leaves join one single branch as solitary leaves.",
		"This branch grows steadily.",
		"Another issue is the presence of (large numbers of) inversions.",
		
		# Average Linkage:
		h3("Average Linkage:", style = style()),
		"There are very few macro-clusters.",
		"Many micro- and mini-clusters join one of the large branches, which grows sequentially to a larger size.",
		"Most leaves form actually micro-clusters, which merge then to mini-clusters.",
		"However, the method seems to sequentially add many of these small clusters to a bigger branch.",
		"It merges very rarely 2 larger clusters! This is easily seen on sub-trees and
		when pruning a tree: a pruning size of over 100 is needed to get a sparser tree.",
		"It may be less suited to classify biologically active peptides and infer similar receptor affinities.",
		
		# Inversions:
		h3("Inversions:", style = style()),
		"Inversions are common with certain linkage types, e.g. with centroid linkage.",
		fluidRow(plotOutput("imgTutor_Inversions")),
		NBSP()
	);
}

### Future-Ideas:
panelTODO = function() {
	b = function(x) tag("b", x);
	# UI:
	fluidRow(
		style = "padding-left:20px;",
		### Clustering:
		h2("Clustering Methods", style = style()),
		"Implement additional robust methods",
		h3("Robust Single Linkage"),
		b("Weighted Single Linkage: "),
		"use a weighted distance based on the centrality of each point. A point Pi would be weighted by the distance to the center/centroid of its corresponding cluster, e.g. D|Pi – Centroid| + 1, where the factor + 1 ensures that the distance will not vanish for points very close to the centroid. The weighted distance between 2 points is then the product between the actual distance and the 2 weights.",
		"The weighted distance for points far apart from the corresponding centers is therefore inflated, minimising the impact of such points. The factor +1 could be tuned, e.g. by using the standard error of the mean in the corresponding cluster.",
		fluidRow(""),
		b("Combined Single & Complete Linkage: "),
		"combine the single linkage with the complete linkage and use the sum between the 2 distances.",
		"The maximum distance will counteract the effect of a single pair of very close points.",
		### Features
		h2("Peptides & Features", style = style()),
		### PTM:
		h3("Post-Translational Modifications of Amino Acids", style = styleB()),
		"Natural peptides contain many post-translational modifications (PTM).",
		"Currently, only the 20 standard amino acids are permitted and used to compute the features.",
		"However, various modifications are common in nature:",
		tags$ul(
		tag("li", "C-terminal amidation;"),
		tag("li", "Hydroxy-proline;"),
		tag("li", "Gamma-carboxy-glutamate;"),
		tag("li", "N-terminal cyclization (pyrD);"),
		tag("li", "Sulfation of tyrosines;"),
		),
		"Such modifications can still be represented with 1 letter codes",
		"and could easily be incorporated into the application.",
		### Disulphide Bonds:
		h3("Disulphide Bonds", style = styleB()),
		"Small peptides usually lack a defined secondary and tertiary structure.",
		"Conversely, disulphide bonds are covalent bonds which can convey rigidity to a structure.",
		"Many conotoxins and spider venoms have one or more such disulphide bonds.",
		"Incorporating these structural details as chemical features",
		"will be quite challenging.",
		### Cyclic Peptides:
		h3("Cyclic Peptides", style = styleB()),
		"Some of the natural peptides are cyclic.",
		"This may require an update to the algorithm which generates the n-grams.",
		"The current application may be easily adapted to accommodate in part such peptides as well.",
		NBSP()
	)
}

