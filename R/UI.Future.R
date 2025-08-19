#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## UI: Future Plans
##
## URL: https://github.com/discoleo/PeptideClassifier
##
## draft v.0.2c


### 
panelTODO = function() {
	b = function(x) tag("b", x);
	# UI:
	fluidRow(
		h2("Clustering Methods"),
		"Implement additional robust methods",
		h3("Robust Single linkage"),
		b("Weighted Single Linkage: "),
		"use a weighted distance based on the centrality of each point. A point Pi would be weighted by the distance to the center/centroid of its corresponding cluster, e.g. D|Pi – Centroid| + 1, where the factor + 1 ensures that the distance will not vanish for points very close to the centroid. The weighted distance between 2 points is then the product between the actual distance and the 2 weights.",
		"The weighted distance for points far apart from the corresponding centers is therefore inflated, minimising the impact of such points. The factor +1 could be tuned, e.g. by using the standard error of the mean in the corresponding cluster.",
		fluidRow(""),
		b("Combined Single & Complete Linkage: "),
		"combine the single linkage with the complete linkage and use the sum between the 2 distances.",
		"The maximum distance will counteract the effect of a single pair of very close points.",
		# Shortcomings:
		h2("Shortcomings"),
		h3("Average Link:"),
		"The method seems to sequentially add small clusters to a big one.",
		"It merges very rarely 2 larger clusters! This is easily seen on sub-trees and
		when pruning a tree: a pruning size of over 100 is needed to get a sparser tree.",
		"It may be less suited to classify biologically active peptides and infer similar receptor affinities."
	)
}

