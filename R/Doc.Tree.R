##########################
##
## Tutorials
## Hierarchical Clustering
##
## Leonard Mada

### Various Tutorials
### Hierarchical Clustering


#' @export
link.centroid.data = function() {
	df = data.frame(
		x = c(0,0.6, 0,0.6, 2, 3,3),
		y = c(0,0, 2,2, 1, 0.5,1.5)
	);
	return(df);
}

### Centroids
#' @export
plot.data.centroid = function(col = c("red", "blue"), pch = 16,
		labels.stem = "C", labels.points = "P",
		dxy = 0.5, dC = 0.25) {
	df = link.centroid.data();
	cc = data.frame(
		x  = c(0.3,0.3,0.3,  3, 8/3, 9.2/7),
		y  = c(0,2,1,   1, 1, 1),
		id = c(1,2,5, 3, 4, 6)
	);
	if(length(dxy) == 1) dxy = c(-dxy, dxy);
	limx = range(df$x) + dxy;
	limy = range(df$y) + dxy;
	plot(df, xlim = limx, ylim = limy);
	points(cc, col=col[[1]], pch=pch);
	if(! is.null(labels.points)) {
		lbls = paste(labels.points, seq(nrow(df)), sep = "");
		text(df$x, df$y + dC, labels = lbls, col = col[[2]]);
	}
	if(! is.null(labels.stem)) {
		lbls = paste(labels.stem, cc$id, sep = "");
		text(cc$x, cc$y + dC, labels = lbls, col = col[[1]]);
	}
	invisible();
}

# Inversions:
plot.tree.centroid = function() {
	d = link.centroid.data();
	tree = hclust(dist(d), method = "centroid");
	tree$call[[2]] = "D = Set w. Inversions";
	plot(tree, main = "Dendrogram with Inversions");
	invisible(tree);
}
