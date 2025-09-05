#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## URL: https://github.com/discoleo/PeptideClassifier
##
## Clustering Tools: Graphics
##
## draft v.0.1d


### Plot

### Plot Tree/SubTree
#' @export
plot.subtree = function(x, mark = TRUE, lwd = 3,
		col = "#FF243680", adj = c(0.25)) {
	node = attr(x, "N0");
	lblN = paste0("Subtree: ", node);
	# xPos:
	xPos   = match(as.character(node), x$labels);
	isNode = ! is.na(xPos);
	height = max(x$height);
	xPos   = which(x$order == xPos);
	# Plot:
	plot(x, xlab = lblN);
	if(mark && isNode) lines(c(xPos, xPos) + adj, c(0, height),
		lwd = 3, col = col);
	invisible(c(xPos, height));
}

### Aggregate Branches
# x  = Tree;
# by = Attribute/Value of leaves to aggregate;
#' @export
aggregate.tree = function(x, by, FUN, ...) {
	UseMethod("aggregate.tree");
}

#' @exportS3Method aggregate.tree hclust
aggregate.tree.hclust = function(x, by, FUN = NULL, ...) {
	vb = by;
	if(! is.null(x$labels)) {
		id = as.numeric(x$labels);
		vb = vb[id];
	}
	# Levels:
	iCL = as.factor(vb);
	lvl = levels(iCL);
	if(is.numeric(vb)) {
		if(is.integer(vb)) {
			lvl = as.integer(lvl);
		} else {
			lvl = as.numeric(lvl);
		}
	}
	nCL = length(lvl);
	iCL = as.integer(iCL);
	x = x$merge;
	n = nrow(x);
	# Result:
	m = matrix(0, nrow = nCL, ncol = n);
	if(n == 0) return(m);
	for(id in seq(n)) {
		idN = x[id, 1];
		if(idN < 0) {
			idC = iCL[- idN];
			m[idC, id] = m[idC, id] + 1;
		} else {
			m[, id] = m[, id] + m[, idN];
		}
		idN = x[id, 2];
		if(idN < 0) {
			idC = iCL[- idN];
			m[idC, id] = m[idC, id] + 1;
		} else {
			m[, id] = m[, id] + m[, idN];
		}
	}
	# Original Value:
	if(! is.null(FUN)) {
		m = apply(m, 2, FUN);
		m = lvl[m];
	}
	attr(m, "col") = vb;
	attr(m, "levels") = lvl;
	return(m);
}

### Edge IDs
# x = Tree;
#' @export
as.edgeId = function(x) {
	nn  = x$merge;
	LEN = nrow(nn);
	lst = list();
	if(LEN == 0) return(lst);
	# Warning: naive implementation; NO bound checks;
	for(id in seq(LEN)) {
		tmp = c();
		lst[[id]] = numeric(0);
		idS = nn[id, 1];
		if(idS < 0) {
			tmp = c(tmp, idS, id);
		} else {
			tmp = c(tmp, lst[[idS]], id);
			lst[[idS]] = numeric(0); # clean-up;
		}
		idS = nn[id, 2];
		# Note: reverted
		# TODO: sometimes has normal orientation;
		if(idS < 0) {
			tmp = c(idS, id, tmp);
		} else {
			tmp = c(lst[[idS]], id, tmp);
			lst[[idS]] = numeric(0);
		}
		lst[[id]] = tmp;
	}
	lst = lst[[LEN]];
	return(lst);
}


### Coloured Dendrogram
# TODO: some rework of Midpoint;
#' @export
as.dendrogramCol = function(x, col, scale = 1) {
	cn  = count.nodes(x);
	nn  = x$merge;
	LEN = nrow(nn);
	lst = list();
	if(LEN == 0) return(lst);
	colN = attr(col, "col");
	# Warning: naive implementation; NO bound checks;
	asLeaf = function(id, col = NULL) {
		tmp = list(id);
		attr(tmp, "label")   = as.character(id);
		attr(tmp, "members") = 1;
		attr(tmp, "height")  = 0;
		attr(tmp, "leaf")    = TRUE;
		if(! is.null(colN)) {
			col = colN[id];
		}
		attr(tmp, "nodePar") = list(col = col);
		attr(tmp, "edgePar") = list(col = col);
		return(tmp);
	}
	for(id in seq(LEN)) {
		tmp = list();
		idS = nn[id, 1];
		if(idS < 0) {
			tmp[[1]] = asLeaf(-idS, col[id]);
			pp1 = (match(- idS, x$order) - 1) * scale;
		} else {
			tmp[[1]]   = lst[[idS]];
			lst[[idS]] = numeric(0);
			pp1 = attr(tmp[[1]], "midID");
		}
		idS = nn[id, 2];
		if(idS < 0) {
			tmp[[2]] = asLeaf(-idS, col[id]);
			pp2 = (match(- idS, x$order) - 1) * scale;
		} else {
			tmp[[2]]   = lst[[idS]];
			lst[[idS]] = numeric(0);
			pp2 = attr(tmp[[2]], "midID");
		}
		if(pp1 > pp2) {
			tmpN = tmp[[1]]; tmp[[1]] = tmp[[2]]; tmp[[2]] = tmpN;
		}
		# MidPoint:
		midL = attr(tmp[[1]], "members");
		midR = attr(tmp[[2]], "members");
		mid1 = if(midL <= 2) 0.5
			else attr(tmp[[1]], "midpoint");
		mid2 = if(midR <= 1) -0.5
			else if(midR == 2) 0
			else attr(tmp[[2]], "midpoint");
		midB = (mid1 + midL + mid2) / 2; # Old: (cn[id] - 1) / 2;
		attr(tmp, "members")  = cn[id];
		attr(tmp, "height")   = x$height[id];
		attr(tmp, "midpoint") = midB;
		attr(tmp, "midID")    = (pp1 + pp2) / 2;
		attr(tmp, "edgePar")  = list(col = col[id]);
		lst[[id]] = tmp;
	}
	lst = lst[[LEN]];
	class(lst) = "dendrogram";
	return(lst);
}


### Identify Solitary Leafs
# x   = Object of type hclust;
# col = Colour used for solitary leaves;
#     or 2 colours for solitary vs 2-Leaf branch;
#' @export
as.dendrogramLeaf = function(x, col, h.leaf = 0) {
	if(inherits(x, "collapsedTree")) {
		# TODO:
		warning("Collapsed Tree not yet implemented!");
		return(as.dendrogram(x));
	}
	nn  = x$merge;
	LEN = nrow(nn);
	lst = list();
	if(LEN == 0) return(lst);
	if(length(col) == 1) col = c(L1 = col, L2 = NA);
	cn = count.nodes(x);
	# Warning: naive implementation; NO bound checks;
	asLeaf = function(id, col = NA) {
		tmp = list(id);
		attr(tmp, "label")   = as.character(id);
		attr(tmp, "members") = 1;
		attr(tmp, "height")  = h.leaf;
		attr(tmp, "leaf")    = TRUE;
		if(! is.na(col)) {
			attr(tmp, "nodePar") = list(col = col);
			attr(tmp, "edgePar") = list(col = col);
		}
		return(tmp);
	}
	for(id in seq(LEN)) {
		tmp = list();
		idS = nn[id, 1];
		isL2 = FALSE;
		if(idS < 0) {
			if(nn[id, 2] < 0) isL2 = TRUE; # 2 Leaves;
			colL = if(isL2) col[2] else col[1];
			tmp[[1]] = asLeaf(-idS, colL);
			pp1  = match(- idS, x$order);
		} else {
			tmp[[1]]   = lst[[idS]];
			lst[[idS]] = numeric(0);
			pp1 = attr(tmp[[1]], "midID");
		}
		idS = nn[id, 2];
		if(idS < 0) {
			colL = if(isL2) col[2] else col[1];
			tmp[[2]] = asLeaf(-idS, colL);
			pp2  = match(- idS, x$order);
		} else {
			tmp[[2]]   = lst[[idS]];
			lst[[idS]] = numeric(0);
			pp2 = attr(tmp[[2]], "midID");
		}
		if(pp1 > pp2) {
			tmpN = tmp[[1]]; tmp[[1]] = tmp[[2]]; tmp[[2]] = tmpN;
		}
		# MidPoint:
		midL = attr(tmp[[1]], "members");
		midR = attr(tmp[[2]], "members");
		mid1 = if(midL <= 2) 0.5
			else attr(tmp[[1]], "midpoint");
		mid2 = if(midR <= 1) -0.5
			else if(midR == 2) 0
			else attr(tmp[[2]], "midpoint");
		midB = (mid1 + midL + mid2) / 2; # Old: (cn[id] - 1) / 2;
		attr(tmp, "members")  = cn[id];
		attr(tmp, "height")   = x$height[id];
		attr(tmp, "midpoint") = midB;
		attr(tmp, "midID")    = (pp1 + pp2) / 2;
		# attr(tmp, "edgePar")  = list(col = col[id]);
		lst[[id]] = tmp;
	}
	lst = lst[[LEN]];
	class(lst) = "dendrogram";
	return(lst);
}
