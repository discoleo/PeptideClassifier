#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## URL: https://github.com/discoleo/PeptideClassifier
##
## Clustering Tools:
## Analysis & Topology of Trees
##
## draft v.0.1d


### Tools

### Tree Topology

# Leaf-Branch = Branch to which at least 1 leaf joins;
# - at least 1 of the direct descendants is a leaf;


### Count Leaves: All
# - on each branch;
# Note: full count of leaves for ALL nodes;
#' @export
count.nodes = function(x) {
	x = x$merge;
	if(nrow(x) == 0) return(numeric(0));
	v = rep(0, nrow(x));
	v[x[,1] < 0] = 1;
	hasLeaf = x[,2] < 0;
	v[hasLeaf] = v[hasLeaf] + 1;
	# Nodes:
	for(i1 in seq(nrow(x))) {
		if(x[i1,1] > 0) v[i1] = v[i1] + v[x[i1,1]];
		if(x[i1,2] > 0) v[i1] = v[i1] + v[x[i1,2]];
	}
	return(v);
}

### Count Leaves
# - only on Leaf-Branches;
# Note: Branch with 2 leaves = 0;
# nc = Count of leaves on a branch;
#' @export
count.jLeaves = function(x, nc = NULL) {
	if(is.null(nc)) {
		nC = count.jLeaves0(x);
		return(nC);
	}
	x = x$merge;
	LEN = nrow(x);
	if(LEN == 0) return(integer(0));
	nL  = integer(0);
	for(id in seq(LEN)) {
		if(x[id,1] < 0) {
			tmp = nc[id] - 1;
			if(tmp == 1) tmp = 0;
			nL = c(nL, tmp);
		} else if(x[id,2] < 0) {
			tmp = nc[id] - 1;
			nL = c(nL, tmp);
		}
	}
	return(nL);
}
count.jLeaves0 = function(x) {
	# One-pass algorithm;
	x = x$merge;
	LEN = nrow(x);
	if(LEN == 0) return(integer(0));
	# Number of Leaves:
	nL = rep(0, LEN);
	hasLeaf = rep(FALSE, LEN);
	for(id in seq(LEN)) {
		n1 = x[id,1]; n2 = x[id,2];
		if(n1 < 0) {
			hasLeaf[id] = TRUE;
			if(n2 > 0) {
				nL[id] = nL[n2] + 1;
			} else nL[id] = 2;
		} else {
			if(n2 < 0) {
				hasLeaf[id] = TRUE;
				nL[id] = nL[n1] + 1;
			} else {
				nL[id] = nL[n1] + nL[n2];
			}
		}
	}
	nL = nL[hasLeaf];
	nL = nL + ifelse(nL == 2, -2, -1);
	return(nL);
}

### Size of Leaf-Branches
# x = Tree;
# counts = Number of leafs on each branch;
#' @export
size.leafBranch = function(x, counts = NULL) {
	if(is.null(counts)) counts = count.nodes(x);
	x = x$merge;
	LEN = nrow(x);
	if(length(counts) != LEN) {
		stop("Number of nodes differs!");
	}
	if(LEN == 0) return(numeric(0));
	#
	idLeaf = which(x[, 1] < 0);
	res    = rep(1, length(idLeaf));
	# Real Branches:
	idPos  = which(x[idLeaf, 2] > 0);
	if(length(idPos) > 0) {
		idBranch   = idLeaf[idPos];
		res[idPos] = counts[x[idBranch, 2]];
	}
	# 2nd Col:
	idLeaf   = which(x[, 2] < 0);
	idBranch = idLeaf[x[idLeaf, 1] > 0];
	res = c(res, counts[x[idBranch, 1]]);
	return(res);
}


### Count Leaves
# - Useful to compute various tree-Indexes
#   using simpler 1 pass algorithms;
# - Enables computations without the merge-data-structure;
# Out = Matrix with 2 rows;
#' @export
count.leavesD = function(x) {
	x   = x$merge;
	LEN = nrow(x);
	if(LEN < 1) return(matrix(0, nrow=2, ncol=0));
	mR = matrix(0, nrow = 2, ncol = LEN);
	for(id in seq(LEN)) {
		ni = x[id,1];
		nc = if(ni < 0) 1 else sum(mR[, ni]);
		mR[1, id] = nc;
		ni = x[id,2];
		nc = if(ni < 0) 1 else sum(mR[, ni]);
		mR[2, id] = nc;
	}
	# TODO: sort values?
	class(mR) = c("DualCount", class(mR));
	return(mR);
}

###############

### Analysis

### Agglomerative Index

### Average Height of "Leafs"
# x = Object of type hclust;
# Out = Simple average of heights of nodes joined by a leaf;
index.agg1 = function(x, h0 = 0) {
	h  = x$height;
	hL = h[x$merge[,1] < 0];
	hL = c(hL, h[x$merge[,2] < 0]);
	hT = sum(hL);
	if(h0 != 0) hT = hT - h0 * length(hL);
	# Height of Root:
	maxH = h[length(h)] - h0;
	hT = hT / (length(hL) * maxH);
	return(hT);
}
### Sum of All Heights
# - this is probably the Agglomerative coefficient;
index.aggAll = function(x, h0 = 0) {
	h  = x$height;
	hT = sum(h);
	if(h0 != 0) hT = hT - h0 * length(h);
	# Height of Root:
	maxH = h[length(h)] - h0;
	hT   = hT / (length(h) * maxH);
	return(hT);
}
### Harmonic-Weighted Sum of Heights
index.aggHarm = function(x, h0 = 0) {
	h  = x$height;
	nc = count.nodes(x);
	hT = sum((h - h0) / nc);
	hT = hT / sum(1/nc);
	# Height of Root:
	maxH = h[length(h)] - h0;
	hT   = hT / maxH;
	return(hT);
}
### Simple Weighted Sum of Heights
index.aggWSum = function(x, h0 = 0) {
	h  = x$height;
	nc = count.nodes(x); # Counts only leaves;
	hT = sum((h - h0) * nc);
	hT = hT / sum(nc);
	# Height of Root:
	maxH = h[length(h)] - h0;
	hT   = hT / maxH;
	return(hT);
}

### Tree Balance

### Leaf Ratios
# - between number of leafs on each branch;
# Note:
# - Centroid: very, very high 3rd Quartile;
# - Average: high to very high 3rd Quartile;
# - Many leaves (or minute branches) are added sequentially
#   to an ever increasing branch!
#
# x = Tree;
# counts = Number of leafs on each branch;
#' @export
branch.ratios = function(x, counts = NULL, rm.onlyLeaves = TRUE) {
	if(is.null(counts)) counts = count.nodes(x);
	x = x$merge;
	LEN = nrow(x);
	if(length(counts) != LEN) {
		stop("Number of nodes differs!");
	}
	if(LEN == 0) return(numeric(0));
	#
	res = rep(0, LEN);
	for(id in seq(LEN)) {
		if(x[id, 1] > 0) {
			if(x[id, 2] < 0) {
				res[id] = counts[id] - 1;
			} else {
				tmp2 = counts[x[id, 2]];
				tmp1 = counts[id] - tmp2;
				if(tmp1 >= tmp2) {
					res[id] = tmp1 / tmp2;
				} else {
					res[id] = tmp2 / tmp1;
				}
			}
		} else if(x[id, 2] > 0) {
			# Note: x[id, 1] < 0!
			res[id] = counts[id] - 1;
		}
	}
	if(rm.onlyLeaves) {
		isOnlyLeaves = res == 0;
		res = res[! isOnlyLeaves];
		attr(res, "onlyLeaves") = sum(isOnlyLeaves);
	}
	class(res) = c("BranchRatios", class(res));
	return(res);
}

# Summary: Branch Ratios
# x    = List of Trees;
# NBr0 = Branches with only leaves;
# Note:
# - Ratio is artificially set to 0, so that they are distinguished easily
#   from other branches with Ratio = 1;
#' @export
summary.branchRatiosTn = function(x, name.Br0 = "NBr0") {
	LEN = length(x);
	smr = sapply(seq(LEN), function(id) {
		rBr = branch.ratios(x[[id]]);
		tmp = summary(rBr);
		nLeafs = attr(rBr, "onlyLeaves");
		tmp = c(nLeafs, tmp);
		names(tmp)[1] = name.Br0;
		return(tmp);
	});
	smr = data.frame(t(smr))
	smr = cbind(names(x), smr);
	names(smr)[c(1,4,7)] = c("Method", "Q1", "Q3");
	return(smr);
}

### Chaining Coefficient
# Note: similar in concept to the Branch Ratio;
#' @export
index.chaining = function(x) {
	x = x$merge;
	LEN = nrow(x);
	if(LEN < 3) return(0);
	# nc = Leaf count;
	nc  = rep(0, LEN);
	dCh = 0;
	for(id in seq(LEN)) {
		ni = x[id, 1];
		n1 = if(ni < 0) 1 else nc[ni];
		ni = x[id, 2];
		n2 = if(ni < 0) 1 else nc[ni];
		dn = abs(n1 - n2);
		dCh    = dCh + dn;
		nc[id] = n1 + n2;
	}
	# Normalisation:
	dCh = (2*dCh) / (LEN - 1) / (LEN - 2);
	return(dCh);
}


### Entropy
# Note:
# - Ward & Complete linkage seem to do a fair job
#   at generating relatively balanced trees;
#   (informational entropy = 40% - 50%)
#' @export
index.entropy = function(x) {
	x = x$merge;
	LEN = nrow(x);
	if(LEN < 1) return(0);
	if(LEN == 1) return(1);
	log2 = log(2);
	# nc = Leaf count;
	nc  = rep(0, LEN);
	ent = 0;
	for(id in seq(LEN)) {
		ni = x[id, 1];
		if(ni < 0) {
			n1 = 1; e1 = log2;
		} else {
			n1 = nc[ni]; e1 = log(n1);
		}
		ni = x[id, 2];
		if(ni < 0) {
			n2 = 1; e2 = log2;
		} else {
			n2 = nc[ni]; e2 = log(n2);
		}
		nci = n1 + n2;
		ent = ent - (n1*e1 + n2*e2) / nci + log(nci);
		nc[id] = nci;
	}
	# Normalisation:
	ent = ent / LEN;
	ent = ent / log2; # Log-Base = 2;
	return(ent);
}
