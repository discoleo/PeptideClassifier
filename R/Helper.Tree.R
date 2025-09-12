#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## URL: https://github.com/discoleo/PeptideClassifier
##
## Clustering Tools
##
## draft v.0.1d


#' @export
as.TreeList = function(x, force = FALSE) {
	if(inherits(x, "TreeList")) return(x);
	if(inherits(x, "hclust")) {
		nm = x$method;
		x  = list(x); names(x) = nm;
	} else if(inherits(x, "list")) {
	} else {
		warning("Object does NOT seem to be a valid list of trees!");
		if(! force) return(NULL);
	}
	class(x) = c("TreeList", class(x));
	return(x);
}


### Tree Topology

# Leaf-Branch = Branch to which at least 1 leaf joins;


### Count Leafs
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
# - on Leaf-Branches;
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

#############
### Sub-Trees

### Extract Subtree
# x = Leaf included in the subtree;
# n = Minimal Number of leafs in subtree;
subtree.nc = function(x, n, tree, debug = FALSE) {
	v  = count.nodes(tree);
	n0 = - x;
	x  = tree$merge;
	# Minimalistic checks:
	LEN   = nrow(x);
	idErr = which(x[,1] >= seq(LEN));
	if(length(idErr) > 0) {
		x[idErr, 1] = -Inf;
	}
	idErr = which(x[,2] >= seq(LEN));
	if(length(idErr) > 0) {
		x[idErr, 2] = -Inf;
	}
	# Start from Leaf:
	id0 = which(x[,1] == n0 | x[,2] == n0);
	nnQ = c(id0); id = id0;
	if(debug) print(id0);
	# Add descendants:
	if(x[id, 1] > 0) nnQ = c(nnQ, x[id, 1]);
	if(x[id, 2] > 0) nnQ = c(nnQ, x[id, 2]);
	nS = 2;
	while(length(nnQ) >= nS) {
		id = nnQ[nS];
		nS = nS + 1;
		if(x[id,1] > 0) nnQ = c(nnQ, x[id, 1]);
		if(x[id,2] > 0) nnQ = c(nnQ, x[id, 2]);
	}
	if(debug) { cat("N0 + Desc: "); print(head(nnQ)); }
	if(v[id0] >= n) {
		subT = subtree.nn(nnQ, tree);
		return(subT);
	}
	subT = function(id) {
		nnQ = c(id); nS = 1;
		while(length(nnQ) >= nS) {
			id = nnQ[nS];
			nS = nS + 1;
			if(x[id,1] > 0) nnQ = c(nnQ, x[id, 1]);
			if(x[id,2] > 0) nnQ = c(nnQ, x[id, 2]);
		}
		return(nnQ);
	}
	idL = id0;
	for(id in seq(id0, nrow(x))) {
		if(x[id,1] == idL) {
			nnQ = c(nnQ, id);
			if(x[id,2] > 0) {
				# Add subtree:
				nnQ = c(nnQ, subT(x[id,2]));
			}
			if(v[id] >= n) break;
			idL = id;
			next;
		}
		if(x[id,2] == idL) {
			nnQ = c(nnQ, id);
			if(x[id,1] > 0) {
				# Add subtree:
				nnQ = c(nnQ, subT(x[id,1]));
			}
			if(v[id] >= n) break;
			idL = id;
			next;
		}
	}
	subT = subtree.nn(nnQ, tree);
	return(subT);
}

### Collapse/Prune Tree
# n = Min-size of branch;
#   where size = count(Leaves);
collapse.tree = function(n, tree) {
	nn = count.nodes(tree);
	x  = tree$merge;
	LEN = nrow(x); MAX = max(tree$order);
	nnQ = LEN; nposQ = 1;
	nnT = c(); nnL = list();
	# Warning: naive implementation; NO bound checks;
	while(nposQ <= length(nnQ)) {
		id = nnQ[nposQ];
		if(nn[id] >= n) {
			nnT = c(id, nnT);
			# Keep Leaves; Add Sub-nodes to Queue;
			# Verify SubTrees:
			hasLeafs = TRUE; nmN = NULL;
			if(x[id,1] > 0) {
				if(nn[x[id,1]] >= n) {
					# nnT = c(nnT, x[id,1]);
					nnQ[nposQ] = x[id,1];
					hasLeafs = FALSE;
				} else {
					nnLeafs = collect.nodes(x[id,1], x);
					MAX = MAX + 1;
					x[id,1] = - MAX;
					nmN = paste0("L", id);
					nnL[[nmN]] = list(Lvl = id, N = nnLeafs);
				}
			}
			if(x[id,2] > 0) {
				if(nn[x[id,2]] >= n) {
					# nnT = c(nnT, x[id,2]);
					if(hasLeafs) {
						nnQ[nposQ] = x[id,2];
						hasLeafs = FALSE;
					} else {
						nnQ = c(nnQ, x[id,2]);
					}
				} else {
					nnLeafs = collect.nodes(x[id,2], x);
					MAX = MAX + 1;
					x[id,2] = - MAX;
					if(is.null(nmN)) {
						nmN = paste0("L", id);
						nnL[[nmN]] = list(Lvl = id, N = nnLeafs);
					} else nnL[[nmN]]$N = c(nnL[[nmN]]$N, nnLeafs);
				}
			}
			if(hasLeafs) nposQ = nposQ + 1; # Next node in Queue;
		} else {
			# Probably the Root-Node < n;
			# TODO: re-design / keep node & make Leaves;
			print("Root fails.")
			nposQ = nposQ + 1;
			nmN   = paste0("L", id);
			nnL[nmN] = list(Lvl = id, N = c());
			if(x[id,1] < 0) {
				nnL[[nmN]]$N = c(nnL[[nmN]]$N, x[id,1]);
			} else {
				allL = collect.nodes(id, x);
				nnL[[nmN]]$N = c(nnL[[nmN]]$N, allL);
			}
			if(x[id,2] < 0) {
				nnL[[nmN]]$N = c(nnL[[nmN]]$N, x[id,2]);
			} else {
				allL = collect.nodes(x[id,2], x);
				nnL[[nmN]]$N = c(nnL[[nmN]]$N, allL);
			}
			MAX = MAX + 2;
			x[id,] = c(-MAX + 1, -MAX);
			nnT = c(nnT, id);
		}
	}
	# Extract Tree:
	nnT = sort(unique(nnT));
	tree$merge = x;
	tree = subtree.fromRoot(nnT, tree);
	# Collapsed Nodes (nnL list);
	# TODO: convert nnL to matrix form;
	attr(tree, "nodes") = nnL;
	class(tree) = c("collapsedTree", class(tree));
	return(tree);
}


### Repair Tree
# - for extracted Subtree;
# Note: may be meaningful to construct de-novo the tree;
# x = Rows to keep from the tree$merge matrix;
subtree.nn = function(x, tree) {
	x  = sort(unique(x));
	id = rep(0, nrow(tree$merge));
	id[x] = 1; id = cumsum(id);
	nn = tree$merge[x,];
	isLf  = nn < 0;
	oldN  = - nn[isLf];
	isRef = nn[,1] > 0; nn[isRef, 1] = id[nn[isRef, 1]];
	isRef = nn[,2] > 0; nn[isRef, 2] = id[nn[isRef, 2]];
	# nn0 = - nn[isLf];
	# idN = seq(order(nn0));
	idN   = seq(sum(isLf));
	nn[isLf] = - idN;
	tree$merge  = nn;
	tree$height = tree$height[x];
	tree$order = order.tree(tree);
	if(! is.null(tree$labels)) {
		tree$labels = tree$labels[oldN];
	}
	return(tree);
}

# Repair pruned SubTree (from root)
# x = Rows to keep from the tree$merge matrix;
subtree.fromRoot = function(pos, tree) {
	nnT = pos; x = tree$merge;
	TOP = length(tree$order);
	T0 = rep(0, nrow(x));
	T0[nnT] = 1;
	T0 = cumsum(T0);
	T0 = c(T0, T0);
	x  = x[nnT,]; # New Tree;
	oldN = - x[x < 0];
	oldN = oldN[oldN <= TOP];
	nLen = sum(x < 0); # Leafs
	x[x < 0] = - seq(nLen);
	x[x > 0] = T0[x[x > 0]];
	# New Tree:
	tree$merge  = x;
	tree$height = tree$height[nnT];
	tree$order = order.tree(tree);
	if(! is.null(tree$labels)) {
		# TODO
		tree$labels = tree$labels[oldN];
		LEN = length(oldN);
		if(LEN < nLen) {
			tree$labels = c(tree$labels, seq(LEN + 1, nLen));
		}
	}
	return(tree);
}

### Order of Nodes in Tree
# - Repair order in Subtree;
order.tree = function(x) {
	nn  = x$merge;
	LEN = nrow(nn);
	nnQ = LEN; nS = 1;
	nOrder = c();
	# Warning: naive implementation; NO bound checks;
	while(nS <= length(nnQ)) {
		id = nnQ[nS];
		hasLeaf = TRUE;
		if(nn[id,1] < 0) {
			nOrder = c(nOrder, - nn[id,1]);
		} else {
			nnQ[nS] = nn[id,1];
			hasLeaf = FALSE;
		}
		if(nn[id,2] < 0) {
			nOrder = c(nOrder, - nn[id,2]);
			if(hasLeaf) nS = nS + 1;
		} else if(hasLeaf) {
			# same nS;
			nnQ[nS] = nn[id,2];
		} else {
			lenQ  = length(nnQ);
			tailQ = if(lenQ <= nS) c() else nnQ[seq(nS+1, lenQ)];
			nnQ = c(nnQ[1:nS], nn[id,2], tailQ);
		}
	}
	return(nOrder);
}

### Aggregate Functions

### Collect all Leaves on a Branch
# x = Matrix of Nodes;
collect.nodes = function(id, x) {
	nQ = c(id); nnL = c();
	nS = 1;
	while(nS <= length(nQ)) {
		id = nQ[nS];
		hasLeaf = TRUE;
		if(x[id,1] < 0) {
			nnL = c(nnL, x[id,1]);
		} else {
			nQ[nS] = x[id,1];
			hasLeaf = FALSE;
		}
		if(x[id,2] < 0) {
			nnL = c(nnL, x[id,2]);
		} else {
			if(hasLeaf) {
				nQ[nS] = x[id,2];
				hasLeaf = FALSE;
			} else {
				# Order does NOT matter;
				nQ = c(nQ, x[id,2]);
			}
		}
		if(hasLeaf) nS = nS + 1;
	}
	return(nnL);
}

###############

### Analysis


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
#   (informational entropy ~ 30%)
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
	return(ent);
}

###############

#########
### Tools

### Remove Inversions
# - Adjust height in order to remove inversions;
# x = Tree (hclust);
# scale, pow = affect how much the inverted branch will descend;
#' @export
adjust.height = function(x, ...) {
	UseMethod("adjust.height");
}

#' @exportS3Method adjust.height hclust
adjust.height.hclust = function(x, scale = 1/16, pow = 1) {
	nn  = x$merge;
	LEN = nrow(nn);
	Hm = diff(range(x$height));
	Hm = Hm * scale;
	hasPow = pow != 1;
	for(id in seq(LEN, 1)) {
		h   = x$height[id];
		idN = nn[id, ];
		idN = idN[idN > 0];
		if(length(idN) > 0) {
			isBig = x$height[idN] > h;
			if(any(isBig)) {
				idN = idN[isBig];
				div = LEN + 1 - id;
				if(hasPow) div = div^pow;
				x$height[idN] = h - Hm / div;
			}
		}
	}
	invisible(x);
}

##############
