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
## draft v.0.1c


### Count Leafs
# - on each branch;
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


### Extract Subtree

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
			# TODO: re-design / keep node & makes Leafs;
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
	# TODO: save collapsed Nodes (nnL list);
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
		attr(res, "PnlyLeaves") = sum(isOnlyLeaves);
	}
	return(res);
}


### Size of Leaf-Branches
# Leaf-Branch = Branch to which a leaf joins;
#
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


##############

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

