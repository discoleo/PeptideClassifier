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
## draft v.0.1a


### Count Leafs
# - on each branch;
count.nodes = function(x) {
	x = x$merge;
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
subtree.nc = function(x, n, tree) {
	v  = count.nodes(tree);
	n0 = - x;
	x  = tree$merge;
	# Start from Leaf:
	id0 = which(x[,1] == n0 | x[,2] == n0);
	nnQ = c(id0); id = id0; print(id0)
	if(x[id, 1] > 0) nnQ = c(nnQ, x[id, 1]);
	if(x[id, 2] > 0) nnQ = c(nnQ, x[id, 2]);
	nS = 2;
	while(length(nnQ) >= nS) {
		id = nnQ[nS];
		nS = nS + 1;
		if(x[id,1] > 0) nnQ = c(nnQ, x[id, 1]);
		if(x[id,2] > 0) nnQ = c(nnQ, x[id, 2]);
	}
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
# Repair Tree
# - for extracted Subtree;
# Note: may be meaningful to construct de-novo the tree;
subtree.nn = function(x, tree) {
	x  = sort(unique(x));
	id = rep(0, nrow(tree$merge));
	id[x] = 1; id = cumsum(id);
	nn = tree$merge[x,];
	oldN = - nn[nn < 0];
	nn[nn[,1] > 0, 1] = id[nn[nn[,1] > 0, 1]];
	nn[nn[,2] > 0, 2] = id[nn[nn[,2] > 0, 2]];
	isN = nn < 0; nn0 = - nn[isN];
	idN = order(nn0);
	nn[isN] = - idN;
	tree$merge  = nn;
	tree$height = tree$height[x];
	tree$order = order.tree(tree);
	if(! is.null(tree$labels)) {
		tree$labels = tree$labels[oldN];
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
