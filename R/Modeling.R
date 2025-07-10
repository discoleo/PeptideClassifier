#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## draft v.0.1c


###############

############
### Data ###

### Classic n-Grams

# Demo: Compute set of n-Grams
ngrams.demo = function(x,
		add.3 = FALSE, add.3u = TRUE, add.length = TRUE,
		breaks = c(0, 9, 19, 29, 100), prefix = "_") {
	# nGrams:
	xg2  = ngrams(x, n = 2);
	xg2u = as.ngram.undirected(xg2, prefix = prefix);
	xall = merge.list(xg2u, xg2);
	# 3-Grams
	xg3  = ngrams(x, n = 3)
	xg3u = as.ngram.undirected(xg3, prefix = prefix);
	if(add.3u) {
		xall = merge.list(xall, xg3u);
	}
	if(add.3) {
		xall = merge.list(xall, xg3);
	}
	# Length:
	if(add.length) {
		xall = merge.list(xall, as.list(len.pp(x, breaks = breaks)));
	}
	invisible(xall);
}

ngrams = function(x, n = 2, directed = TRUE) {
	FUN = function(x) {
		len = nchar(x);
		if(len == n) return(x);
		if(len < n) return("");
		len = len - n + 1;
		ngr = character(len);
		for(npos in seq(len)) {
			ngr[npos] = substr(x, npos, npos + n - 1);
		}
		return(ngr);
	}
	ngr = lapply(x, FUN);
	if(! directed) {
		ngr = as.ngram.undirected(ngr);
	}
	return(ngr);
}

### nGrams: Undirected
# e.g. "AC" & "CA" => "_AC";
as.ngram.undirected = function(x, prefix = "_") {
	len = length(x);
	if(len == 0) return(x);
	n  = nchar(x[[1]][1]);
	aa = LETTERS; len.aa = length(aa);
	# Bi-Grams
	if(n == 2) {
		ngr1 = rep(aa[-len.aa], seq(len.aa-1, 1));
		id2  = lapply(seq(len.aa - 1), \(LEN) {
				seq(LEN+1, len.aa);
			});
		ngr2 = aa[unlist(id2)];
		ngrA = paste0(ngr2, ngr1);
		ngrB = paste0(ngr1, ngr2);
		for(id in seq(len)) {
			idRev = match(x[[id]], ngrA);
			isRev = ! is.na(idRev);
			x[[id]][isRev] = ngrB[idRev[isRev]];
		}
		if(! is.null(prefix)) {
			for(id in seq(len)) {
				x[[id]] = paste0(prefix, x[[id]]);
			}
		}
		return(x);
	}
	# 3-Grams
	if(n == 3) {
		ngr1 = rep(aa[-len.aa], seq(len.aa-1, 1));
		id2  = lapply(seq(len.aa - 1), \(LEN) {
				seq(LEN+1, len.aa);
			});
		ngr2 = aa[unlist(id2)];
		ngrA = paste0(ngr2, ngr1);
		for(id in seq(len)) {
			tmp  = x[[id]];
			str1 = substr(tmp, 1, 1);
			str3 = substr(tmp, 3, 3);
			isRev = str3 < str1;
			str1  = str1[isRev]; str3 = str3[isRev];
			x[[id]][isRev] = paste0(str3, substr(tmp[isRev], 2, 2), str1);
		}
		if(! is.null(prefix)) {
			for(id in seq(len)) {
				x[[id]] = paste0(prefix, x[[id]]);
			}
		}
		return(x);
	} else if(n == 4) {
		# 4-Grams
		lst = as.ngram.u4(x, prefix=prefix);
		return(lst);
	}
	# NOT yet;
	warning("Not yet!");
}
as.ngram.u4 = function(x, prefix = "_") {
	len = length(x);
	if(len == 0) return(x);
	for(id in seq(len)) {
		tmp = x[[id]];
		s1 = substr(tmp, 1, 1); s4 = substr(tmp, 4, 4);
		s2 = substr(tmp, 2, 2); s3 = substr(tmp, 3, 3);
		# Case: "BXXA"
		isRev  = s4 < s1;
		tmpRev = tmp[isRev];
		s2 = substr(tmpRev, 2, 2); s3 = substr(tmpRev, 3, 3);
		x[[id]][isRev] = paste0(s4[isRev], s3, s2, s1[isRev]);
		# Case: "sBAs", s = same;
		isRev  = s4 == s1;
		idRev  = which(isRev);
		tmpRev = tmp[idRev];
		s2 = substr(tmpRev, 2, 2); s3 = substr(tmpRev, 3, 3);
		isRev2 = s3 < s2;
		idRev  = idRev[isRev2];
		x[[id]][idRev] = paste0(s4[idRev], s3[isRev2], s2[isRev2], s1[idRev]);
	}
	if(! is.null(prefix)) {
		for(id in seq(len)) {
			x[[id]] = paste0(prefix, x[[id]]);
		}
	}
	return(x);
}

### Peptide Length
# Encode Length as n-Gram;
len.pp = function(x, breaks = c(0, 9, 19, 29, 100),
		prefix = "", sufix = "L") {
	len = nchar(x);
	len = cut(len, breaks = breaks);
	len = as.numeric(len);
	len = paste0(prefix, len, sufix);
	return(len)
}

### Charge
# Compute naive/simple charge of the n-Grams;
# Note:
# - is inherently undirected;
# - for actual charge, use:
#   ngrams.charge( _SEQ_ , breaks = NULL, prefix = NULL);
ngrams.charge = function(x, n = 4, breaks = 5, prefix = "+-") {
	x = lapply(x, function(x) {
		as.numeric(charToRaw(x)) - 64;
	});
	lst = ngrams.charge.numeric(x, n=n, breaks=breaks, prefix=prefix);
	return(lst);
}
ngrams.charge.numeric = function(x, n = 4, breaks = 5, prefix = "+-") {
	x = lapply(x, function(x) {
		len = length(x);
		z = rep(0, len);
		z[x == 4  | x ==  5] = -1;
		z[x == 11 | x == 18] = 1;
		z[x == 8] = 1; # His
		z[1] = z[1] + 1; z[len] = z[len] - 1;
		if(n == 1) return(z);
		if(n >= len) return(sum(z));
		# Proper n-Grams:
		sapply(seq(len - n + 1), function(npos) {
			sum(z[seq(npos, npos + n - 1)]);
		});
	});
	if(! is.null(breaks)) {
		if(length(breaks == 1)) {
			br = range(unlist(x));
			dx = diff(br);
			Ln = breaks + 1;
			br = seq(br[1], br[2], length.out = Ln);
			br[1] = br[1] - 1; br[Ln] = br[Ln] + 1;
			breaks = br;
		}
		x = lapply(x, function(x) {
			as.numeric(cut(x, breaks = breaks));
		});
	}
	if(! is.null(prefix)) {
		x = lapply(x, function(x) {
			paste0(prefix, x);
		});
	}
	return(x);
}


# Helper:

merge.list = function(x, y) {
	len = length(x);
	if(length(y) != len) stop("Length Mismatch!");
	if(len == 0) return(x);
	tmp = lapply(seq(len), function(id) {
		c(x[[id]], y[[id]]);
	});
	return(tmp);
}

###################
###################

#' @export
dtm = function(x, min.len = 1) {
	x.str = sapply(x, \(x) paste0(x, collapse = " "));
	tm::DocumentTermMatrix(x.str,
		control = list(
			stemming = FALSE, stopwords = FALSE, tolower = FALSE,
			minWordLength = min.len,
			removeNumbers = FALSE, removePunctuation = FALSE));
}

#' @export
dim.dtm = function(x) {
	dim(x);
}

#' @export
freq.term = function(x, data) {
	id = match(x, data$dimnames$Terms);
	if(is.na(id)) return(NA);
	sum(data$j == id);
}


### TF IDF
tf.idf = function(x) {
	meanTF = tapply(x$v / row_sums(x)[x$i], x$j, mean);
	tf_idf = meanTF * log2(tm::nDocs(x) / col_sums(x > 0))
}


### Filter Terms:
filter.dtm = function(x, tf.idf, lim = 0.1) {
	x = x[, tf.idf >= lim[1]];
	x = x[row_sums(x) > 0, ];
	invisible(x);
}

###################

### Topic Modelling


# data = dtm;
model.lda = function(data, n, control = NULL, seed = NULL) {
	control = control.seed(control, seed=seed);
	topicmodels::LDA(data, k = n, control = control);
}
model.lda.fixed = function(data, n, control = NULL, seed = NULL) {
	if(is.null(control)) control = list();
	control$estimate.alpha = FALSE;
	control = control.seed(control, seed=seed);
	LDA(data, k = n, control = control);
}
model.lda.gibbs = function(data, n, iter = 1000, burn.in = 1000, seed = NULL) {
	control = list(iter = iter, burnin = burn.in, thin = 100);
	control = control.seed(control, seed = seed);
	LDA(data, k = n, method = "Gibbs", control = control);
}
model.ctm = function(data, n, tol.var = 10^-4, tol.em = 10^-3, seed = NULL) {
	control = list(var = list(tol = tol.var), em = list(tol = tol.em));
	control = control.seed(control, seed = seed);
	CTM(data, k = n, control = control);
}
# Demo: Multiple Models
model.demo = function(data, n, seed = NULL, verbose = TRUE) {
	resVEM = model.lda(data, n = n, seed = seed);
	if(verbose) cat("Finished VEM model.\n");
	fixVEM = model.lda.fixed(data, n = n, seed = seed);
	if(verbose) cat("Finished fixed VEM model.\n");
	Gibbs  = model.lda.gibbs(data, n = n, seed = seed);
	if(verbose) cat("Finished Gibbs model.\n");
	resCTM = model.ctm(data, n = n, seed = seed);
	if(verbose) cat("Finished CTM model.\n");
	#
	tmp = list(
		VEM    = resVEM,
		fixVEM = fixVEM,
		Gibbs  = Gibbs,
		CTM    = resCTM );
	return(tmp);
}
model.byType = function(n, dtm,
		type = c("VEM", "fixVEM", "Gibbs", "CTM", "All"),
		SEED = NULL) {
	type = match.arg(type);
	if(type == 'VEM') {
		list(VEM = model.lda(dtm, n = n, seed = SEED));
	} else if(type == 'fixVEM') {
		list(fixVEM = model.lda.fixed(dtm, n = n, seed = SEED));
	} else if(type == 'Gibbs') {
		list(Gibbs = model.lda.Gibbs(dtm, n = n, seed = SEED));
	} else if(type == 'CTM') {
		list(CTM = model.CTM(dtm, n = n, seed = SEED));
	} else {
		model.demo(tmp.dtm, n = n, seed = SEED);
	}
}

# Helper:
control.seed = function(x = NULL, seed = NULL) {
	if(! is.null(seed)) {
		if(is.null(x)) {
			x = list(seed = seed);
		} else {
			x$seed = seed;
		}
	}
	return(x);
}

##########################

### Analysis

#' @export
topicSpread = function(x) {
	spread = apply(posterior(x)$topics, 1, function(z) - sum(z * log(z)));
	mean(spread);
}

#' @export
analyseTerm = function(x, data) {
	lst = table(topic.term(x, data=data));
	list(Total = sum(lst), Topics = lst);
}

#' @export
topic.term = function(x, data) {
	id = match(x, data@terms);
	if(is.na(id)) return(NA);
	idDoc = data@wordassignments$i[data@wordassignments$j == id];
	topics(data, 1)[idDoc];
}


### Topic Separation:
# Diff between Top 2 topics:
# Note: gamma = proportion of each topic;
#' @export
diffTopics = function(x) {
	t2 = topics(x, 2);
	nTop = ncol(x@gamma);
	nDoc = nrow(x@gamma);
	id1 = nDoc * (t2[1,] - 1) + seq(nDoc);
	id2 = nDoc * (t2[2,] - 1) + seq(nDoc);
	gmd = abs(x@gamma[id1] - x@gamma[id2]);
	return(gmd);
}

