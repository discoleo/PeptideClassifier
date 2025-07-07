


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
	}
	# 4-Grams
	# NOT yet;
	warning("Not yet!");
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
	x.str = sapply(xpu, \(x) paste0(x, collapse = " "));
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

