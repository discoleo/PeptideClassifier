#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## Modelling: DTM


### DTM

# Generate DTM
#' @export
dtm = function(x, min.len = 1) {
	# Skip re-Tokenization: How?
	x.str = sapply(x, \(x) paste0(x, collapse = " "));
	tm::DocumentTermMatrix(x.str,
		control = list(
			# tokenize = function(x) x,
			stemming = FALSE, stopwords = FALSE, tolower = FALSE,
			wordLengths   = c(min.len, Inf),
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
# - Mean(TF.IDF) per all documents;
tf.idf = function(x) {
	meanTF = tapply(x$v / row_sums(x)[x$i], x$j, mean);
	tf_idf = meanTF * log2(tm::nDocs(x) / col_sums(x > 0))
}

### Retained Terms
# x = DTM;
terms.doc = function(x) {
	tapply(x$v, x$i, sum);
}

### Terms for document[n];
# x = DTM;
terms.byDoc = function(n, x) {
	# print(x$j[x$i == n]);
	x$dimnames[[2L]][x$j[x$i == n]];
}


### Filter Terms:
# Note:
# - does NOT seem a perfect technique;
# - Suppose: N documents each with 16 terms,
#   one word dominant in 2 clusters out of 6 clusters (of similar size);
#   => tf.idf = (log2(6) - 1) / 16 = 0.099;
filter.dtm = function(x, tf.idf, lim = 0.1) {
	x = x[, tf.idf >= lim[1]];
	idRm = which(row_sums(x) == 0);
	# Filter out Docs with NO terms left;
	if(length(idRm) > 0) {
		x = x[- idRm, ];
		attr(x, "idDocRm") = idRm;
	}
	invisible(x);
}

table.term.idf = function(x, tf.idf, lim = 0.1, as.df = TRUE) {
	isFlt = tf.idf < lim[1];
	ids = which(isFlt);
	if(length(ids) == 0) {
		if(as.df) return(data.frame(Term = character(0), Freq = numeric(0)));
		return(table(numeric(0)));
	}
	# x$j == Term id;
	tbl = table(x$j[x$j %in% ids]);
	trm = x$dimnames$Terms[which(isFlt)];
	names(tbl) = trm;
	if(as.df) {
		tbl = as.data.frame(tbl);
		names(tbl)[1] = "Term";
	}
	return(tbl);
}

table.term = function(term, x) {
	idT = match(term, x$dimnames$Terms);
	tbl = table(x$j[x$j %in% idT]);
	names(tbl) = term;
	tbl = as.data.frame(tbl);
	names(tbl)[1] = "Term";
	return(tbl);
}

# Filtered Terms:
which.term.idf = function(x, tf.idf, lim = 0.1) {
	isFlt = tf.idf < lim[1];
	x$dimnames$Terms[which(isFlt)];
}

# Filtered Docs:
which.doc.idf = function(x, tf.idf, lim = 0.1) {
	x = x[, tf.idf >= lim[1]];
	which(row_sums(x) == 0);
}


