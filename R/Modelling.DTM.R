#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## Modelling: DTM


### DTM

### Note:
# DTM in tm:
# $i = id Doc;
# $j = id Term;
# $v = Term frequency;
# $dimnames = list(Docs, Terms);


### Generate DTM
#' @export
dtm = function(x, min.len = 1) {
	# Skip re-Tokenization: How?
	# tokenize.NOP = function(x) x;
	x.str = sapply(x, \(x) paste0(x, collapse = " "));
	tm::DocumentTermMatrix(x.str,
		control = list(
			# tokenize = tokenize.NOP,
			stemming = FALSE, stopwords = FALSE, tolower = FALSE,
			wordLengths   = c(min.len, Inf),
			removeNumbers = FALSE, removePunctuation = FALSE));
}


### Dim
#' @export
dim.dtm = function(x) {
	dim(x);
}

### Document Names/Labels
#' @export
labels.dtm = function(x) {
	x$dimnames[[1]];
}

# x = DTM;
# dtmFlt = Filtered DTM;
summary.dtm = function(x, dtmFlt, tf.idf) {
	# x = dtm;
	tbl = summary(col_sums(x));
	# tbl = data.frame(as.list(tbl), check.names = FALSE);
	# Table: in 1 column;
	tbl = unclass(tbl);
	tbl = data.frame(Stat = names(tbl), DTM = tbl);
	tbl = cbind(tbl, "DTM.Filtered" = unclass(summary(col_sums(dtmFlt))));
	tbl = cbind(tbl, "TF.IDF" = unclass(summary(tf.idf)));
	tbl = cbind(tbl, "Terms"  = unclass(summary(terms.doc(x))));
	rownames(tbl) = NULL;
	# Dim:
	dim1 = dim(x);
	dim2 = dim(dtmFlt);
	df2  = data.frame(c("Docs", "Terms"), dim1, dim2,
			c(0,0), dim1);
	names(df2) = names(tbl);
	tbl = rbind(tbl, df2, make.row.names = FALSE);
	return(tbl);
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
# - TF-IDF: does NOT seem a perfect technique;
# - Suppose: N documents each with 16 terms,
#   one word dominant in 2 clusters out of 6 clusters (of similar size);
#   => tf.idf = (log2(6) - 1) / 16 = 0.099;
filter.idf.dtm = function(x, tf.idf, lim = 0.1) {
	x = x[, tf.idf >= lim[1]];
	invisible(filter.docs.dtm(x));
}

# Filter DTM by Term Frequency
filter.freq.dtm = function(x, lim = 1) {
	if(lim < 1) return(x);
	tbl = table(x$j);
	idT = as.integer(names(tbl));
	isFlt = tbl <= lim;
	x = x[, ! isFlt];
	invisible(filter.docs.dtm(x));
}
# Filter out Docs with NO terms left;
filter.docs.dtm = function(x) {
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

# Frequency: Only terms in "term"
table.term = function(term, x) {
	idT = match(term, x$dimnames$Terms);
	tbl = table(x$j[x$j %in% idT]);
	idT = as.numeric(names(tbl));
	names(tbl) = x$dimnames$Terms[idT];
	tbl = as.data.frame(tbl);
	names(tbl)[1] = "Term";
	tbl$Term = as.character(tbl$Term);
	return(tbl);
}
# Frequency: All Terms
table.terms = function(x) {
	tbl = table(x$j);
	idT = as.numeric(names(tbl));
	names(tbl) = x$dimnames$Terms[idT];
	tbl = as.data.frame(tbl);
	names(tbl)[1] = "Term";
	tbl$Term = as.character(tbl$Term);
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


