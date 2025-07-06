


### Classic n-Grams

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
}
