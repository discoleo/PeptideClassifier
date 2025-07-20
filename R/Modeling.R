#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## Modelling: Data Preparation
##
## draft v.0.1e


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
ngrams.select = function(x, type, prefix = c("_", "+", "="),
		breaks = c(5)) {
	types = c("2", "2u", "3", "3u", "4", "4u", "Len",
		"ch3.tot", "ch3.aa", "ch4.tot", "ch4.aa", "ch5.tot", "ch5.aa");
	iType = match(type, types);
	isNA  = is.na(iType);
	if(any(isNA)) warning("Some types are NOT yet implemented: ", type[isNA]);
	type  = types[iType[! isNA]];
	# n-Grams:
	xall  = NULL;
	isGrS = "2" %in% type; isGrU = "2u" %in% type;
	if(isGrS || isGrU) {
		tmp.gr = ngrams(x, n = 2);
		if(isGrS) xall = merge.list(xall, tmp.gr);
		if(isGrU) xall = merge.list(xall,
			as.ngram.undirected(tmp.gr, prefix = prefix[1]));
	}
	isGrS = "3" %in% type; isGrU = "3u" %in% type;
	if(isGrS || isGrU) {
		tmp.gr = ngrams(x, n = 3);
		if(isGrS) xall = merge.list(xall, tmp.gr);
		if(isGrU) xall = merge.list(xall,
			as.ngram.undirected(tmp.gr, prefix = prefix[1]));
	}
	isGrS = "4" %in% type; isGrU = "4u" %in% type;
	if(isGrS || isGrU) {
		tmp.gr = ngrams(x, n = 3);
		if(isGrS) xall = merge.list(xall, tmp.gr);
		if(isGrU) xall = merge.list(xall,
			as.ngram.undirected(tmp.gr, prefix = prefix[1]));
	}
	# Length of PP:
	if("Len" %in% type) {
		tmp.gr = as.list(len.pp(x, breaks = breaks[1]));
		xall   = merge.list(xall, tmp.gr);
	}
	# Charges:
	isCh = grepl("^ch", type);
	if(! any(isCh)) return(xall);
	type = type[isCh];
	# Total Charge:
	pp.charge  = as.charges(x);
	isChargeTt = c("ch3.tot", "ch4.tot", "ch5.tot") %in% type;
	idChargeTt = which(isChargeTt);
	if(length(idChargeTt) > 0) {
		n = c(3,4,5);
		for(id in idChargeTt) {
			tmp  = ngrams.charge.numeric(pp.charge, n = n[id], prefix = prefix[2]);
			xall = merge.list(xall, tmp);
		}
	}
	# Charged AA:
	isChargeAA = c("ch3.aa", "ch4.aa", "ch5.aa") %in% type;
	idChargeAA = which(isChargeAA);
	if(length(idChargeAA) > 0) {
		n = c(3,4,5);
		for(id in idChargeAA) {
			tmp  = ngrams.charged.numeric(pp.charge, n = n[id], prefix = prefix[3]);
			xall = merge.list(xall, tmp);
		}
	}
	return(xall);
}

### Simple n-Grams
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
# Prefix: enables distinction between directed & undirected n-Grams;
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
ngrams.charge = function(x, n = 4, breaks = 5, prefix = "+") {
	x   = as.charges(x);
	lst = ngrams.charge.numeric(x, n=n, breaks=breaks, prefix=prefix);
	return(lst);
}
ngrams.charge.numeric = function(x, n = 4, breaks = 5, prefix = "+") {
	if(length(x) == 0) return(x);
	lst = lapply(x, function(x) {
		nch  = length(x);
		# N & C-Terminus:
		x[1] = x[1] + 1; x[nch] = x[nch] - 1;
		if(nch <= n) {
			chAA = sum(x);
			return(chAA);
		}
		LEN = nch - n + 1;
		lst = rep(0, LEN);
		for(npos in seq(LEN)) {
			lst[npos] = sum(x[seq(npos, npos + n - 1)]);
		}
		return(lst);
	})
	# TODO: breaks;
	if(! is.null(prefix)) {
		# n-Gram Length needed to differentiate;
		prefix = paste0(n, prefix);
		lst = lapply(lst, function(x) {
			paste0(prefix, x);
		});
	}
	return(lst);
}
ngrams.charge.old = function(x, n = 4, breaks = 5, prefix = "+-") {
	x = lapply(x, function(x) {
		as.numeric(charToRaw(x)) - 64;
	});
	lst = ngrams.charge.numericChar(x, n=n, breaks=breaks, prefix=prefix);
	return(lst);
}
# [Old approach]
ngrams.charge.numericChar = function(x, n = 4, breaks = 5, prefix = "+-") {
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

### Charged AA per n-Gram
# TODO:
# - test which method is faster: charToRaw vs as.charges;
ngrams.charged = function(x, n = 4, breaks = 5, prefix = "=") {
	x = as.charges(x);
	lst = ngrams.charged.numeric(x, n=n, breaks=breaks, prefix=prefix);
	return(lst);
}
ngrams.charged.numeric = function(x, n = 4, breaks = 5, prefix = "=") {
	if(length(x) == 0) return(x);
	lst = lapply(x, function(x) {
		nch = length(x);
		if(nch <= n) {
			chAA = sum(x != 0);
			return(chAA);
		}
		LEN = nch - n + 1;
		lst = rep(0, LEN);
		for(npos in seq(LEN)) {
			lst[npos] = sum(x[seq(npos, npos + n - 1)] != 0)
		}
		return(lst);
	})
	# TODO: breaks;
	if(! is.null(prefix)) {
		# n-Gram Length needed to differentiate;
		prefix = paste0(n, prefix);
		lst = lapply(lst, function(x) {
			paste0(prefix, x);
		});
	}
	return(lst);
}

### Convert to Seq of Charges:
# Note: N/C-Terminal AA are NOT marked as extra-charged;
as.charges = function(x) {
	len = length(x);
	if(len == 0) return(x);
	sSz = nchar(x);
	tmp = lapply(seq(len), function(id) {
		szSeq = sSz[[id]];
		if(szSeq == 0) return(numeric(0));
		seqAA = x[[id]];
		seqCh = rep(0, szSeq);
		for(npos in seq(szSeq)) {
			ch = substr(seqAA, npos, npos);
			if(ch == 'D' || ch == 'E') {
				seqCh[npos] = -1;
			} else if(ch == 'H' || ch == 'K' || ch == 'R') {
				seqCh[npos] = +1;
			}
		}
		return(seqCh);
	});
	return(tmp);
}

# Helper:

merge.list = function(x, y) {
	if(is.null(x)) return(y);
	len = length(x);
	if(is.null(y)) return(x);
	if(length(y) != len) stop("Length Mismatch!");
	if(len == 0) return(x);
	tmp = lapply(seq(len), function(id) {
		c(x[[id]], y[[id]]);
	});
	return(tmp);
}

###################
###################
