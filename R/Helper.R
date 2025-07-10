#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## draft v.0.1a


### Digrams of AA
digrams = function(x, ordered = TRUE) {
	# aa = LETTERS[1:25];
	# "B" = "Asx", "Z" = "Glx", "J" = "Xle"; "U" = "Se-Cys"
	# aa = aa[ - c(2, 10)];
	aa = seq(25)[ - c(2, 10)];
	a2 = expand.grid(aa, aa);
	n2 = a2[,2] * 32 + a2[,1];
	# a2 = rbind(a2[,2], a2[,1]);
	n2 = sort(n2);
	#
	# "A" = 65; # TODO: UPPER
	y  = lapply(x, function(x) as.numeric(charToRaw(x)) - 64);
	y2 = lapply(y, function(x) {
		len = length(x);
		tmp = x[seq(len-1)] * 32 + x[seq(2, len)];
		id  = match(tmp, n2);
		return(id);
	});
	if(! ordered) {
		isU = a2[,2] > a2[,1];
		idU = seq(length(n2));
		idU[isU] = match(a2[isU, 1]*32 + a2[isU, 2], n2);
		y2 = lapply(y2, function(id) idU[id]);
	}
	invisible(y2);
}
as.digrams.char = function(x) {
	# x1 = x %/% 32;
	# x2 = x %% 32;
	# aa = paste0(aa[x1], aa[x2]);
	aa = LETTERS[1:25][ - c(2, 10)];
	# "B" = "Asx", "Z" = "Glx", "J" = "Xle";
	a2 = expand.grid(aa, aa);
	a2 = paste0(a2[,2], a2[,1]);
	aa = a2[x];
	return(aa);
}

# x = list with n-grams;
table.ngram = function(x) {
	len = length(x);
	if(len == 0) return(NULL);
	x = lapply(seq(len), function(id) {
		unique(x[[id]]);
	});
	x   = unlist(x);
	tbl = table(x);
	return(tbl);
}

####################

### Charge / Charges

charge.pp = function(x) {
	len = length(x);
	if(len == 0) {
		tmp = data.frame(
			Charge    = numeric(0),
			ChargedAA = numeric(0));
		return(tmp);
	}
	sSz = nchar(x);
	tmp = sapply(seq(len), function(id) {
		szSeq = sSz[[id]];
		if(szSeq == 0) return(c(0,0));
		seqAA = x[[id]];
		totCh = 0; nChAA = 0;
		for(npos in seq(szSeq)) {
			ch = substr(seqAA, npos, npos);
			if(ch == 'D' || ch == 'E') {
				totCh = totCh - 1; # Total Charge
				nChAA = nChAA + 1; # Number of Charges
			} else if(ch == 'H' || ch == 'K' || ch == 'R') {
				totCh = totCh + 1;
				nChAA = nChAA + 1;
			}
		}
		return(c(totCh, nChAA));
	})
	# Note:
	# - Terminal charges do NOT matter;
	# - Cyclo-peptides & C-terminal Amides are not yet covered;
	tmp = data.frame(
			Charge    = tmp[1,],
			ChargedAA = tmp[2,]);
	invisible(tmp);
}

