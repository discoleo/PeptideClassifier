
### IO

### Any File
read.pp = function(file, path = "", sep = ",", sep.fasta = "[|]") {
	if(nchar(path) > 0) {
		file = paste0(path, "/", file);
	}
	if(length(file) > 1) {
		tmp = lapply(file, read.pp, path = path,
			sep = sep, sep.fasta = sep.fasta);
		tmp = do.call(rbind, tmp);
		return(tmp);
	}
	# Read File:
	isCsv   = grepl("(?i)\\.csv$", file);
	isFasta = grepl("(?i)\\.fasta$|\\.fa$", file);
	if(isCsv) {
		x = read.csv(file);
	} else if(isFasta) {
		x = read.fasta2(file, sep = sep.fasta);
	} else {
		warning("Not yet implemented!");
		return(NULL);
	}
	x = x[x$hasSeq, c("Seq", "IDP", "Len"), drop = FALSE];
	names(x)[2] = "Type"; # IDP => Type
	type  = tolower(x$Type);
	# TODO: pmatch broken?
	# type  = pmatch(type, c("positive", "negative"));
	type  = match(type, c("pos", "neg"))
	x$Type = strType[type];
	return(x);
}

strType = c("Pos", "Neg");

filter.byType = function(flt, data) {
	x = data;
	if(flt == "All") {
		# All
	} else if(flt == "Positive") {
		x = x[x$Type == strType[1], ];
	} else if(flt == "Negative") {
		x = x[x$Type == strType[2], ];
	} else {
		warning("Invalid Filter!");
	}
	return(x);
}


### Read FASTA file
# - Variant of FASTA: pipe-separated;
# Output:
# - Start = first line with AA;
# - End   = last line with AA;
# Side-Effects:
# - print.lines: Which lines of AA-seq to print; can be logical;
read.fasta2 = function(file, path = "", token = 3, sep = "[|]") {
	fn = file;
	if(nchar(path) > 0) fn = paste0(path, "/", file);
	x  = readLines(fn);
	# Names:
	id = which(grepl("^>", x));
	sNms = x[id];
	sPr  = split.names.fasta2(sNms, token = token, sep = sep);
	lst  = data.frame(Seq = "", IDP = sPr, Start = id + 1, End = c(id[-1] - 1, length(x)));
	# AA-Seq:
	lst$hasSeq = (lst$End - lst$Start >= 0);
	# Seq-Length
	nposPP = lst[lst$hasSeq, c("Start", "End")];
	nr = nrow(nposPP);
	if(nr == 0) {
		lst$Len = 0;
		return(lst);
	}
	pp = sapply(seq(nr), function(id) {
		tmp = nposPP[id, ];
		pp  = paste0(trim(x[seq(tmp$Start, tmp$End)]), collapse = "");
	})
	lst$Seq[lst$hasSeq] = pp;
	lst$Len = nchar(lst$Seq);
	return(lst)
}

# x = Names;
split.names.fasta = function(x) {
	extract.reg(x, "^>(?:sp|tr)[|][^| ]*+", offset = c(4,0));
}

# Variant: Pipe-separator
split.names.fasta2 = function(x, token = 3, sep = "[|]",
		collapse = ".") {
	nms = strsplit(x, sep);
	if(length(token) == 1) {
		nms = sapply(nms, function(x) { x[token]; });
		return(nms);
	}
	nms = sapply(nms, function(x) {
		paste0(x[token], collapse = collapse);
	});
	return(nms);
}

### Extract String

# Regex:
extract.reg = function(x, pattern, offset = c(0,0), perl = TRUE) {
	npos   = regexpr(pattern, x, perl=perl);
	hasStr = npos > 0;
	sTk  = rep("", length(x));
	if( ! any(hasStr)) return(sTk);
	posS = npos[hasStr];
	posE = attr(npos, "match.length")[hasStr];
	posE = posE + posS - 1;
	#
	sTk[hasStr] = substr(x[hasStr], posS + offset[1], posE + offset[2]);
	return(sTk);
}

trim = function(x) {
	sub("^[ \t\n\r]+", "",
		sub("[ \t\n\r]+$", "", x));
}
