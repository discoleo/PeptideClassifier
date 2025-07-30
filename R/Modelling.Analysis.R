

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

# x    = Topic Model;
# data = Original Data;
summary.topics = function(x, data) {
	idTopic = topics(x, 1);
	idDocs  = as.integer(x@documents);
	dfDocs  = data;
	dfDocs$Topic = 0;
	dfDocs$Topic[idDocs] = idTopic;
	dfDocs = dfDocs[dfDocs$Topic > 0, ];
	LQ3All = quantile(dfDocs$Len, 0.75);
	tbl = tapply(dfDocs, dfDocs$Topic, function(x) {
		len = range(x$Len);
		LQ2 = median(x$Len);
		LG3 = sum(x$Len >= LQ3All);
		data.frame(Topic = x$Topic[1], N = nrow(x),
			Charge    = mean(x$Charge),
			ChargedAA = mean(x$ChargesN),
			LMin = len[1], LM = LQ2, LMax = len[2], nG3 = LG3);
	}, simplify = FALSE);
	tbl = do.call(rbind, tbl);
	return(tbl);
}


### Simple Summary: Frequencies
# x = List of TM;
table.topics = function(x) {
	idTopic   = topics(x[[1]], 1);
	tblTopics = table(idTopic);
	tblTopics = as.data.frame(tblTopics);
	len = length(x);
	if(len > 1) {
		# Compact table
		# Note: direct comparison NOT meaningful;
		nT = x[[1]]@k;
		for(id in seq(2, len)) {
			idTopic = topics(x[[id]], 1);
			tblTi   = table(idTopic);
			lenTi   = length(tblTi);
			if(lenTi < nT) tblTi = c(tblTi, rep(0, nT - lenTi));
			tblTopics = cbind(tblTopics, as.numeric(tblTi));
		}
		names(tblTopics) = c("Topic", names(x));
	}
	return(tblTopics);
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
	idT  = topics(x, 2);
	nTop = ncol(x@gamma);
	nDoc = nrow(x@gamma);
	id1 = nDoc * (idT[1,] - 1) + seq(nDoc);
	id2 = nDoc * (idT[2,] - 1) + seq(nDoc);
	gmd = abs(x@gamma[id1] - x@gamma[id2]);
	return(gmd);
}
