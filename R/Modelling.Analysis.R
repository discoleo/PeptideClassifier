

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
	tbl = tapply(dfDocs, dfDocs$Topic, function(x) {
		len = range(x$Len);
		data.frame(Topic = x$Topic[1], N = nrow(x),
			Charge    = mean(x$Charge),
			ChargedAA = mean(x$ChargesN),
			LMin = len[1], LMax = len[2]);
	}, simplify = FALSE);
	tbl = do.call(rbind, tbl);
	return(tbl);
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
