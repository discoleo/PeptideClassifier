

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
