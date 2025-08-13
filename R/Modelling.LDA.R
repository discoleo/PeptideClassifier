

### Topic Modelling

### LDA Variants

# data = dtm;
model.lda = function(data, n, control = NULL, seed = NULL) {
	control = control.seed(control, seed=seed);
	topicmodels::LDA(data, k = n, control = control);
}
model.lda.fixed = function(data, n, control = NULL, seed = NULL) {
	if(is.null(control)) control = list();
	control$estimate.alpha = FALSE;
	control = control.seed(control, seed=seed);
	topicmodels::LDA(data, k = n, control = control);
}
model.lda.gibbs = function(data, n, iter = 1000, burn.in = 1000, seed = NULL) {
	control = list(iter = iter, burnin = burn.in, thin = 100);
	control = control.seed(control, seed = seed);
	topicmodels::LDA(data, k = n, method = "Gibbs", control = control);
}
### Correlated Topic Model
# tol.em = Tolerance for the EM-algorithm;
model.ctm = function(data, n, tol.var = 10^-4, tol.em = 10^-3,
		iter = 500, iter.em = 1000, seed = NULL) {
	control = list(
		var = list(tol = tol.var, iter.max = iter),
		em  = list(tol = tol.em,  iter.max = iter.em));
	control = control.seed(control, seed = seed);
	topicmodels::CTM(data, k = n, control = control);
}
# Demo: Multiple Models
model.demo = function(data, n, seed = NULL, verbose = TRUE, ...) {
	if(verbose) cat("Starting to compute the Topic Models ...\n");
	resVEM = model.lda(data, n = n, seed = seed);
	if(verbose) cat("Finished VEM model.\n");
	fixVEM = model.lda.fixed(data, n = n, seed = seed);
	if(verbose) cat("Finished fixed VEM model.\n");
	Gibbs  = model.lda.gibbs(data, n = n, seed = seed, ...);
	if(verbose) cat("Finished Gibbs model.\n");
	resCTM = model.ctm(data, n = n, seed = seed, ...);
	if(verbose) cat("Finished CTM model.\n");
	#
	tmp = list(
		VEM    = resVEM,
		fixVEM = fixVEM,
		Gibbs  = Gibbs,
		CTM    = resCTM );
	return(tmp);
}
model.byType = function(n, dtm,
		type = c("VEM", "fixVEM", "Gibbs", "CTM", "All"),
		SEED = NULL, ...) {
	type = match.arg(type);
	if(type == 'VEM') {
		list(VEM = model.lda(dtm, n = n, seed = SEED));
	} else if(type == 'fixVEM') {
		list(fixVEM = model.lda.fixed(dtm, n = n, seed = SEED));
	} else if(type == 'Gibbs') {
		list(Gibbs = model.lda.gibbs(dtm, n = n, seed = SEED, ...));
	} else if(type == 'CTM') {
		list(CTM = model.ctm(dtm, n = n, seed = SEED, ...));
	} else {
		model.demo(dtm, n = n, seed = SEED, ...);
	}
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
