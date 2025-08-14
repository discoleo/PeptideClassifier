#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## URL: https://github.com/discoleo/PeptideClassifier
##
## Maths: Generalised Beta to 2D
##
## draft v.0.1a


### Beta Function: Generalised to 2D
math.Beta = function(input, tol, options = NULL) {
	if(missing(tol)) {
		if(is.null(options)) { tol = 1E-13; }
		else tol = options$tol;
	}
	# Warning: a little bit dangerous;
	px = eval(parse(text = input$inBeta_xPow));
	py = eval(parse(text = input$inBeta_yPow));
	nx = eval(parse(text = input$inBeta_RxPow));
	ny = eval(parse(text = input$inBeta_RyPow));
	k  = eval(parse(text = input$inBeta_RPow));
	k  = - k; # 1 / (...)^k in exact Formula;
	### Generalised Beta: Exact formula;
	# see formulas on GitHub:
	# https://github.com/discoleo/R/blob/master/Math/Integrals.Double.Radicals.R
	div = (py+1)*nx-(px+1)*ny;
	if(abs(div) <= 1E-12) {
		# TODO: Derivative;
		print("Int: Not yet implemented!");
	}
	vv  = c((px+1)/nx, (py+1)/ny, 1-k);
	if(any(vv < 0)) {
		res = gamma((px+1)/nx) * gamma(1-k) / gamma((px+1)/nx + 1-k) +
			- gamma((py+1)/ny) * gamma(1-k) / gamma((py+1)/ny + 1-k);
	} else {
		res = beta((px+1)/nx, 1-k) - beta((py+1)/ny, 1-k);
	}
	# Exact result:
	res = res / div;
	# Test: Compute empirically
	FUN.test = function(tol) {
		integrate(\(x) sapply(x, \(y) integrate(\(x)
			x^px * y^py / (1 - x^nx*y^ny)^k, 0, 1, rel.tol=tol)$value),
			0, 1, rel.tol=tol);
	}
	resNum = try(FUN.test(tol));
	if(inherits(resNum, "try-error")) {
		# TODO: something;
	} else {
		resNum = paste0(c("Value = ", "Error = "),
			c(resNum$value, resNum$abs.error), ";");
	}
	# Result
	lst = list(Exact = res, Numerical = resNum);
	return(lst);
}

