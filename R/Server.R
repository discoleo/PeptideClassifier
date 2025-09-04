#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## URL: https://github.com/discoleo/PeptideClassifier
##
## draft v.0.2c


getServer = function(x) {
	shiny::shinyServer(server.app)
}

server.app = function(input, output, session) {
	# Global Options
	options = list(
		fltGlobalLen    = c(6, 60), # Default value for Global Length;
		sep = ",",         # csv Separator
		sepCsvCor = ",",   # csv Separator: Cor-Matrix;
		# Regex & Other Options:
		reg.Data  = TRUE,  # Regex for Data-Table
		reg.PP    = TRUE,  # Regex for Epitopes-Table
		highlight = TRUE,  # Highlight Search Term
		# Topic Models:
		seedTM    = NULL,  # Seed: 123
		### Trees:
		# Solitary Leaf vs L2-Leaf;
		colTreeLeaves = c(L1 = "#FF886490", L2 = "#6432FFB8"),
		# TODO
		NULL
	);
	
	const = list(
		Info	= "Enter ...",
		NULL
	);
	
	asChar = function(x) {
		if(is.null(x)) return("");
		return(as.character(x));
	}
	asCollapsed = function(x, collapse = "; ") {
		if(is.null(x)) return("");
		x = x[! is.na(x)];
		if(length(x) == 0) return("");
		x = paste0(x, collapse = collapse);
		return(x);
	}
	
	### Init:
	updateNumericInput(session, "fltGlobalLen",
		value = options$fltGlobalLen);
	updateTextInput(session, "inTMSeed",
		value = asChar(options$seedTM));
	updateTextInput(session, "inTreeColLeaf",
		value = asCollapsed(options$colTreeLeaves));
	
	# Dynamic variable
	values = reactiveValues(
		fullData  = NULL,   # initial Data
		dfGlData  = NULL,   # Globally filtered Data
		dfFltData = NULL,   # Data filtered in Table
		dfDTMData = NULL,   # Filtered Data used for DTM
		dataDTM   = NULL,   # Filtered Seq-Data for Topic Models
		dtmData   = NULL,   # DTM (Un-Filtered)
		dtmFlt    = NULL,   # Filtered DTM
		tf.idf    = NULL,   # TF-IDF
		termsAll  = NULL,   # All Terms (non-filtered)
		# Filters:
		fltGlobalLen  = NULL,   # is set automatically
		reg.Data  = options$reg.Data,
		fltType   = NULL,
		fltCols   = NULL,
		fltDTMSeq  = "Positive",
		fltDocFreq = 1, # Term Frequency in > Docs;
		# Filtering: TF-IDF
		termsFlt  = NULL,   # Terms removed following TF-IDF
		idDocRm   = NULL,   # Docs removed following TF-IDF
		isDocRm   = FALSE,  # Are there any Docs Removed?
		# Categories:
		brkLen    = c(0, 9, 19, 29, 39, 100), # PP-Length
		# Topic Models
		nClusters   = 0,
		seedTM      = NULL,  # Seed for TMs;
		tmResult    = NULL,  # the TM Models;
		summaryTM   = NULL,  # Summary of all TMs;
		idModel     = 1, # when Multiple Models
		topTopics   = 1, # Download Top Topics
		# Explore / Analyse Topics
		idTopic     = 1,
		# Hierarchical Clusters
		clustResult  = NULL,
		clustSubTree = NULL,
		colTreeLeaves  = options$colTreeLeaves,
		pruneTreeSize  = 0, # Prune Tree: Min Size of Branches;
		adjTreeHeight  = TRUE, # Remove Inversions
		# Clustering: Diagnostics
		clustResultAll = NULL, # List with All Trees
		clustBrRatios  = NULL, # Branch Ratios (all branches)
		corClusters    = NULL,
		# Maths
		optMath = list(tol = 1E-13),
		NULLARG = NULL
	);
	
	
	### File Input
	observeEvent(input$file, {
		file1 = input$file;
		if (is.null(file1)) 
			return(NULL);
		#
		x = read.pp(file1$datapath, sep = options$sep);
		x = augment.pp(x);
		#
		values$fullData = x;
	})
	
	### Filters
	
	observeEvent(input$fltType, {
		if(! is.null(values$fltType) && input$fltType == values$fltType) return();
		values$fltType = input$fltType;
		filter.df();
	})
	observeEvent(input$fltGlobalLen, {
		len = input$fltGlobalLen;
		values$fltGlobalLen = len;
		filter.df();
	})
	observeEvent(input$chkRegex, {
		isReg = input$chkRegex;
		if(values$reg.Data != isReg) {
			values$reg.Data = isReg;
		}
	})
	observeEvent(values$fullData, {
		filter.df();
	})
	observeEvent(input$tblData_rows_all, {
		filter.byTable();
	})
	filter.df = function() {
		# TODO
		if(is.null(values$fullData)) return();
		# Type: Positive vs Test PP
		fltType = values$fltType;
		x = values$fullData;
		x = filter.byType(fltType, data = x);
		# Length (Global Filter)
		len = values$fltGlobalLen;
		x = x[x$Len >= len[1] & x$Len <= len[2], ];
		#
		values$dfGlData = x;
		cat("Rows: ", nrow(x), "\n");
	}
	filter.byTable = reactive({
		# Analysis is performed on the filtered data;
		print("Filter Table:");
		id = input$tblData_rows_all;
		values$dfFltData = values$dfGlData[id, ];
	})
	
	### Regex
	
	option.regex = function(x, varia = NULL, caseInsens = TRUE) {
		opt = list(search = list(regex = x, caseInsensitive = caseInsens),
			searchHighlight = options$highlight);
		if( ! is.null(varia)) opt = c(opt, varia);
		return(opt);
	}
	
	### Tables
	
	buildTable = function(x, filter = 'top', dom = NULL) {
		# dom = 'tip', 't', 'tp', NULL;
		if(is.logical(filter)) {
			filter = if(filter) 'top' else NULL;
		}
		hasDom = ! is.null(dom);
		# DT Table:
		if(is.null(filter)) {
			if(hasDom) {
				DT::datatable(x, options = list(dom = dom));
			} else {
				DT::datatable(x);
			}
		} else {
			if(hasDom) {
				opt = option.regex(values$reg.Data, varia = list(dom = dom));
			} else {
				opt = option.regex(values$reg.Data);
			}
			DT::datatable(x, filter = filter, options = opt);
		}
	}
	
	# Table: Primary Data
	dataTable = function() {
		if(is.null(values$dfGlData)) return(NULL);
		flt = values$fltCols;
		if(! is.null(flt)) flt = list(searchCols = flt);
		DT::datatable(values$dfGlData, filter = 'top',
			options = option.regex(values$reg.Data, varia = flt)) |>
		formatRound("ChargesN", 2);
	}
	
	output$tblData = DT::renderDT(dataTable())
	
	output$tblDataSummary = DT::renderDT({
		xdf = values$dfFltData;
		if(is.null(xdf)) return(NULL);
		len = cut(xdf$Len, breaks = values$brkLen);
		tbl = as.data.frame(table(len));
		names(tbl)[1] = "Length";
		tbl = rbind(tbl, data.frame(Length = "Total", Freq = sum(tbl$Freq)));
		DT::datatable(tbl, options = list(dom = 't'));
	})
	
	output$txtDataSummary = renderText({
		if(is.null(values$dfFltData)) {
			return("Data summary: No Data loaded.");
		}
		return("Data summary:");
	})
	
	# Jump to Page:
	observeEvent(input$btnDataGoTo, {
		pp = input$txtDataPP;
		if(is.null(pp) || nchar(pp) == 0) return();
		x = values$dfFltData;
		if(is.null(x)) return();
		#
		if(grepl("^[0-9]", pp)) {
			# Page number:
			id = as.integer(pp);
			if(is.na(pp) || pp == 0) return();
			setPage(pp, 'tblData');
			return();
		}
		# Peptide Seq:
		id = which(grepl(pp, x$Seq));
		if(length(id) == 0) return();
		if(length(id) > 1) {
			if(length(id) <= 10) print(id);
			id = id[1];
		}
		# Page:
		pg = (id - 1) %/% 10 + 1; # TODO: Items per page;
		setPage(pg, 'tblData');
	})
	setPage = function(pg, tbl) {
		proxy = dataTableProxy(tbl);
		selectPage(proxy, pg);
	}
	
	### Modelling: DTM
	
	filterSeq = reactive({
		xdf = values$dfFltData;
		if(is.null(xdf)) return();
		# Filter data:
		print("Filtering Seq for DTM");
		fltType = values$fltDTMSeq;
		if(fltType == 'Positive') {
			xdf = xdf[xdf$Type == 'Pos', ];
		} else if(fltType == 'Negative') {
			xdf = xdf[xdf$Type == 'Neg', ];
		}
		# Filter Length:
		fltLen = input$fltDTMLen;
		isData = xdf$Len >= fltLen[1] & xdf$Len <= fltLen[2];
		# Filtered Data:
		xdf = xdf[isData, ];
		values$dfDTMData = xdf;
		### Seq Data:
		xdt = xdf$Seq;
		values$dataDTM = xdt;
		### n-Grams:
		tGr = input$chkNGrams; # Types of n-Grams;
		xgr = ngrams.select(xdt, tGr);
		# xgr = ngrams.demo(xdt);
		### Terms
		values$termsAll = unique(unlist(xgr));
		### DTM
		tmp.dtm = dtm(xgr);
		# Flt: Documents per Term > lim;
		tmp.dtm = filter.freq.dtm(tmp.dtm, lim = values$fltDocFreq);
		values$dtmData = tmp.dtm;
		values$dtmFlt  = tmp.dtm; # Not yet filtered;
		values$tf.idf  = tf.idf(tmp.dtm);
		# Excluded Terms:
		termsEx = setdiff(values$termsAll, tmp.dtm$dimnames[[2L]]);
		cat("Excluded Terms: ", length(termsEx), "\n");
		# print(termsEx);
	})
	
	observeEvent(input$btnDTM, {
		filterSeq();
	})
	observeEvent(input$fltDTMLen, {
		filterSeq();
	})
	observeEvent(input$fltDTMSeq, {
		values$fltDTMSeq = input$fltDTMSeq;
		filterSeq();
	})
	observeEvent(input$btnDTMInspectPP, {
		dtm = values$dtmFlt;
		if(is.null(dtm)) {
			output$txtDTM_PP = renderText("Nothing to inspect: no DTM!");
			output$txtDTM_PP_Terms = NULL;
			return();
		}
		# Doc ID:
		id = input$fltDTMDocID;
		# Reset:
		if(is.null(id) || nchar(id) == 0) {
			output$txtDTM_PP = NULL;
			output$txtDTM_PP_Terms = NULL;
			return();
		}
		isNum = grepl("^\\-?[0-9]", id);
		if(! isNum) {
			id = match(id, values$dfDTMData$Seq);
			if(is.na(id)) {
				cat("Invalid Seq!\n");
				output$txtDTM_PP_Terms = NULL;
				return();
			}
			# doc = id; # TODO: avoid duplicated code;
			# TODO: is id already == match(...)?
			id = match(as.character(id), labels.dtm(dtm));
		} else {
			id = as.integer(id);
		}
		if(is.na(id) || id < 1) {
			warning("Invalid Doc id!");
			id = 1;
		}
		dimDtm = dim.dtm(dtm);
		if(id > dimDtm[1L]) {
			warning("Invalid Doc id: beyond last doc!");
			id = dimDtm[1L];
		}
		# Docs & Terms:
		doc = as.numeric(dtm$dimnames[[1L]][id]);
		sPP = values$dfDTMData$Seq[doc];
		sTm = terms.byDoc(id, dtm);
		sTm = paste(sTm, collapse = " ");
		# AA-Sequence:
		sSeq = paste("Inspect:", sPP, collapse = " ");
		output$txtDTM_PP = renderText(sSeq);
		output$txtDTM_PP_Terms = renderText(sTm);
	})
	
	observeEvent(input$btnDTMFilter, {
		lim = input$fltTF;
		dtm = values$dtmData;
		# TF-IDF:
		tfIDF  = values$tf.idf;
		dtmFlt = filter.idf.dtm(dtm, tfIDF, lim = lim);
		values$termsFlt = which.term.idf(dtm, tfIDF, lim = lim);
		values$idDocRm = attr(dtmFlt, "idDocRm");
		values$dtmFlt  = dtmFlt;
		# Excluded Terms:
		# trmEx = setdiff(values$termsAll, dtmFlt$dimnames[[2L]]);
		# cat("Excluded terms: "); print(trmEx);
	})
	
	# Original DTM:
	output$tblDTMSummary = DT::renderDT({
		dtm = values$dtmData;
		if(is.null(dtm)) return();
		# Summary:
		tbl = summary.dtm(dtm, values$dtmFlt, values$tf.idf);
		DT::datatable(tbl, options = list(dom = 't')) |>
			formatRound(c("DTM", "DTM.Filtered", "TF.IDF", "Terms"),
				c(1,1,3,1));
	})
	
	# Terms Retained:
	output$tblDTMRetainedTerms = DT::renderDT({
		xdf = values$dtmFlt;
		if(is.null(xdf)) {
			output$txtTermsRetained = NULL;
			return();
		}
		output$txtTermsRetained = renderText("Terms retained:");
		tbl = table.terms(xdf);
		buildTable(tbl, 'top', dom = 'tip');
	})
	
	# Terms Filtered:
	# output$tblDTMFltTerms = DT::renderDT({
	output$tblDTMRemovedTerms = DT::renderDT({
		xTerms = values$termsFlt;
		xdf    = values$dtmData;
		if(is.null(xTerms) || is.null(xdf)) {
			output$txtTermsRemoved = NULL;
			return();
		}
		output$txtTermsRemoved = renderText("Terms removed:");
		tbl = table.term(xTerms, xdf);
		buildTable(tbl, 'top', dom = 'tip');
	})
	
	# Removed Docs:
	output$tblDTMRemovedDocs = DT::renderDT({
		xdf = values$dfDTMData;
		idDocs = values$idDocRm;
		if(is.null(idDocs) || is.null(xdf)) {
			values$isDocRm = FALSE;
			return();
		}
		xdf = xdf[idDocs, ];
		values$isDocRm = TRUE;
		# Table:
		buildTable(xdf, 'top', dom = 'tip') |>
		formatRound("ChargesN", 2);
	})
	# Title & Messages:
	output$txtDTM_RemovedDocsTitle = renderText({
		isDocRm = values$isDocRm;
		txt = if(isDocRm) {
			"Removed Polypeptides by TF-IDF:";
		} else {
			NULL;
		}
	})
	output$txtDTM_RemovedDocs = renderText({
		isDocRm = values$isDocRm;
		txt = if(isDocRm) {
			NULL;
		} else {
			"Filter TF-IDF: No Polypeptides removed."
		}
	})
	
	#####################
	### Clusters / Topics
	
	# TM Model: valid ID;
	getModel = function(nMod, idModel = 1) {
		idModel = if(nMod > 1) min(nMod, values$idModel)
			else idModel;
		return(idModel);
	}
	getModelTM = function() {
		tm  = values$tmResult;
		if(is.null(tm)) return(NULL);
		idModel = getModel(length(tm));
		tm  = tm[[idModel]];
		return(tm);
	}
	
	observeEvent(input$numTMClusters, {
		values$nClusters = input$numTMClusters;
	})
	
	observeEvent(input$inTMSeed, {
		x = input$inTMSeed;
		if(is.null(x) || nchar(x) == 0) {
			values$seedTM = NULL;
			return();
		}
		x = as.integer(x);
		if(x <= 0) x = NULL;
		values$seedTM = x;
	})
	
	observeEvent(input$btnModelTopics, {
		n = values$nClusters;
		if(n < 2) return();
		cat("Started TM:\n");
		res.tm = modelTopics(n, dtm = values$dtmFlt);
		cat("Finished TM.\n");
		values$tmResult = res.tm;
	})
	
	output$downloadTMSummary = downloadHandler(
		filename = function() {
			n = values$nClusters;
			type = input$fltTMType;
			paste("TMSummary.", type, "_", n, ".csv", sep = "");
		},
		content = function(file) {
			x = values$summaryTM;
			if(is.null(x)) return(NULL);
			write.csv(x, file, row.names = FALSE);
		}
	)
	
	output$downloadTMTopics = downloadHandler(
		filename = function() {
			n = values$nClusters;
			type = input$fltTMType;
			paste("TMTopics.", type, "_", n, ".csv", sep = "");
		},
		content = function(file) {
			x = values$tmResult;
			if(is.null(x)) return(NULL);
			# Top Topics:
			n = values$topTopics;
			if(n == 0) n = values$nClusters;
			xdf = docTopic(x, n = n);
			write.csv(xdf, file, row.names = FALSE);
		}
	)
	
	### Full Models:
	output$downloadTMFull = downloadHandler(
		filename = function() {
			n = values$nClusters;
			type = input$fltTMType;
			paste("TMFull.M", type, "_", n, ".rds", sep = "");
		},
		content = function(file) {
			x = values$tmResult;
			if(is.null(x)) return(NULL);
			seed = values$seedTM;
			x = list(TM = x, Seed = seed);
			saveRDS(x, file);
		}
	)
	observeEvent(input$loadTMFull, {
		ff = input$loadTMFull;
		if(is.null(ff)) 
			return(NULL);
		#
		x = readRDS(ff$datapath);
		if(is.null(x$TM)) return();
		cat("Finished Loading TMs.\n");
		# Set parameters:
		updateNumericInput (session, "numTMClusters",
			value = x$TM[[1]]@k);
		seed = x$Seed;
		if(is.null(seed)) seed = "";
		updateTextInput (session, "inTMSeed",
			value = asChar(seed));
		values$tmResult = x$TM;
	})
	
	observeEvent(input$inTMDownloadTop, {
		x = input$inTMDownloadTop;
		if(is.null(x) || nchar(x) == 0) {
			x = 1;
		} else {
			x = as.integer(x);
			if(is.na(x)) { x = 1; }
			else if(x < 0) {
				x = values$nClusters + x;
				if(x < 0) x = 0;
			}
		}
		values$topTopics = x;
	})
	
	# Compute TM:
	modelTopics = function(n, dtm) {
		type = input$fltTMType;
		SEED = values$seedTM;
		iter = as.numeric(input$numTM_Iterations);
		lst  = model.byType(n=n, dtm=dtm, type=type, SEED = SEED,
			iter = iter);
		return(lst);
	}
	
	### Explore / Analyse Topics
	
	# Topics: Basic Summary:
	output$tblTopicInfo = DT::renderDT({
		xtm = values$tmResult;
		if(is.null(xtm)) return();
		# Which Model:
		nMod  = length(xtm);
		idModel = getModel(nMod);
		dfDocs  = values$dfFltData;
		tbl = summary.topics(xtm[[idModel]], dfDocs);
		# Multiple Models:
		if(nMod > 1) {
			allT = table.topics(xtm);
			# CTM: sometimes less clusters;
			nrM = nrow(tbl); nr.diff = nrow(allT) - nrM;
			if(nr.diff > 0) tbl = rbind0(tbl, nr.diff);
			tbl  = cbind(allT, tbl[, -1]);
		}
		values$summaryTM = tbl;
		#
		DT::datatable(tbl, options = list(dom = 'tp')) |>
		formatRound(c("Charge", "ChargedAA"), c(2,2));
	})
	
	### Topic Terms
	output$tblTopicTerms = DT::renderDT({
		xtm = values$tmResult;
		if(is.null(xtm)) return();
		# nClusters has changed: wait for the TM;
		nClusters = values$nClusters;
		if(xtm[[1]]@k != nClusters) return();
		#
		idModel = getModel(length(xtm));
		top = input$numTermsTM;
		termsT = terms(xtm[[idModel]], top);
		termsT = termsT[, seq(nClusters)];
		DT::datatable(termsT, options = list(dom = 'tp'));
	})
	
	# Explore Model: Multiple Models
	observeEvent(input$inModelID, {
		values$idModel = as.integer(input$inModelID);
	})
	
	# Explore Specific Topic
	
	observeEvent(input$inTopicID, {
		id = round(as.numeric(input$inTopicID));
		if(is.na(id) || id <= 0 || id > values$nClusters) {
			warning("Invalid Topic id!");
			return();
		}
		values$idTopic = id;
	})
	
	# PP belonging to specific Topic / Cluster
	output$tblPPTopic = DT::renderDT({
		xdf = values$dfDTMData;
		xtm = values$tmResult;
		if(is.null(xdf) || is.null(xtm)) return();
		# Note: Docs could be extracted from the TM-Result as well;
		idDocRm = values$idDocRm;
		if(! is.null(idDocRm)) xdf = xdf[- idDocRm, ];
		id0 = values$idTopic;
		idModel = getModel(length(xtm));
		idTopic = topics(xtm[[idModel]], 1); # Only Top Topic;
		tblPP   = xdf[idTopic == id0, ];
		tblPP$nCov = coverQ(xtm[[idModel]], prob = 0.5, filter = idTopic == id0);
		DT::datatable(tblPP, filter = 'top',
			options = option.regex(values$reg.Data, varia = list(dom = "tip"))) |>
		formatRound("ChargesN", 2);
	})
	
	### Hierarchical Clustering
	
	observeEvent(input$btnTreeBuild, {
		typeSource = input$fltTreeDataSource;
		if(typeSource == "DTM") {
			# Based on DTM-Data!
			# - which is a filtered version of values$dfFltData;
			dtm = values$dtmFlt;
			if(is.null(dtm)) return();
		} else if(typeSource == "TM") {
			tm  = getModelTM();
			if(is.null(tm)) return();
			dtm = tm@gamma;
			rownames(dtm) = tm@documents;
			print(str(dtm));
		} else if(typeSource == "Tn") {
			tm  = getModelTM();
			if(is.null(tm)) return();
			id  = values$idTopic;
			dtm = topics.byMainTopic(id, tm);
		} else if (typeSource == "DTM.Tn") {
			tm  = getModelTM();
			if(is.null(tm)) return();
			id  = values$idTopic;
			idD = docs.byTopic(id, tm);
			dtm = values$dtmFlt;
			dtm = dtm[idD, ];
		} else if(typeSource == "Raw") {
			# dtm = values$dfDTMData; # Actual/Raw Data
			# TODO
		}
		type   = input$fltTreeType;
		distM  = dist(dtm);
		# Note: attr(distM, "labels") == NULL!
		hClust = hclust(distM, method = type);
		values$clustResult = hClust;
		print(str(hClust));
	})
	
	# Prune Tree:
	observeEvent(input$btnTreePruneSize, {
		size = input$inTreeMinSize;
		if(is.null(size) || nchar(size) == 0) {
			values$pruneTreeSize = 0;
			return();
		}
		size = as.integer(size);
		if(is.na(size)) size = 0;
		values$pruneTreeSize = size;
	})
	
	output$tblClusters = DT::renderDT({
		return(NULL);
	})
	output$imgTree = renderPlot({
		hClust = values$clustResult;
		if(is.null(hClust)) return();
		# Remove Inversions:
		if(values$adjTreeHeight) {
			hClust = adjust.height.hclust(hClust);
		}
		size = values$pruneTreeSize;
		if(size != 0) {
			hClust = collapse.tree(size, tree = hClust);
		}
		# Plot:
		orientH   = input$fltTreePlotOrientation;
		colLeaves = values$colTreeLeaves;
		if(! is.null(colLeaves) && ! all(is.na(colLeaves))) {
			hClust = as.dendrogramLeaf(hClust, col = colLeaves);
		} else {
			hClust = as.dendrogram(hClust);
		}
		plot(hClust, horiz = orientH);
	})
	
	### Extract SubTree:
	observeEvent(input$btnSubtree, {
		resetST = function() values$clustSubTree = NULL;
		x = values$clustResult;
		if(is.null(x)) { resetST(); return(); }
		# Include Leaf:
		strN = input$txtSubTree_Node;
		node = as.integer(strN);
		if(is.na(node)) {
			node = which(grepl(strN, values$dfDTMData$Seq));
			if(length(node) == 0) {
				cat("Found 0 sequences!\n");
				node = 0;
			} else {
				# TODO
				if(length(node) > 1) {
					cat("Found multiple sequences!\n");
					node = node[1];
				}
				# dfDTMData => DTM => FilteredDTM => Clustering
				# Note: matching is done (again) below;
				# node = which(values$clustResult$labels == node);
				# node = as.integer(node);
				cat("Node: ", node, "\n");
			}
		}
		if(node == 0) {
			resetST(); return();
		}
		iNode = node;
		if(iNode > 0) {
			# Match Seq-Label:
			iNode = match(iNode, x$labels);
			if(is.na(iNode)) { resetST(); return(); }
		}
		# Size
		size = input$txtSubTree_Size;
		size = as.integer(size);
		if(is.na(size) || size == 0) { resetST(); return(); }
		subT = subtree.nc(iNode, n = size, x);
		attr(subT, "N0") = node;
		values$clustSubTree = subT;
	})
	
	# Plot SubTree:
	output$imgSubTree = renderPlot({
		x = values$clustSubTree;
		if(is.null(x)) return();
		# Plot SubTree:
		pos = plot.subtree(x);
		print(pos);
	});
	
	readFromTree = function(add.id = FALSE) {
		x = values$clustSubTree;
		if(is.null(x)) return(NULL);
		ids = as.integer(x$labels)[x$order];
		# TODO: check if Data-source correct!
		doc = labels.dtm(values$dtmFlt)[ids];
		doc = as.integer(doc);
		dfx = values$dfDTMData;
		dfx = dfx[doc, ];
		if(add.id) dfx = cbind(ID = ids, dfx);
		return(dfx);
	}
	
	# SubTree: PP Details
	observeEvent(input$btnSubT_Details, {
		x = readFromTree(add.id = TRUE);
		if(is.null(x)) return();
		output$tblSubTree = DT::renderDT({
			DT::datatable(x, filter = 'top',
				options = option.regex(values$reg.Data, varia = list(dom = "tip"))) |>
			formatRound("ChargesN", 2);
		});
	})
	
	# SubTree: Save Data
	output$downloadSubT_Data = downloadHandler(
		filename = function() {
			node = input$txtSubTree_Node;
			paste("Data.SubTree.P", node, ".csv", sep = "");
		},
		content = function(file) {
			x = readFromTree();
			if(is.null(x)) return();
			write.csv(x, file, row.names = FALSE);
		}
	)
	# Tree: Save full Tree
	output$downloadTree = downloadHandler(
		filename = function() {
			type = input$fltTreeType;
			srcT = input$fltTreeDataSource; srcSep = "_";
			if(srcT == "DTM") { srcT = ""; srcSep = ""; }
			paste("Tree.Full.", srcT, srcSep,
				"M_", type, ".rds", sep = "");
		},
		content = function(file) {
			x = values$clustResult;
			if(is.null(x)) return();
			saveRDS(x, file = file);
		}
	)
	observeEvent(input$loadTree, {
		ff = input$loadTree;
		if(is.null(ff)) 
			return(NULL);
		#
		x = readRDS(ff$datapath);
		values$clustResult = x;
	})
	
	# Messages:
	output$txtTreeWarn = renderText({
		hasDTM = ! is.null(values$dfDTMData);
		msg = if(hasDTM) {
			"";
		} else {
			"No DTM!";
		}
		return(msg);
	})
	output$txtTreeInfo = renderText({
		hasDTM = ! is.null(values$dfDTMData);
		msg = if(hasDTM) {
			paste("Press the *SubTree* button",
				"to extract a sub-tree containing at least n Leaves",
				"and the node specified in the Node-field.");
		} else {
			paste("Note: A DocumentTermMatrix needs to be generated first:",
				"see the DTM tab.");
		}
		return(msg);
	})
	
	### Colours:
	observeEvent(input$btnTreeColLeaves, {
		val = input$inTreeColLeaf;
		if(is.null(val) || nchar(val) == 0) {
			values$colTreeLeaves = c(NA, NA);
			return();
		}
		val = strsplit(val, "[; ]+");
		val = val[[1]];
		LEN = length(val);
		if(LEN == 0) {
			values$colTreeLeaves = c(NA, NA);
			return();
		}
		isNum = grepl("^[0-9]+$", val);
		if(all(isNum)) val = as.integer(val);
		if(length(val) == 1) val = c(L1 = val, L2 = NA);
		values$colTreeLeaves = val;
	})
	
	### Clustering: Diagnostics
	
	clusterAll = function(xdf, verbose = TRUE) {
		# Clustering: All Methods
		d = dist(xdf);
		if(verbose) cat("Finished Dist; starting Clustering!\n");
		clustMethods = c("ward.D", "single", "complete", "average",
			"mcquitty", "median", "centroid", "ward.D2");
		lstClust = list();
		for(i in seq_along(clustMethods)) {
			clustX = hclust(d, method = clustMethods[i]);
			lstClust[[i]] = clustX;
		}
		if(verbose) cat("Finished Clustering!\n");
		# Add names:
		names(lstClust) = clustMethods;
		lstClust = as.TreeList(lstClust);
		return(lstClust);
	}
	
	### Cluster: 1 Method
	observeEvent(input$btnDx_BranchSummary_1, {
		x = values$clustResult;
		if(is.null(x)) return();
		x = list(x);
		names(x) = x[[1]]$method;
		values$clustResultAll = x;
		values$clustBrRatios = summary.branchRatiosTn(x);
	})
	### Clusters: All Methods
	observeEvent(input$btnDx_BranchSummary, {
		xdf = values$dtmFlt; # Based on DTM;
		if(is.null(xdf)) return();
		if(is.null(values$clustResultAll)) {
			x = clusterAll(xdf);
			values$clustResultAll = x;
		} else {
			x = values$clustResultAll;
		}
		# Branch Ratios:
		values$clustBrRatios = summary.branchRatiosTn(x);
	})
	
	# Trees: Save All Trees
	output$downloadTrees_All = downloadHandler(
		filename = function() {
			paste("Tree.All.rds", sep = "");
		},
		content = function(file) {
			x = values$clustResultAll;
			if(is.null(x)) return();
			saveRDS(x, file = file);
		}
	)
	# Trees: Load All Trees
	observeEvent(input$loadTrees_Dx, {
		ff = input$loadTrees_Dx;
		if(is.null(ff)) 
			return(NULL);
		#
		x = readRDS(ff$datapath);
		x = as.TreeList(x);
		values$clustResultAll = x;
		values$clustBrRatios  = summary.branchRatiosTn(x);
	})
	
	output$tblDx_BranchSummary = DT::renderDT({
		x = values$clustBrRatios;
		if(is.null(x) || nrow(x) == 0) return();
		DT::datatable(x, filter = 'top',
			options = option.regex(values$reg.Data, varia = list(dom = "tip"))) |>
		formatRound(names(x)[-1], 2);
	})
	
	### Plot: Branch Ratios
	output$imgDx_HistBranchRatios = renderPlot({
		x = values$clustResultAll;
		if(is.null(x) || length(x) == 0) return();
		LEN  = length(x);
		nCol = ceiling(LEN / 2);
		nRow = if(LEN == 1) 1 else 2;
		par.old = par(mfrow = c(nRow, nCol));
		on.exit(par(par.old));
		#
		for(id in seq(LEN)) {
			rBr = branch.ratios(x[[id]]);
			hist(rBr, main = names(x)[[id]], xlab = "Branch Ratio");
		}
	})
	
	###############
	
	### Clustering: Correlations
	
	observeEvent(input$btnTreeCor, {
		# xdf = values$dfDTMData;
		xdf = values$dtmFlt; # Based on DTM;
		if(is.null(xdf)) return();
		# Required Packages:
		pkgs = c("dendextend", "corrplot");
		if(!(require(pkgs[1], character.only = TRUE) &
			require(pkgs[2], character.only = TRUE))) {
			return();
		}
		# Clustering: All Methods
		lst = clusterAll(xdf);
		lstClust = dendlist();
		LEN = length(lst); print(str(lst))
		if(LEN == 0) return();
		for(i in seq(LEN)) {
			lstClust = dendlist(lstClust, as.dendrogram(lst[[i]]));
		}
		names(lstClust) = names(lst);
		cat("Finished Clustering!\n");
		cat("Starting cor.dendlist: this takes time!\n");
		### Correlations:
		corClusters = cor.dendlist(lstClust);
		cat("Finished cor.dendlist!\n");
		values$corClusters = corClusters;
	})
	
	# Plot: Correlations between Trees
	output$imgTreeCor = renderPlot({
		x = values$corClusters;
		if(is.null(x)) return();
		if(! require("corrplot", character.only = TRUE)) {
			return();
		}
		typeCor  = input$fltDxCorTypes;
		orderCor = input$fltDxCorOrder;
		# Plot:
		corrplot::corrplot(x, method = typeCor, type = "lower",
			order = orderCor);
	})
	
	# Correlation between Trees: Save Values
	output$downloadTreeCor = downloadHandler(
		filename = function() {
			"Tree.Corr.All.csv";
		},
		content = function(file) {
			x = values$corClusters;
			if(is.null(x)) return();
			write.csv(x, file, row.names = FALSE);
		}
	)
	observeEvent(input$loadTreeCor, {
		ff = input$loadTreeCor;
		if(is.null(ff)) 
			return(NULL);
		#
		x = read.csv(ff$datapath, sep = options$sepCsvCor);
		nms = names(x);
		x   = as.matrix(x);
		dimnames(x) = list(nms, nms);
		values$corClusters = x;
	})
	
	# Messages:
	output$txtTreeDx_Warn = renderText({
		hasDTM = ! is.null(values$dfDTMData);
		msg = if(hasDTM) {
			"Warning: Takes quite some time to compute!";
		} else {
			"No DTM!";
		}
		return(msg);
	})
	output$txtTreeDx_Info = renderText({
		hasDTM = ! is.null(values$dfDTMData);
		msg = if(hasDTM) {
			paste("Plot the correlations between the trees",
				"generated using the various hierarchical clustering methods.",
				"Either press the *Correlation* button to compute these values;",
				"or load a previously saved matrix of correlations."
				);
		} else {
			paste("Note: A DocumentTermMatrix needs to be generated first:",
				"see the DTM tab.");
		}
		return(msg);
	})
	
	### Generalised Beta Distribution
	
	observeEvent(input$btnBetaCompute, {
		res = math.Beta(input, options = values$optMath);
		output$txtBetaResult  = renderText(res$Exact);
		output$txtBetaNumeric = renderText(res$Numerical);
	})
}