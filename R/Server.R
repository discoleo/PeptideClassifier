
getServer = function(x) {
	shiny::shinyServer(server.app)
}

server.app = function(input, output, session) {
	# Global Options
	options = list(
		fltGlobalLen    = c(6, 60), # Default value for Global Length;
		sep = ",",         # csv Separator
		# Regex & Other Options:
		reg.Data  = TRUE,  # Regex for Data-Table
		reg.PP    = TRUE,  # Regex for Epitopes-Table
		highlight = TRUE,  # Highlight Search Term
		# TODO
		# Protein Graph:
		col.Pr    = "#FF0032A0",
		border.Pr = "#640000A0",
		NULL
	);
	
	const = list(
		Info	= "Enter ...",
		NULL
	);
	
	### Init:
	updateNumericInput(session, "fltGlobalLen",
		value = options$fltGlobalLen);
	
	# Dynamic variable
	values = reactiveValues(
		fullData  = NULL,   # initial Data
		dfGlData  = NULL,   # Globally filtered Data
		dfFltData = NULL,   # Data filtered in Table
		dfDTMData = NULL,   # Filtered Data used for DTM
		dataDTM   = NULL,   # Filtered Seq-Data for Topic Models
		dtmData   = NULL,   # DTM
		dtmFlt    = NULL,   # Filtered DTM
		tf.idf    = NULL,   # TF-IDF
		# Filters:
		fltGlobalLen  = NULL,   # is set automatically
		reg.Data  = options$reg.Data,
		fltType   = NULL,
		fltCols   = NULL,
		fltDTMSeq = "Positive",
		# Filtering: TF-IDF
		termsFlt  = NULL,   # Terms removed following TF-IDF
		idDocRm   = NULL,   # Docs removed following TF-IDF
		isDocRm   = FALSE,  # Are there any Docs Removed?
		# Categories:
		brkLen    = c(0, 9, 19, 29, 39, 100), # PP-Length
		# Topic Models
		nClusters   = 0,
		seedTopics  = NULL,
		tmResult    = NULL,
		# Explore / Analyse Topics
		idTopic     = 1,
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
	
	# Data
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
		nGr = input$chkNGrams;
		xgr = ngrams.select(xdt, nGr);
		# xgr = ngrams.demo(xdt);
		### DTM
		tmp.dtm = dtm(xgr);
		values$dtmData = tmp.dtm;
		values$dtmFlt  = tmp.dtm; # Not yet filtered;
		values$tf.idf  = tf.idf(tmp.dtm);
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
		id = as.integer(id);
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
		dtmFlt = filter.dtm(dtm, tfIDF, lim = lim);
		values$termsFlt = which.term.idf(dtm, tfIDF, lim = lim);
		values$idDocRm = attr(dtmFlt, "idDocRm");
		values$dtmFlt  = dtmFlt;
	})
	
	# Original DTM:
	output$tblDTMSummary = DT::renderDT({
		dtm = values$dtmData;
		if(is.null(dtm)) return();
		tbl = summary(col_sums(dtm));
		# tbl = data.frame(as.list(tbl), check.names = FALSE);
		# Table: in 1 column;
		tbl = unclass(tbl);
		tbl = data.frame(Stat = names(tbl), DTM = tbl);
		tbl = cbind(tbl, "DTM.Filtered" = unclass(summary(col_sums(values$dtmFlt))));
		tbl = cbind(tbl, "TF.IDF" = unclass(summary(values$tf.idf)));
		tbl = cbind(tbl, "Terms"  = unclass(summary(terms.doc(dtm))));
		rownames(tbl) = NULL;
		# Dim:
		dim1 = dim(dtm);
		dim2 = dim(values$dtmFlt);
		df2  = data.frame(c("Docs", "Terms"), dim1, dim2,
			c(0,0), dim1);
		names(df2) = names(tbl);
		tbl = rbind(tbl, df2, make.row.names = FALSE);
		DT::datatable(tbl, options = list(dom = 't')) |>
			formatRound(c("DTM", "DTM.Filtered", "TF.IDF", "Terms"),
				c(1,1,3,1));
	})
	
	# Filtered Terms:
	output$tblDTMFltTerms = DT::renderDT({
		xTerms = values$termsFlt;
		xdf    = values$dtmData;
		if(is.null(xTerms) || is.null(xdf)) return();
		tbl = table.term(xTerms, xdf);
		DT::datatable(tbl, options = list(dom = 'tip'));
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
		DT::datatable(xdf, filter = 'top',
			options = option.regex(values$reg.Data, varia = list(dom = "tip"))) |>
		formatRound("ChargesN", 2);
	})
	output$txtDTM_RemovedDocs = renderText({
		isDocRm = values$isDocRm;
		txt = if(isDocRm) {
			"Removed Polypeptides by TF-IDF:";
		} else {
			"Filter TF-IDF: No Polypeptides removed."
		}
	})
	
	#####################
	### Clusters / Topics
	
	observeEvent(input$numClusters, {
		values$nClusters = input$numClusters;
	})
	
	observeEvent(input$btnModelTopics, {
		n = values$nClusters;
		if(n < 2) return();
		res.tm = modelTopics(n, dtm = values$dtmFlt);
		values$tmResult = res.tm;
	})
	
	modelTopics = function(n, dtm) {
		type = input$fltTMType;
		SEED = values$seedTopics;
		lst  = model.byType(n=n, dtm=dtm, type=type, SEED = SEED);
		return(lst);
	}
	
	### Explore / Analyse Topics
	
	# Topics:
	output$tblTopics = DT::renderDT({
		res.tm = values$tmResult;
		if(is.null(res.tm)) return();
		idTopic   = topics(res.tm[[1]], 1);
		tblTopics = table(idTopic);
		tblTopics = as.data.frame(tblTopics);
		len = length(res.tm);
		if(len > 1) {
			# Compact table
			# Note: direct comparison NOT meaningful;
			for(id in seq(2, len)) {
				idTopic = topics(res.tm[[id]], 1);
				tblTi   = table(idTopic);
				tblTopics = cbind(tblTopics, as.numeric(tblTi));
				# names(tblTopics)[id + 1] = names(res.tm)[id];
			}
			names(tblTopics) = c("idTopic", names(res.tm));
		}
		print(str(tblTopics))
		#
		DT::datatable(tblTopics, options = list(dom = 'tp'));
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
		idTopic = topics(xtm[[1]], 1); # Only Top Topic;
		tblPP   = xdf[idTopic == id0, ];
		DT::datatable(tblPP, filter = 'top',
			options = option.regex(values$reg.Data, varia = list(dom = "tip"))) |>
		formatRound("ChargesN", 2);
	})
}