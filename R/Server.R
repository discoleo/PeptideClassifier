
getServer = function(x) {
	shiny::shinyServer(server.app)
}

server.app = function(input, output, session) {
	# Global Options
	options = list(
		fltUNK    = 0.55, # Default value for Rank-Filter;
		sep = ",",        # csv Separator
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
	updateNumericInput(session, "fltUNK", value = options$fltUNK);
	
	# Dynamic variable
	values = reactiveValues(
		fullData  = NULL,   # initial Data
		dfGlData  = NULL,   # Globally filtered Data
		dfFltData = NULL,   # Data filtered in Table
		dataDTM   = NULL,   # Filtered Seq-Data for Topic Models
		dtmData   = NULL,   # DTM
		dtmFlt    = NULL,   # Filtered DTM
		tf.idf    = NULL,   # TF-IDF
		reg.Data  = options$reg.Data,
		fltUnk    = NULL,   # is set automatically
		fltType   = NULL,
		fltCols   = NULL,
		fltDTMSeq = "Positive",
		brkLen    = c(0, 9, 19, 29, 39, 100),
		# Topic Models
		nClusters   = 0,
		seedTopics  = NULL,
		tmResult    = NULL,
		NULLARG = NULL
	);
	
	
	### File Input
	observeEvent(input$file, {
		file1 = input$file;
		if (is.null(file1)) 
			return(NULL);
		#
		x = read.pp(file1$datapath, sep = options$sep);
		x = cbind(x, charge.pp(x$Seq));
		x$ChargesN = x$Charged / x$Len;
		#
		values$fullData = x;
	})
	
	### Filters
	
	observeEvent(input$fltType, {
		if(! is.null(values$fltType) && input$fltType == values$fltType) return();
		values$fltType = input$fltType;
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
		fltType = values$fltType;
		x = values$fullData;
		x = filter.byType(fltType, data = x);
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
	
	### Modeling: DTM
	
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
		fltLen = input$fltLen;
		isData = xdf$Len >= fltLen[1] & xdf$Len <= fltLen[2];
		xdt = xdf$Seq[isData];
		### Seq Data:
		values$dataDTM = xdt;
		### n-Grams:
		xgr = ngrams.demo(xdt);
		### DTM
		tmp.dtm = dtm(xgr);
		values$dtmData = tmp.dtm;
		values$dtmFlt  = tmp.dtm; # Not yet filtered;
		values$tf.idf  = tf.idf(tmp.dtm);
	})
	
	observeEvent(input$btnDTM, {
		filterSeq();
	})
	observeEvent(input$fltLen, {
		filterSeq();
	})
	observeEvent(input$fltDTMSeq, {
		values$fltDTMSeq = input$fltDTMSeq;
		filterSeq();
	})
	
	observeEvent(input$btnDTMFilter, {
		lim = input$fltTF;
		dtm = values$dtmData;
		tfIDF = values$tf.idf;
		# print(summary(tfIDF));
		values$dtmFlt = filter.dtm(dtm, tfIDF, lim = lim);
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
		tbl = cbind(tbl, "TF.Filtered" = unclass(summary(col_sums(values$dtmFlt))));
		tbl = cbind(tbl, "TF.IDF" = unclass(summary(values$tf.idf)));
		rownames(tbl) = NULL;
		# Dim:
		dim1 = dim(dtm);
		dim2 = dim(values$dtmFlt);
		df2  = data.frame(c("Docs", "Terms"), dim1, dim2, c(0,0));
		names(df2) = names(tbl);
		tbl = rbind(tbl, df2, make.row.names = FALSE);
		DT::datatable(tbl, options = list(dom = 't')) |>
			formatRound(c("DTM", "TF.Filtered", "TF.IDF"), 3);
	})
	
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
}