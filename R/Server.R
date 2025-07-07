
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
		reg.Data  = options$reg.Data,
		fltUnk    = NULL,   # is set automatically
		fltType   = NULL,
		fltCols   = NULL,
		brkLen    = c(0, 9, 19, 29, 39, 100),
		NULLARG = NULL
	);
	
	
	### File Input
	observeEvent(input$file, {
		file1 = input$file;
		if (is.null(file1)) 
			return(NULL);
		#
		x = read.pp(file1$datapath, sep = options$sep);
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
		flt = values$fltCols;
		if(! is.null(flt)) flt = list(searchCols = flt);
		DT::datatable(values$dfGlData, filter = 'top',
			options = option.regex(values$reg.Data, varia = flt));
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
}