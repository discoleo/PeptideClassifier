
### DTM Tab

panelDTM = function() {
	sidebarLayout(
		sidebarPanel(
			selectInput("fltDTMSeq",
				label = "Filter sequences:", selected = "Positive",
				choices = list("All" = "All", "Positive" = "Positive",
					"Negative" = "Negative") ),
			sliderInput(inputId = "fltDTMLen", label = "Length",
				value = c(6, 40), min = 0, max = 100, step = 1),
			sliderInput(inputId = "fltTF", label = "TF-IDF",
				# more conservative value:
				value = 0.0125, min = 0, max = 0.5, step = 0.0125),
			actionButton("btnDTM", "Build DTM"),
			actionButton("btnDTMFilter", "Filter DTM"),
			actionButton("btnDTMInspectPP", "Inspect"),
			# n-Grams:
			fluidRow(HTML("&nbsp;")),
			checkboxNGrams(),
			fluidRow(
			column(6,
				# Inspect: Doc ID
				textInput(inputId = "fltDTMDocID", label = "Inspect PP",
					value = "1", width = 150) ),
			column(6,
				# Show Removed Docs:
				checkboxInput(inputId = "fltDTMDocRem", label = "Show Removed PPs",
					value = TRUE) )
			),
		),
		### Main Panel
		mainPanel(
			fluidRow(column(7,
			fluidRow(tag("h1", "Document Term Matrix")),
			fluidRow(DT::DTOutput("tblDTMSummary")),
			),
			column(5, DT::DTOutput("tblDTMFltTerms"))
			),
			# Inspected PP:
			fluidRow(HTML("&nbsp;")),
			fluidRow(textOutput("txtDTM_PP")),
			fluidRow(textOutput("txtDTM_PP_Terms")),
			# Removed Docs:
			fluidRow(h3(textOutput("txtDTM_RemovedDocsTitle"))),
			fluidRow(textOutput("txtDTM_RemovedDocs")),
			fluidRow(DT::DTOutput("tblDTMRemovedDocs")),
			fluidRow(
			column(6,
				fluidRow(textOutput("txtTermsRetained")),
				DT::DTOutput("tblDTMRetainedTerms")),
			column(6,
				fluidRow(textOutput("txtTermsRemoved")),
				DT::DTOutput("tblDTMRemovedTerms")),
			),
		)
	)
}

checkboxNGrams = function(width = 360, inline = TRUE) {
	checkboxGroupInput("chkNGrams", "n-Grams",
		inline = inline, width = width,
		choices = c(
			"2 Ord" = "2", "2 UnOrd" = "2u",
			"3 Ord" = "3", "3 UnOrd" = "3u",
			"4 Ord" = "4", "4 UnOrd" = "4u",
			"PP Length" = "Len",
			"3 Gr-Charge" = "ch3.tot", "3 Gr-Charged AA" = "ch3.aa",
				"3 H-Donor" = "ch3.hd",
			"4 Gr-Charge" = "ch4.tot", "4 Gr-Charged AA" = "ch4.aa",
				"4 H-Donor" = "ch4.hd",
			"5 Gr-Charge" = "ch5.tot", "5 Gr-Charged AA" = "ch5.aa",
				"5 H-Donor" = "ch5.hd"),
		selected = c("2", "2u", "3u", "Len",
			"ch4.tot", "ch4.aa", "ch4.hd")
	)
}
