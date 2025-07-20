
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
				value = 0.1, min = 0, max = 2, step = 0.025),
			actionButton("btnDTM", "Build DTM"),
			actionButton("btnDTMFilter", "Filter DTM"),
			actionButton("btnDTMInspectPP", "Inspect"),
			# n-Grams:
			fluidRow(HTML("&nbsp;")),
			checkboxGroupInput("chkNGrams", "n-Grams",
				inline = TRUE, width = 360,
				choices = c("2 Ord" = "2", "2 UnOrd" = "2u",
					"3 Ord" = "3", "3 UnOrd" = "3u",
					"4 Ord" = "4", "4 UnOrd" = "4u",
					"PP Length" = "Len",
					"3 Gr-Charge" = "ch3.tot", "3 Gr-Charged AA" = "ch3.aa",
					"4 Gr-Charge" = "ch4.tot", "4 Gr-Charged AA" = "ch4.aa",
					"5 Gr-Charge" = "ch5.tot", "5 Gr-Charged AA" = "ch5.aa"),
				selected = c("2", "2u", "3u", "Len", "ch4.tot", "ch4.aa")
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
		)
	)
}
