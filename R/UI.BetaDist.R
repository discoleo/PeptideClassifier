#####################
##
## Peptide Classifier
##
## Leonard Mada
##
## URL: https://github.com/discoleo/PeptideClassifier
##
## draft v.0.1a

### Generalised Beta Distribution

panelBetaDist = function() {
	sidebarLayout(
		sidebarPanel(
			# R = Radical;
			fluidRow(tag("h3", "Generalised Beta Distribution:")),
			fluidRow("Double Int(  x^p * y^q * (1 - x^m*y^n)^k ) over [0,1]^2"),
			fluidRow(
			column(6, textInput("inBeta_xPow", "x^p", "sqrt(2)", width = 150)),
			column(6, textInput("inBeta_yPow", "y^q", "sqrt(3)", width = 150)),
			),
			fluidRow(
			column(6, textInput("inBeta_RxPow", "x^m", "sqrt(2) - 2/3", width = 150)),
			column(6, textInput("inBeta_RyPow", "y^n", "sqrt(5)", width = 150)),
			),
			fluidRow(
			column(6, textInput("inBeta_RPow", "()^k", "1 / sqrt(13)", width = 150)),
			),
			actionButton("btnBetaCompute", "Compute"),
		),
		### Main Panel
		mainPanel(
			# Topic Models:
			fluidRow(tag("h1", "Result:")),
			fluidRow(textOutput("txtBetaResult")),
			fluidRow(textOutput("txtBetaNumeric")),
		)
	)
}