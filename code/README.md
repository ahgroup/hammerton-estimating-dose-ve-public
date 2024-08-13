This folder contains all the code used to conduct our analysis and generate results. The code produces output, tables, and figures, which are saved in the `results` folder.

  * The `processing-code.R` script uses the raw HAI titer data to estimate protection against influenza, primarily used for EDA
  * The `statistical-analysis.R` script performs most of the primary analysis, estimating VEs and DVEs with bootstrapped 95% confidence intervals
  * The `fig-gen.R` and `table-gen.R` scripts contain all the code used to generate tables and figures presented in the manuscript
  * The `objects.R` script contains basic repeatedly used objects such as manual color and shape scales for plots
  * The `functions.R` script contains repeatedly used functions that calculate VE and DVE with their 95% CIs

Both the `objects.R` and `functions.R` scripts are called at the beginning of each other script, so they do not need to be opened and run, but the `functions.R` script particularly contains much of the code actually used to calculate results, so any alterations may need to be implemented within that script. 
