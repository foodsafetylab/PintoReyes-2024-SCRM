# Pinto-Reyes-2024-SCRM

## Overview
### Introduction
The Supply Chain Risk Model (SCRM) is a tool that was designed to evaluate how practices applied under produce safety systems could affect microbial contamination, and, our outcome measure of risk: the risk of a test being positive. Produce supply chains can generally be described as having five main process stages from farm to fork: Primary Raw Material Production, Harvest, Processing, Presentation to Consumer (e.g., retail, restaurant, etc.) and Consumer Handling. Model steps can be added to represent each of these supply chain stages, where any of the following model step types could be selected: increase or reduction, contamination or removal, product test, or risk output test. Each step has a given probability of occurrence.

Our initial test case for this model was a leafy green supply chain contaminated with Shiga-toxin-producing E. coli (STEC). The parameters of this supply chain are described in the paper and details around the parametrized values and distributions can also be found in the Supplemental Materials. We described two contamination scenarios (high and low variability) and two interventions added to the baseline system (improved process controls and additional product testing) to give a robust comparison between contamination variability and the tradeoffs between food safety interventions under these contamination scenarios. 

To facilitate use of our model the SCRM-Lite was developed as an interactive webpage for users to explore the results of the scenarios described in our leafy green supply chain scenario. Here you will find a description of how to run the SCRM-Lite to visualize the results presented in the manuscript.

## Usage
### Setup
This app was developed and initialized in R 4.3.3 using Shiny. The packages which are required to run this app include: shiny, shinyjs, bslib, shinyWidgets, and ggplot2

### Running
Below describes how to run the SCRM-Lite. You will need to access the files in the SCRM-Lite folder of this GitHub.

#### Step 1.
Open the R Project file titled "SCRM-Lite"

#### Step 2.
Open the R file titled "app"

#### Step 3.
Install the required packages (if needed). RStudio can detect this.

### Step 4.
Click "Run App" (top right corner)

### Step 5.
You will be launched into the SCRM-Lite interactive webpage. There are three slider bars in total: Contamination Event Variability, Additional Product Testing, and Process Wash. The default selection is high variability baseline contamination scenario described in the manuscript. The following selections in the SCRM-Lite tool match the scenarios in the manuscript: 

SCRM-Lite Options  | Manuscript Scenario
------------- | -------------
High Variability; None; Standard |  High Variability Baseline
High Variability; Some, Standard  | High Variability Additional Product Testing
High Variability; None; Improved |  High Variability Improved Process Controls
High Variability; Some, Improved  | N/A, combination of practices
Low Variability; None; Standard |  Low Variability Baseline
Low Variability; Some, Standard  | Low Variability Additional Product Testing
Low Variability; None; Improved |  Low Variability Improved Process Controls
Low Variability; Some, Improved  | N/A, combination of practices


## Authors
You can view the list of authors in the [AUTHORS](/AUTHORS) file.

## Contact
Corresponding author: Matthew J. Stasiewicz<br>
103 Agricultural Bioprocess Lab<br>
1302 W. Pennsylvania<br>
Urbana, IL, 1361801<br>
USA<br>
+1-217-265-0963<br>
[mstasie@illinois.edu](mailto:mstasie@illinois.edu)

## Citation
TBD

## Funding
Funding for the project was made possible by The Center for Produce Safety project 2023CPS08. Any opinions, findings, conclusions, or recommendations expressed in this publication are those of the authors and do not necessarily reflect the view of The Center for Produce Safety. (https://www.centerforproducesafety.org/).
