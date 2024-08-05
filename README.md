# Pinto-Reyes-2024-SCRM

## Overview
### Introduction
The Supply Chain Risk Model (SCRM) is a tool that was designed to evaluate how practices applied under produce safety systems could affect microbial contamination, and, our outcome measure of risk: the risk of a test being positive. Produce supply chains can generally be described as having five main process stages from farm to fork: Primary Raw Material Production, Harvest, Processing, Presentation to Consumer (e.g., retail, restaurant, etc.) and Consumer Handling. Model steps can be added to represent each of these supply chain stages, where any of the following model step types could be selected: increase or reduction, contamination or removal, product test, or risk output test. Each step has a given probability of occurrence.

Our initial test case for this model was a leafy green supply chain contaminated with Shiga-toxin-producing *E. coli* (STEC). The parameters of this supply chain are described in the paper and details around the parametrized values and distributions can also be found in the Supplemental Materials folder of this GitHub. We described two contamination scenarios (high and low variability) and two interventions added to the baseline system (improved process controls and additional product testing) to give a robust comparison between contamination variability and the tradeoffs between food safety interventions under these contamination scenarios. The base code for the SCRM described in the manuscript can be found here under the SCRM folder.

To facilitate use of our model the SCRM-Lite was developed as an interactive webpage for users to explore the results of the scenarios described in our leafy green supply chain scenario. Here you will find a description of how to run the SCRM-Lite to visualize the results presented in the manuscript.

## Usage - SCRM Lite
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

#### Step 4.
Click "Run App" (top right corner)

#### Step 5.
You will be launched into the SCRM-Lite interactive webpage. There are three slider bars in total: Contamination Event Variability, Additional Product Testing, and Process Wash. The default selection is high variability baseline contamination scenario described in the manuscript. 

**Contamination Event Variability options:** Low, High

**Additional Product Testing options:** None, Some

**Process Wash options:** Standard, Improved

The following selections in the SCRM-Lite tool match the scenarios in the manuscript: 

SCRM-Lite Options  | Manuscript Scenario
------------- | -------------
High; None; Standard |  High Variability, Baseline
High; Some, Standard  | High Variability, Additional Product Testing
High; None; Improved |  High Variability, Improved Process Controls
High; Some, Improved  | N/A, combination of practices
Low; None; Standard |  Low Variability, Baseline
Low; Some, Standard  | Low Variability, Additional Product Testing
Low; None; Improved |  Low Variability, Improved Process Controls
Low; Some, Improved  | N/A, combination of practices

## Usage - SCRM
### Setup
This app was developed and initialized in R 4.3.3 using Shiny. The packages which are required to run this app include: parallel, lhs, mc2d, bslib, bsicons, shiny, shinyalert, shinyjs, shinyWidgets, DT, beepr

### Running
Below describes how to run the SCRM base code. You will need to access the files in the SCRM folder of this GitHub.

#### Step 1.
Open the R Project file titled "SCRM-Shiny"

#### Step 2.
Open the R file titled "app - V2.1"

#### Step 3.
Install the required packages (if needed). RStudio can detect this.

#### Step 4.
Click "Run App" (top right corner)

#### Step 5.
You will be launched into the SCRM Shiny app. To load and run the scenarios described in the manuscript go to the "Model" tab.

#### Step 6.
Select from the list of Preset Scenarios. Preset Scenarios default to the high variability contamination system. To make steps appear hit the "Load Scenario" button. The following table describes the preset scenario as it relates to description and any modifications which should need to be made to reproduce the results: 
Preset Scenario Name  | Modifications Needed | Manuscript Scenario | 
------------- | ------------- | -------------
Baseline |  N/A | High Variability, Baseline
FPT | N/A |  High Variability, Additional Product Testing
Improved Process Controls |  N/A | High Variability, Improved Process Controls
Baseline |  Change "Initial Contamination" step (Step 01) SD to value of 0.2 | Low Variability, Baseline
FPT | Change "Initial Contamination" step (Step 01) SD to value of 0.2 | Low Variability, Additional Product Testing 
Improved Process Controls | Change "Initial Contamination" step (Step 01) SD to value of 0.2 | Low Variability, Improved Process Controls 

#### Step 7.
The default field mass will = 160,000; the default mass unit will = lb; the default iterations will = 100,000. Leave these defaults to reproduce the results of the manuscript.

#### Step 8. (optional)
You can check the "beep" box if you would like your computer to make a sound once the model has completed running all iterations

#### Step 9.
Click "Run". Depending on the scenario, you will be notified that there are step(s) present after at risk output test at retail (current model output). This was coded for future use of the model to force the user to check for errors. We parametrized a Consumer Handling stage in our manuscript, but this does not currently affect our decided location of the risk output test. However, these values are still pulled into our baseline. One could simply uncheck the Consumer Handling model step in the SCRM Shiny app to avoid this pop up error, or simply click "Continue." 

A small pop-up will appear at the bottom right of the app once the model has begun running. This pop up will let you know the model is running and provide a time stamp for the start.

### Step 10.
After successfully running the model, a second small pop-up will appear at the bottom right of the app to notify you the model has completed running. Your results will appear in all 3 tables.

### Step 11.
You can save results from a single run by hitting "Download Results" under each table, or by loading them into slots (equivalent to table rows) so that they appear on the Results Comparisons tab. If you are running multiple scenarios, we suggest saving the results to slots after running each scenario. You can save these by selecting your desired slot and clicking the "Save to Slot" button. This will populate your results into two tables on the Results Comparisons tab.


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
