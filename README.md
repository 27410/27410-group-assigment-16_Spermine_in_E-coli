[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-c66648af7eb3fe8bc4f294546bfd86ef473780cde1dea487d3c4ff354943c9ae.svg)](https://classroom.github.com/online_ide?assignment_repo_id=9067129&assignment_repo_type=AssignmentRepo)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/27410/27410-group-assigment-16_Spermine_in_E-coli/main)

# 27410 - Group assignment - Group 16 - Spermine production in growing *Escherichia coli*

## Project summary
With the emergence of prion disease, it has become crucial for modern biomedicine to arm itself against this elusive menace. Multiple polyamines are to be assessed as prion disease treatment notably by playing the role of protein chaperone, we choose Spermine as our target product. The purpose of this project was to lay the foundation for a potential cell factory design to produce Spermine in vivo. We computationally went through all relevant models hosted on BiGG and determined that iAF1260b was the best fit for our combined objective of growth and production. Plotting some phenotypic phase planes, we determined that in order to produce the maximum of Spermine the model needed a medium rich in D-Glucose and O2. Once we established the medium composition, with the algorithm described in the report we established a series of reaction modifications to increase the flux of Spermine up to 99.9% of the theoretical Spermine yield for this model. Some up-regulations alone could increase it up to 98%, and most of them were not expected to have such an effect.
Our achievements during this project:
* The computational assessment of the maximum value for the quadratic objective of growth and Spermine production for 102 BiGG models
* We determined the conditions needed for the cell factory to produce and grow in the best conditions
* We ran multiple FVA of the full iAF1260b model
* We established a list of modifications and ranked them by their capacity to increase Spermine flux 

## Project overview
* **Txt files**: “spmd_models.txt” is an input file where we initially copy pasted the list of all BiGG models containing Spermidine, “Optimization startegy.txt” is the output file listing the reactions regulations found by our algorithm
* **CSV files**: Both FVA_BIOMASS.csv and FVA_Q_OBJECTIVE.csv are data files that contain respectively the result of the FVA for the Wild type model and the Theoretical optimal Strain 
* **GSM folder**: a folder containing all spermidine GSM, only iAF1260b.xml is needed for later analysis

### Python scripts: 
* **GSM_assessment.py** performs the analysis of the BiGG models maximum Spermine flux 
* **Heterologuous_pathway_addition.py** allows the addition of all needed metabolites and reactions to a model (main() function) and also verify * that requirements for Spermine synthesis are present in the model (model_assertion() function)
* **Phenotipic_phase_planing.py** performs phenotypic phase plane analysis and establishes the optimal medium composition
* **FVA_script.py** performs FVA for WT and OS and stores results in a CSV file
* **Optimization.py** contains the algorithm that returns a ranked reaction modifications list and stores them in “Optimization startegy.txt”

