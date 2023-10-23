# DeepRMSF
## Introduction
DeepRMSF is an automated deep-learning based approach for ‘imaging’ the dynamics of RNA at atomic resolution. Starting with a given PDB-formatted structure, DeepRMSF, it is first translated to density map, structure feature in density map is translated by structure. The maps are segmented into a series of density boxes, which served as input for the model. Finally, predicted RMSF subboxes were then merged into an RMSF map as the RMSF prediction map for this RNA.
## Files
### rna_util.py
This file undergoes the most basic data processing, such as generating simulated maps.
### rna_to_input.py
This file divides the maps into boxes as input to the model.
### model_util.py
Functions required for model training.
### rmsf_model.py
DeepRMSF model.
### rna_five_fold.py
We train and test model using 5-fold cross-validation.
### main.py
## Usage
Download the code file and run **main.py**. python main.py [-h] [--ori_dir ORI_DIR] [--data_dir DATA_DIR] [--box_file BOX_FILE] [--log_dir LOG_DIR]  
You can enter the following parameters,
### --ori_dir
The folder for saving PDBs of RNAs.
### --data_dir
The folder where you want to save the simulated maps.
### --box_file
The folder where you want to save the box_files.
### --log_dir
The folder where you want to save predicted data.
### -h or --help
You can consult the help.
## Input
PDBs of RNAs.
## Output
PDBs with predictive normalized RMSF values which replace B-factor values. The file names are "{pdbid}_pre_nor.pdb". These PDBs can be visualized by Chimera and so on.
## Example
You can view the **test** folder to learn about the DeepRMSF prediction process.
## Supporting softwares
### x3DNA-DSSR
To obtain secondary structure.
### Chimera
To obtain simulated maps.
