# Data
## Folders
### RNA_data
RNA_data contains original data of RNAs, including PDBs and .dat&_ref.pdb which can generate RMSF list through "rna_get_rmsf_list" in rna_util.py, indexed according to PDB ID.
### model
This folder contains the model trained by 5-fold cross-validation.
### example output
This folder contains **3dj2_output.pdb**, **rmsf.png** and **3dj2_input.pdb**.  
**3dj2_output.pdb**: the output PDB with predicted RMSF values which are normalized.  
**rmsf.png**: visualizing output PDB using Chimera.  
**3dj2_input.pdb**: the original PDB.