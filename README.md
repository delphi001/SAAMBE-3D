#  SAAMBE-3D: Predicting Effect of Mutations on Protein–Protein Interactions.
# About
SAAMBE-3D is a tool written in python 3, it is tested to work with python3.7. SAAMBE-3D implements a structure based mathod for predicting the binding free energy change for protein-protein due to a point-mutation. It requires the 3D structure in PDB of the protein-protein complex as input and information about the mutation to make the prediction.
##### SAMMBE 3D can be used to :
- i) predict the binding free energy change caused by single mutation, and 
- ii) predict whether a particular mutation is disruptive or non-disruptive.

# Dependepcies
SAAMBE-3D requires following packages and mentioned versions to be installed in the python (Python 3.11.5) to work
- numpy (1.26.0)
- prody (2.4.0)
- xgboost (1.7.6)
  
It has been tested to work on a machine running Ubunutu 22.04.3. The python is mainted using conda 23.7.2. We created a python environment for SAAMBE-3D using:
```sh
conda create -n py311_saambe3d python=3.11
```
Then if any conda environment is already active deactivate it using:
```sh
conda deactivate
```
Activate the conda python environemnt created above using:
```sh
conda activate py311_saambe3d
```
All the above versions of packages are installed using following commands.
```sh
conda install -c conda-forge numpy=1.26.0
conda install -c conda-forge prody=2.4.0
conda install -c conda-forge numpy=1.7.6
```

# Prediction Models
The method provides two different models for predictions:
- regression.model (for predicting: binding free energy change)
- classification.model (for prediction: if the mutation is disruptive or non-disruptive)

# How to run
To have a look at and familiarize yourself with all the options supported use:
```sh
python saambe-3d.py --help
```
The SAAMBE-3D can be run for a single point mutation by executing following command from inside of the SAAMBE-3D directory:
```sh
python saambe-3d.py -i PDBfile -c Chain -r Resid -w wild -m mutation -d model 
```
To make predictionsof a list of point-mutations execute:

```sh
saambe-3d.py -i PDBfile -f mutation_list -d model
```
# Example
##### Single point-mutation: 
For example if user want to predict binding free energy for protein-protein complex with PDB ID 1A22 due to mutation from Cysteine (C) to Alanine (A) at res ID 182 in chain A, they should run the script by typing

```sh
python saambe-3d.py -i 1A22.pdb -c A -r 182 -w C -m A -d 1
```
After running it check the file *`output.out`* or the file name if you provided any with ```sammbe-3d.py -o option``` for the output.

#### List of point-mutations. 
If user want to get multiple predictions for many single mutations at the same complex, user can provide a file
say `mutations_list.txt`, which should be formatted as follows for the above example.
The columns in the mutation list file are: `ChainID` `resID` `wildtype_residue` `mutant_residue`
```sh
A 182 C A   
```
and the command will be: 
```sh
python saambe-3d.py -i 1A22.pdb -f mutations_list.txt -d 1
```
After running it check the file *`output.out`* or the file name if you provided any with ```sammbe-3d.py -o option``` for the output.

#
Similarly. for predicting disruptive/non-disruptive mutation, user can type the same command, just need to change the model from 1 to 0
```sh
python saambe-3d.py -i 1A22.pdb -c A -r 182 -w C -m A -d 0 
```
After running it check the file *`output.out`* or the file name if you provided any with ```sammbe-3d.py -o option``` for the output.

OR
```sh
python saambe-3d.py -i 1A22.pdb -f mutations_list.txt -d 0
```
After running it check the file *`output.out`* or the file name if you provided any with ```sammbe-3d.py -o option``` for the output.
The 'mutations_list.txt' will be exactly same.

# Citation
Pahari, S., Li, G., Murthy, A. K., Liang, S., Fragoza, R., Yu, H., & Alexov, E. (2020). SAAMBE-3D: predicting effect of mutations on protein–protein interactions. *International journal of molecular sciences, 21(7)*, 2563.
