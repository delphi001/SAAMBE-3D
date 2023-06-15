#  
#SAMMBE 3D can be used to i) calculate the binding free energy caused by single mutation and ii) predict whether a particular mutation is disruptive or non-disruptive.
#
#User need to install Python3.x and numpy,pandas and xgboost libraries to run the script

#User need to download the script and two models. Regression.model is the model for predicting binding free energy due to mutation and classification.model is the model
#to identify disruptive/non-disruptive mutation
#
#The script can be run in the local computer from the terminal by writing:
#
#python saambe-3d.py -i PDBfile -c Chain -r Resid -w wild -m mutation -d model OR saambe-3d.py -i PDBfile -f mutation_list -d model
#   
#There are two models:0 and 1. 1 is used to predict the binding free energy due to mutation and 0 is used for identifying if the mutation is disruptive or non-disruptive
#
#For example if user want to predict binding free energy for protein-protein complex with PDB ID 1A22 due to mutation from Cysteine (C) to Alanine (A) at res ID 182
#in chain A, they should run the script by typing
#
#python saambe-3d.py -i 1A22.pdb -c A -r 182 -w C -m A -d 1
#
#This way, user can run the prediction for a single mutation. If user want to get multiple predictions for many single mutations at the same complex, user can provide a file
#called 'mutations_list.txt', which should look like the following for the above example
#
#A 182 C A   (###ChainID resID wildtype_residue mutant_residue)
#
#And the command will be: python saambe-3d.py -i 1A22.pdb -f mutations_list.txt -d 1
#
#Similarly. for predicting disruptive/non-disruptive mutation, user can type the same command, just need to change the model from 1 to 0
#
#python saambe-3d.py -i 1A22.pdb -c A -r 182 -w C -m A -d 0 OR python saambe-3d.py -i 1A22.pdb -f mutations_list.txt -d 0
#
#The 'mutations_list.txt' will be exactly same.
