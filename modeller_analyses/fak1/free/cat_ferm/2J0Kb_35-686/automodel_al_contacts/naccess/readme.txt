#NB the selected models after visual inspection with pymol to remove collisions
Model2 
Model55
Model84

# we copy them with a new name so that naccess can recognize them 

cp ../../../../../../../modeller/fak1/free/cat_ferm/2J0Kb_35-686/automodel_al_yl_native/model2/FAK1_2j0kB_yl_2.B99990001.pdb model2.pdb 
cp ../../../../../../../modeller/fak1/free/cat_ferm/2J0Kb_35-686/automodel_al_yl_native/model55/FAK1_2j0kB_yl_55.B99990001.pdb model55.pdb
cp ../../../../../../../modeller/fak1/free/cat_ferm/2J0Kb_35-686/automodel_al_yl_native/model84/FAK1_2j0kB_yl_84.B99990001.pdb model84.pdb

./run_naccess.sh #model names should be changed in the sh script

#we need only the .rsa file / .log and .asa can be removed

#we would need to add a R-script to extract the Side Chain (REL) accessibility
#on the residues in the region 509-524 and their average as a reference for each model 
