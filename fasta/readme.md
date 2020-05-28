#collect the fasta file for the project
#UNIPROT fasta downloaded from the corresponding entries in UNIPROT
#Q00944_FAK1_CHICK and Q05397_FAK1_HUMAN (fak1_chicken.txt  fak1_human.txt) 
#fasta files for the PDB structure to be extracted with pdb2fasta.pl
#create a symbolic link to the PDB files you need to use
ln -s ../pdbs/2J0J.pdb
ln -s ../pdbs/2J0K_A.pdb
ln -s ../pdbs/2J0K_B.pdb

#2J0J.pdb has only chain A so we can use it as it is
#2J0K.pdb has two identical chains (homodimer)

#run pdb2fasta.pl for each of them in this way:
pdb2fasta.pl (press ENTER)
write the PDB file name
#it will give the fasta in output on the screen - we can copy it in a txt file
#to call with the PDB file name, for example, 2J0J_fasta.txt



