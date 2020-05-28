####repository of PDB files that we use in the project #### 
#to download PDB files  with PDB entry name
./fetch.pdb.sh 2J0J
./fetch.pdb.sh 2J0K

#check with pymol if we have one or more chains
#if some of them have more chains we save also the PDB of each of them, such as 2J0K_A.pdb and 2J0K_B.pdb

#to get a WHATIF report use 
#https://swift.cmbi.umcn.nl/servers/html/index.html Build/Check/Repair model -> Template Structure Check. and save for each pdb the pdbout.txt, as for example, 2J0J_pdbout.txt - we need a copy in these folders.
#N.B. the analysis takes few minutes
#you will download you file locally on your laptop then you need to scp in the whatif subfolder here, example
#to do from your local computer to the server
scp pdbout_2J0J.txt elena@kbf-bioinfo01.cancer.dk:/data/user/shared_projects/fak1_sno_paper_2020/pdbs/whatif  

#for 2J0K better to do the whaif for each of the two chains, i.e. pdbout_2J0K_A.txt and pdbout_2J0K_B.txt

