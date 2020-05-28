##Script for modeling with automodel in a selected region
from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()

env.io.atom_files_directory = ['.', '../atom_files']

# selected atoms do not feel the neighborhood - default is one but this is a first step in which we do not want to change the surroundings so we put it to 2
#env.edat.nonbonded_sel_atoms = 2

#if there are HETATOM to keep into account you can use this command
#env.io.hetatm = True

#set the subclass of automodel
class MyModel(automodel):
   #set actions to renumber residues and label chains after building the models 
   def user_after_single_model(self):
       self.rename_segments(segment_ids=('A'), renumber_residues=[35])
   #define the atoms to be moved during the modelling/optimization, in our case select residues from 506 to 526 from the model numbering
   def select_atoms(self):
       return selection(self.residue_range('328:', '361:'))

#use the subclass defined before
a = MyModel(env,
           alnfile  = 'alignment.ali',      # alignment filename
           knowns   = '2j0j_model1_al',               # codes of the templates
           sequence = 'FAK1_2j0j_yl_1')               # code of the target -> same name of the alignment file
            
a.starting_model= 1                 # index of the first model
a.ending_model  = 1              # index of the last model 


a.make()                            # do modeling and efinement
