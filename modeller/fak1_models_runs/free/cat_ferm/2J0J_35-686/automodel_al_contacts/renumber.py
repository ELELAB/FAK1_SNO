from Bio.PDB.PDBParser import *
from Bio.PDB import PDBIO

name=raw_input("PDB name:")
shift=raw_input("New start number:")
parser = PDBParser()

structure = parser.get_structure(name.replace(".pdb",""), name)
header = parser.get_header()
trailer = parser.get_trailer()

model = structure[0]
chain = model['B']
i=int(shift)
for residue in chain:
    residue.id = (' ', i, ' ')
    i += 1

w = PDBIO()
w.set_structure(structure)
w.save('renumbered_'+name)

