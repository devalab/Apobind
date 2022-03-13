import rdkit
import rdkit.Chem
from rdkit import Chem
import rdkit.Chem.Descriptors
from rdkit.Chem.Descriptors import ExactMolWt,HeavyAtomMolWt
from clean_pdb import clean_pdb
import os
count=0
count_prot=0
f=open('exclude_structs.txt','w')
for dir in sorted(os.listdir('/scratch/rishal/v2019-other-PL')):
    if len(dir)!=4:
        continue
    file=dir+'_protein.pdb'
    mol = Chem.MolFromMolFile(
        os.path.join('/scratch/rishal/v2019-other-PL', dir, file.replace('protein.pdb', 'ligand.sdf')), sanitize=False)
    if mol is None:
        f.write(dir + '\n')
        count += 1
        continue
    mol_wt = HeavyAtomMolWt(mol)
    if mol_wt > 1000:
        f.write(dir + '\n')
        count += 1
        continue
    try:
        clean_pdb(os.path.join('/scratch/rishal/v2019-other-PL', dir, file),
                  os.path.join('/scratch/rishal/v2019-other-PL', dir, file.replace('protein', 'protein_nowat')))
    except:
        f.write(dir+'\n')
        count+=1
        continue


print(count)