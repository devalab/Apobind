from prody import *
import argparse
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import AllChem
from pymol import cmd
import numpy as np
from scipy.spatial.distance import cdist
import os
from multiprocessing import Pool, cpu_count

def get_apo(id):
    dir='/scratch/rishal/v2019-other-PL/'+id+'/'
    protein_file=dir+id+'_protein.pdb'
    ligand_file=dir+id+'_ligand.sdf'
    pathPDBFolder('/scratch/rishal/apo_Structs/')
    apo_file=open('apo_structs.txt','a')
    blast_fail=open('blast_fail.txt','a')
    #load protein file, get final chain seq, save to fasta and get entire protein length

    protein=parsePDB(protein_file)
    chains=sorted(list(set(protein.getChids())))
    seq=protein[chains[-1]].getSequence()
    fasta=open(dir+'fasta.txt','w')
    fasta.write('>'+id+'\n')
    fasta.write(seq)
    fasta.close()
    protein_seq=''
    for chain in chains:
        if chain==' ':
            continue
        else:
            protein_seq+=protein[chain].getSequence()
    protein_len=len(protein_seq)

    #load ligand file and get ligand coords
    mol = Chem.MolFromMolFile(ligand_file,sanitize=False)
    c=mol.GetConformer()
    ligand_coords=c.GetPositions()
    
    #blast sequence and output results
    try:
        os.system('/scratch/rishal/ncbi-blast-2.11.0+/bin/blastp -query '+dir+'fasta.txt -db /scratch/rishal/pdbaa -out '+dir+'blast.txt -outfmt 6')
        blast_results=open(dir+'blast.txt','r').readlines()
    except:
        blast_fail.write(id+'\n')
        return    
    #iterate through hits, full alignment with pymol and check if ligands are close
    checked_pdbs=[]
    for line in blast_results:
        print(line)
        pdb_id=line.split()[1].split('_')[0].lower()
        align_identity=float(line.split()[2])
        if pdb_id in checked_pdbs:
            continue
        if align_identity<80:
            continue
        if pdb_id==id:
            continue
        checked_pdbs.append(pdb_id)
        pdb=fetchPDB(pdb_id,tp='http')
        if pdb is None:
            continue
        cmd.load(pdb)
        cmd.load(protein_file)
        cmd.align( pdb_id+'////CA', id+'////CA',object='aln')
        cmd.align( pdb_id+'////CA', id+'////CA',cycles=0, transform=0,object='aln')
        raw_aln = cmd.get_raw_alignment('aln')
        idx2resn = {}
        idx2chain = {}
        cmd.iterate('aln', 'idx2resn[model, index] = resn', space={'idx2resn': idx2resn})
        cmd.iterate('aln', 'idx2chain[model, index] = chain', space={'idx2chain': idx2chain})
        count_identity=0
        chains_hit=[]
        for idx1, idx2 in raw_aln:
            if idx2resn[idx1]==idx2resn[idx2]:
                count_identity+=1
            chains_hit.append(idx2chain[idx2])
        chains_hit=sorted(list(set(chains_hit)))
        print(chains_hit)
        chain_seq=''
        for chain in chains_hit:
            chain_seq+= 'chain '+chain+' or '
        chain_seq=' '.join(chain_seq.split()[:-1])
        cmd.save(dir+'apo_temp.pdb',pdb_id)
        cmd.delete(pdb_id)
        cmd.delete(id)
        apo_prot=parsePDB(dir+'apo_temp.pdb')
        apo_seq=''
        for chain in chains_hit:
            apo_seq+=apo_prot[chain].getSequence()
        hit_len=len(apo_seq)
        # sequence identity and coverage thresholded at 0.8
        if count_identity/hit_len < 0.8 or count_identity/protein_len < 0.8:
            os.remove(dir+'apo_temp.pdb')
            continue
        ligand_apo = apo_prot.select('not protein and not water and not ion')
        if ligand_apo is None:
            apo_file.write(id+'\n')
            apo_prot=apo_prot.select(chain_seq + ' or ion')
            apo_prot=apo_prot.select('protein or ion and not water')
            #writePDB(dir+id+'_apo.pdb',apo_prot)
            #os.remove(dir+'apo_temp.pdb')
            break
        #check if hit has ligand in same binding pocket
        ligand_apo_coords=ligand_apo.getCoords()
        dist=cdist(ligand_coords,ligand_apo_coords)
        if np.min(dist) <=4:
            os.remove(dir+'apo_temp.pdb')
            continue
        else:
            apo_file.write(id+'\n')
            apo_prot=apo_prot.select(chain_seq + ' or ion')
            apo_prot=apo_prot.select('protein or ion and not water')
            #writePDB(dir+id+'_apo.pdb',apo_prot)
            #os.remove(dir+'apo_temp.pdb')
            break


'''parser = argparse.ArgumentParser()
parser.add_argument('-s', '--string', type=str, help="comma seperated pdb ids for apo query", required=True)
opt = parser.parse_args()
ids=opt.string.split(',')
ncpus=cpu_count()
pool = Pool(ncpus)
pool.map(get_apo,ids)
'''
#for id in ids:
#    get_apo(id)

get_apo('2aq7')
