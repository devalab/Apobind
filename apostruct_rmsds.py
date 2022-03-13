import argparse
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import AllChem
from prody import *
from Bio import pairwise2
from scipy.spatial.distance import cdist
import numpy as np
import os
import operator
import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count,Queue
from functools import partial
import pandas as pd

def get_chain_from_index(index,protein_chain2seqind):
    diff=np.inf
    chain=None
    for key in protein_chain2seqind.keys():
        if index >= protein_chain2seqind[key]:
            ind_diff=index-protein_chain2seqind[key]
            if ind_diff<diff:
                diff=ind_diff
                chain=key
    return chain
def get_indices(alignment):
    indices=[]
    count=-1
    for i in alignment:
        if i =='-':
            indices.append(i)
        else:
            count+=1
            indices.append(count)
    return indices

def get_chains_sequence(protein,protein_chains):
    protein_seq=[]
    protein_chainseq=[]
    protein_chain2seqind={}
    #protein_chains=sorted(list(set(protein.getChids())))
    for chain in protein_chains:
        if chain==' ': #ion
            continue
        else:
            try: #sometimes pymol gives chains not really involved in alignment
                protein_chainseq += [chain]*len(protein[chain].getSequence())
            except:
                continue
            protein_chain2seqind[chain] = len(protein_seq)
            protein_seq+=protein[chain].getSequence()
    assert len(protein_chainseq)==len(protein_seq)
    protein_pdbchains=protein.getChids()
    protein_residues=protein.getResindices()
    protein_chain2resind={}
    for chain,resind in zip(protein_pdbchains,protein_residues):
        if chain==' ': #ion
            continue
        first_chain=chain
        protein_chain2resind[first_chain] = resind
        break
    chain_current=first_chain
    for chain,resind in zip(protein_pdbchains,protein_residues):
        if chain==' ': #ion
            continue
        if chain != chain_current:
            chain_current=chain
            protein_chain2resind[chain_current] = resind
    return protein_seq,protein_chainseq,protein_chain2seqind,protein_chain2resind,protein.getCoords(),protein_residues,protein_pdbchains

def get_apo_bind_indices(holo_binding_resind,holo_binding_chains,holo_align_index,holo_seq,apo_align_index,apo_seq,holo_chain2resind,holo_chain2seqind,apo_chain2resind,apo_chain2seqind):
    apo_bind_indices=[]
    holo_bind_indices=[]
    apo_bind_indices_all = []
    holo_bind_indices_all = []
    holo_bind_res=[]
    apo_bind_res=[]
    for holo_chain,holo_res_ind in zip(holo_binding_chains,holo_binding_resind):
        holo_seq_ind=holo_res_ind - holo_chain2resind[holo_chain] + holo_chain2seqind[holo_chain]
        holo_aln_index=holo_align_index.index(holo_seq_ind)
        holo_res = holo_seq[holo_seq_ind]
        holo_bind_res.append(holo_res)
        holo_bind_indices_all.append(holo_res_ind)
        apo_seq_index=apo_align_index[holo_aln_index]
        if apo_seq_index=='-':
            apo_bind_res.append('-')
            apo_bind_indices_all.append('-')
            continue
        apo_res=apo_seq[apo_seq_index]
        apo_chain = get_chain_from_index(apo_seq_index, apo_chain2seqind)
        apo_res_ind = apo_seq_index - apo_chain2seqind[apo_chain] + apo_chain2resind[apo_chain]
        apo_bind_res.append(apo_res)
        apo_bind_indices_all.append(apo_res_ind)
        if apo_res==holo_res:
            apo_bind_indices.append(apo_res_ind)
            holo_bind_indices.append(holo_res_ind)
    return holo_bind_indices,apo_bind_indices,holo_bind_indices_all,apo_bind_indices_all,holo_bind_res,apo_bind_res
    
def get_selection(protein,protein_binding_indices):
    sel_str= 'resindex '
    for i in protein_binding_indices:
        sel_str+= str(i)+' or resindex '
    sel_str=' '.join(sel_str.split()[:-2])
    return protein.select(sel_str)

def get_tmscore(holo_alignment,apo_alignment,holo_align_index,apo_align_index,holo_chainseq,apo_chainseq,holo_chain2seqind,apo_chain2seqind,holo_chain2resind,apo_chain2resind,holo_prody,apo_prody):
    holo_calphas=[]
    apo_calphas=[]
    #formula from wikipedia
    Ltarget=len(holo_chainseq)
    d0=1.24 * (Ltarget - 15) ** (1.0 / 3.0) - 1.8
    d02=d0 ** 2
    assert  len(holo_alignment)==len(apo_alignment)
    for i in range(len(holo_alignment)):
        holo_res=holo_alignment[i]
        apo_res=apo_alignment[i]
        if holo_res=='-':
            continue
        if holo_res==apo_res:
            # go from align index to seq index to res index
            holo_seq_ind=holo_align_index[i]
            apo_seq_ind=apo_align_index[i]
            holo_res_ind=holo_seq_ind-holo_chain2seqind[holo_chainseq[holo_seq_ind]]+holo_chain2resind[holo_chainseq[holo_seq_ind]]
            apo_res_ind=apo_seq_ind-apo_chain2seqind[apo_chainseq[apo_seq_ind]]+apo_chain2resind[apo_chainseq[apo_seq_ind]]
            holo_ca=holo_prody.select('resindex '+str(holo_res_ind)+' and name CA')
            apo_ca = apo_prody.select('resindex ' + str(apo_res_ind) + ' and name CA')
            try: #some residues have CA missing surprisingly - wont affect Tm-score too much by excluding them
                apo_calphas_coords=apo_ca.getCoords()
                #sometimes more than 1 calpha present, take care by taking only the first occurence
                if apo_calphas_coords.shape[0]>1:
                    apo_calphas_coords=np.expand_dims(apo_calphas_coords[0],axis=0)
                apo_calphas.append(apo_calphas_coords)
            except:
                continue
            try:
                holo_calphas_coords = holo_ca.getCoords()
                # sometimes more than 1 calpha present, take care by taking only the first occurence
                if holo_calphas_coords.shape[0] > 1:
                    holo_calphas_coords = np.expand_dims(holo_calphas_coords[0], axis=0)
                holo_calphas.append(holo_calphas_coords)
            except:
                apo_calphas.pop(-1)
                continue
    assert len(holo_calphas)==len(apo_calphas)
    holo_calphas=np.array(holo_calphas)
    apo_calphas=np.array(apo_calphas)
    #apply tmscore formula
    d1s=np.sum(np.sum(np.power(holo_calphas - apo_calphas, 2), axis=1),axis=1)
    tmarray= 1/(1+np.divide(d1s,d02))
    tmscore=(1/Ltarget) * np.sum(tmarray)
    return tmscore
def get_rmsds(id,holo_chains,apo_chains):
    dir='/scratch/rishal/v2019-other-PL/'+id+'/'
    protein_file=dir+id+'_protein_nowat.pdb'
    ligand_file=dir+id+'_ligand.sdf'
    apo_file=dir+id+'_apo_added.pdb'
    #load files
    holo_prody=parsePDB(protein_file)
    holo_prody=holo_prody.select('protein')
    apo_prody=parsePDB(apo_file)
    apo_prody=apo_prody.select('protein')
    # get holo and apo sequences, chain to sequence index mapping, coordinates and residue indices from pdb file
    holo_seq,holo_chainseq,holo_chain2seqind,holo_chain2resind,holo_coords,holo_resin,holo_chains=get_chains_sequence(holo_prody,holo_chains)
    apo_seq,apo_chainseq,apo_chain2seqind,apo_chain2resind,apo_coords,apo_resin,apo_chains=get_chains_sequence(apo_prody,apo_chains)
    #alignments and alignment to sequence mapping
    alignment=pairwise2.align.globalxx(''.join(holo_seq), ''.join(apo_seq))[0]
    holo_alignment=alignment[0]
    holo_align_index=get_indices(list(holo_alignment))
    apo_alignment=alignment[1]
    apo_align_index=get_indices(list(apo_alignment))
    #get tmscore
    tmscore=get_tmscore(holo_alignment,apo_alignment,holo_align_index,apo_align_index,holo_chainseq,apo_chainseq,holo_chain2seqind,apo_chain2seqind,holo_chain2resind,apo_chain2resind,holo_prody,apo_prody)
    #load ligand and get holo binding residue indices from PDB
    mol = Chem.MolFromMolFile(ligand_file,sanitize=False)
    c=mol.GetConformer()
    ligand_coords=c.GetPositions()
    ligand_dist=cdist(ligand_coords,holo_coords)
    binding_indices=np.where(np.any(ligand_dist<=6, axis=0))
    # get binding residues and corresponding chains
    holo_binding_resind=holo_resin[binding_indices]
    holo_binding_chains=holo_chains[binding_indices]
    new_holo_binding_chains=[]
    new_holo_binding_resind=[]
    for chain,resind in zip(holo_binding_chains,holo_binding_resind):
        if chain==' ': #ion
            continue
        if resind in new_holo_binding_resind:
            continue
        new_holo_binding_resind.append(resind)
        new_holo_binding_chains.append(chain)
    holo_binding_resind=new_holo_binding_resind
    holo_binding_chains=new_holo_binding_chains
    assert len(holo_binding_resind) == len(holo_binding_chains)
    #map to alignment and sequence and back to apo_binding residues
    holo_binding_indices,apo_binding_indices,holo_bind_indices_all,apo_bind_indices_all,holo_bind_res,apo_bind_res=get_apo_bind_indices(holo_binding_resind,holo_binding_chains,holo_align_index,holo_seq,apo_align_index,apo_seq,holo_chain2resind,holo_chain2seqind,apo_chain2resind,apo_chain2seqind)
    sc_dist=np.array([])
    for holo_id,apo_id in zip(holo_binding_indices,apo_binding_indices):
        holo_bind=holo_prody.select('resindex '+str(holo_id)+' and not element H and not name CA and not name C and not name O and not name N ')
        apo_bind=apo_prody.select('resindex '+str(apo_id)+' and not element H and not name CA and not name C and not name O and not name N ')
        #side chains should be present
        if holo_bind is None or apo_bind is None:
            continue
        holo_bind_names=holo_bind.getNames()
        enumerate_holo = enumerate(holo_bind_names)
        holo_sorted_pairs = sorted(enumerate_holo, key=operator.itemgetter(1))
        apo_bind_names=apo_bind.getNames()
        #check for missing atoms/non standard residues
        if len(holo_bind_names)!=len(apo_bind_names):
            continue
        enumerate_apo = enumerate(apo_bind_names)
        apo_sorted_pairs = sorted(enumerate_apo, key=operator.itemgetter(1))
        holo_sorted_indices = []
        for index, element in holo_sorted_pairs:
            holo_sorted_indices.append(index)
        apo_sorted_indices = []
        for index, element in apo_sorted_pairs:
            apo_sorted_indices.append(index)
        sc_dist=np.concatenate((sc_dist,(np.sum(np.power(holo_bind.getCoords()[holo_sorted_indices]-apo_bind.getCoords()[apo_sorted_indices],2),axis=1)))) 
    sc_rmsd=np.power(np.mean(sc_dist),0.5)
    return sc_rmsd,tmscore,holo_bind_indices_all,apo_bind_indices_all,holo_bind_res,apo_bind_res
'''if sc_rmsd>15:
        print(id)
    sc_rmsds.append(sc_rmsd)'''
    
bb_rmsds=[]
sc_rmsds=[]
reject=[]
direcs=[]
'''for direc in os.listdir("/scratch/rishal/v2019-other-PL/"):
    if os.path.exists('/scratch/rishal/v2019-other-PL/'+direc+'/'+direc+'_apo.pdb'):
        if not os.path.exists('/scratch/rishal/v2019-other-PL/'+direc+'/'+direc+'_protein_nowat.pdb'):
            continue
        get_rmsds(bb_rmsds,sc_rmsds,reject,direc)'''
'''function1=partial(get_rmsds,bb_rmsds,sc_rmsds,reject)
ncpus=cpu_count()
pool = Pool(ncpus)
pool.map(function1,direcs)'''
'''print(reject)
print(len(reject))
f=open('reject_apo.txt','w')
f.write('\n'.join(reject))
np.save('bb_rmsds_icml.npy',np.array(bb_rmsds))
np.save('sc_rmsds_icml.npy',np.array(sc_rmsds))
ax=sns.distplot(bb_rmsds,kde=False)
plt.rcParams.update({'font.size': 12})
ax.set(ylabel="Protein count", xlabel="C$_\alpha$ RMSD")
plt.savefig('bb_rmsds.png')
plt.close()
plt.figure()
ax=sns.distplot(sc_rmsds,kde=False,color='orange')
plt.rcParams.update({'font.size': 12})
ax.set(ylabel="Protein count", xlabel="Side Chain RMSD")
plt.savefig('sc_rmsds.png')
'''
apo_file=open('apo_structs_sc_tm.csv','w')
apo_file.write('holo_id,tmscore,holo_bind_res,holo_bind_indices,apo_bind_res,apo_bind_indices,side_chain_rmsd\n')
for direc in os.listdir("/scratch/rishal/v2019-other-PL/"):
    if os.path.exists('/scratch/rishal/v2019-other-PL/'+direc+'/'+direc+'_apo.pdb'):
        if not os.path.exists('/scratch/rishal/v2019-other-PL/'+direc+'/'+direc+'_protein_nowat.pdb'):
            continue
        print(direc)
        apo_csv=pd.read_csv('apo_structs_ids.csv')
        id_csv=apo_csv[apo_csv['holo_id']==direc]
        holo_chains=id_csv['holo_chains'].values
        apo_chains = id_csv['apo_chains'].values
        holo_chains=holo_chains[0].split()
        apo_chains=apo_chains[0].split()
        sc_rmsd,tmscore,holo_bind_indices_all,apo_bind_indices_all,holo_bind_res,apo_bind_res=get_rmsds(direc,holo_chains,apo_chains)
        apo_file.write(direc+','+str(tmscore)+','+' '.join(holo_bind_res)+','+' '.join(map(str, holo_bind_indices_all))+','+' '.join(apo_bind_res)+','+' '.join(map(str, apo_bind_indices_all))+','+str(sc_rmsd)+'\n')