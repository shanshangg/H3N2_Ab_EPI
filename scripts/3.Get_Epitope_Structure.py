"""
This script is adapted from the original script `extract_epitope_from_complex.py` 
available at the CE_BLAST GitHub repository: 
https://github.com/baddtongji/CE_BLAST

All dependencies required for this script are also obtained from the same repository.
"""

import os, re, itertools
from dist import atom_dist
from Bio.PDB import Select
from basic_pdb_related_func import preprocessFile, getPdbStruct, writePdbStruct2File

class ChainSelect(Select):
    def __init__(self, chn_lst):
        self.chn_lst = chn_lst
    def accept_chain(self, chain):
        if chain.get_id() in self.chn_lst:
            return 1
        else:
            return 0

class ResiSelect(Select):
    def __init__(self, resi_lst):
        self.resi_lst = resi_lst
    def accept_residue(self, residue):
        res_id = str(residue.get_id()[1])
        res_uni_id = res_id
        if res_uni_id in self.resi_lst:
            return 1
        else:
            return 0

def creatAbAgPdb(pdb_path, ag_chn, ab_chn_lst):
    pdb_path = preprocessFile(pdb_path)
    pdb_doc = os.path.dirname(pdb_path)
    pdb_basename = os.path.basename(pdb_path)
    struct = getPdbStruct(pdb_path)   
    ag_chn_lst = [ag_chn]
    if ab_chn_lst is None:
        ab_chn_lst = [chn.get_id() for chn in struct.get_chains() if chn.get_id()!=ag_chn]
    else:
        ab_chn_lst = ab_chn_lst
    ag_path = os.path.join(pdb_doc, pdb_basename.split('.')[0] + '_' + '_'.join(ag_chn_lst) + '.temp')
    ab_path = os.path.join(pdb_doc, pdb_basename.split('.')[0] + '_' + '_'.join(ab_chn_lst) + '.temp')
    getChainFromFile(pdb_path, ag_chn_lst, ag_path)
    getChainFromFile(pdb_path, ab_chn_lst, ab_path)
    return ag_path, ab_path

def getChainFromFile(pdb_path, chn_lst, out_path):
    struct = getPdbStruct(pdb_path)
    if type(chn_lst) != list:
        chn_lst = list(chn_lst)
    writePdbStruct2File(struct, out_path, ChainSelect(chn_lst))

def getAgEpitopeResidue(ag_path, ab_path, cutoff=4):
    ag_epitope_resi_lst = []
    ag_struct = getPdbStruct(ag_path)
    ab_struct = getPdbStruct(ab_path)
    agab_resi_product = itertools.product(ag_struct.get_residues(), ab_struct.get_residues())
    for ag_res, ab_res in agab_resi_product:
        atom_pair_product = itertools.product(ag_res.get_list(), ab_res.get_list())
        atom_pair_record = 0
        for atom1, atom2 in atom_pair_product:
            if atom_dist(atom1, atom2) <= cutoff:
                atom_pair_record = 1
                break
        if atom_pair_record == 1:
            res_id = str(ag_res.get_id()[1])
            res_uni_id = res_id
            ag_epitope_resi_lst.append(res_uni_id)
    return ag_epitope_resi_lst

def createAgEpitopeFile(ag_path, ag_epi_resi_lst, ag_epitope_path):
    struct = getPdbStruct(ag_path)
    writePdbStruct2File(struct, ag_epitope_path, ResiSelect(ag_epi_resi_lst))
    
def extractEpitopeFromComplex(pdb_path, ag_chn, epitope_path, ab_chn_lst=None, cutoff=5):
    ag_path, ab_path = creatAbAgPdb(pdb_path, ag_chn, ab_chn_lst=ab_chn_lst)
    ag_epi_resi_lst = getAgEpitopeResidue(ag_path, ab_path, cutoff=cutoff)
    createAgEpitopeFile(ag_path, ag_epi_resi_lst, epitope_path)
    os.remove(ag_path)
    os.remove(ab_path)

def findAllPDBFile(base):
    for root, ds, fs in os.walk(base):
        for f in fs:
            if re.match(r'.*.pdb', f):
                fullname = f
                yield fullname

if __name__ == '__main__':

    base = '/The path of the PDB files/' # change me
    output = '/output path/' # change me
    for i in findAllPDBFile(base):
        ii = i[:-4]
        pdb_info = ii.split('_')
        print(ii)
        pdb_path = base+i
        ag_chn = pdb_info[2]
        epitope_path = output + pdb_info[0] + "_" + pdb_info[1] + ".pdb"
        extractEpitopeFromComplex(pdb_path, ag_chn, epitope_path)
        print 'Done!'
