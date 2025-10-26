from pymol import cmd
import os,re

class PDBparase():
    def __init__(self,pdbid,pdbfilepath,out_path=os.path.abspath('.')):
        self.pdbid = pdbid
        self.path = pdbfilepath
        self.outpath = open(os.path.join(out_path,'Bindsite.txt'),'a')
        self.outpath2 = open(os.path.join(out_path,'ChainSequence.txt'),'a')
        self.aadic = {
         'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E',
         'PHE':'F', 'GLY':'G', 'HIS':'H', 'LYS':'K',
         'ILE':'I', 'LEU':'L', 'MET':'M', 'ASN':'N',
         'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S',
         'THR':'T', 'VAL':'V', 'TYR':'Y', 'TRP':'W'}

    def extractChain(self):
        '''
        extract pdb file Chain list
        :return:
        '''
        cmd.delete('all')
        cmd.load(self.path)
        cmd.remove('solvent')
        cmd.select('all_','all')
        prot_atom = cmd.get_model('all_').atom
        chain_prot = [ atom.chain for atom in prot_atom]
        chain_lig = sorted(list(set(chain_prot)))
        return chain_lig

    def extractinterchain(self):
        chain_ls = self.extractChain()
        prot_chain = []
        for c in chain_ls:
            cmd.select('singleseq', '  (byres polymer & name CA) in chain %s' % c)
            singleseq_atom = cmd.get_model('singleseq').atom
            seq = ''
            aanum = []
            for a in singleseq_atom:
                if a.resn in self.aadic:
                    if a.resi not in aanum:
                        seq += self.aadic[a.resn]
                        aanum.append(a.resi)
            if aanum:
                print('>%s\t%s\n%s'%(self.pdbid,c,seq), file=self.outpath2)
                print('>%s\t%s\n%s'%(self.pdbid,c,seq))
                prot_chain.append(c)
            cmd.delete('singleseq')
        return prot_chain

    def calpip(self,cutoff=3.5):
        inter_chain = self.extractinterchain()
        for c in inter_chain:
            cmd.select('rcpseq', '  (byres polymer & name CA) in chain %s' % c)
            for d in inter_chain:
                if d != c:
                    cmd.select('interpip', '  (byres polymer & name CA) in chain %s' % d)
                    cmd.select('protsite',' byres interpip around %s in rcpseq'%cutoff)
                    bsnum_ls = []
                    site_ls = []
                    for a in cmd.get_model('protsite').atom:
                        if a.resn in self.aadic:
                            aazong = self.aadic[a.resn] + a.resi
                            if a.resi not in bsnum_ls:
                                bsnum_ls.append(a.resi)
                            if aazong not in site_ls:
                                site_ls.append(aazong)
                    seq_atom = cmd.get_model('protsite').atom
                    seq = ''
                    aanum = []
                    for a in seq_atom:
                        if a.resn in self.aadic:
                            if a.resi not in aanum:
                                if a.resi in bsnum_ls:
                                    seq += self.aadic[a.resn].lower()
                                else:
                                    seq += self.aadic[a.resn]
                                aanum.append(a.resi)
                    if site_ls:
                        print('%s\t%s\t%s\t%s\t%s\t%s' % (self.pdbid, d,'PIP', c, ' '.join(site_ls), seq),
                              file=self.outpath)
                        print('%s\t%s\t%s\t%s\t%s\t%s' % (self.pdbid, d,'PIP', c, ' '.join(site_ls), seq))
                        
def findAllPDBFile(base):
    for root, ds, fs in os.walk(base):
        for f in fs:
            if re.match(r'.*.pdb', f):
                fullname = f
                yield fullname

if __name__ == '__main__':

    base = '/The path of the PDB files/' # change me
    for i in findAllPDBFile(base):
        ii = i[:-4]
        print(ii)
        pdb = PDBparase(ii,i)
        pdb.extractinterchain()
        pdb.calpip()
       
    pdb.outpath.close()
    pdb.outpath2.close()
