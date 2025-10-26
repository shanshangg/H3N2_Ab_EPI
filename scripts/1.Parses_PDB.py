import os
import csv
from Bio import PDB

class ChainSplitter:
    def __init__(self, out_dir=None):
        self.parser = PDB.PDBParser()
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "chain_PDBs")
        self.out_dir = out_dir
        
    def make_pdb(self, pdb_path, chain_letters_ab, chain_letters_ag, overwrite=True, struct=None):
        chain_letters = chain_letters_ab + chain_letters_ag
        (pdb_dir, pdb_fn) = os.path.split(pdb_path)
        out_name = "%s_%s_%s.pdb" % (pdb_id, "".join(chain_letters_ab),"".join(chain_letters_ag))
        out_path = os.path.join(self.out_dir, out_name)
        print ("OUT PATH:",out_path)
        plural = "s" if (len(chain_letters) > 1) else ""
        if (not overwrite) and (os.path.isfile(out_path)):
            print("Chain%s %s of '%s' already extracted to '%s'." %
                    (plural, ", ".join(chain_letters), pdb_id, out_name))
            return out_path
        print("Extracting chain%s %s from %s..." % (plural,
                ", ".join(chain_letters), pdb_fn))
        if struct is None:
            struct = self.parser.get_structure(pdb_id, pdb_path)
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=SelectChains(chain_letters))
        return out_path

class SelectChains(PDB.Select):
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)

if __name__ == "__main__":

    pdb_textfn = "../data/Hemagglutinin_related_sabdab_summary.tsv"
    splitter = ChainSplitter("/output path/")  # Change me.
    with open(pdb_textfn) as f:
        pdb_textfile=csv.reader(f,dialect=csv.excel_tab)
        for row in pdb_textfile:
            pdb_id = row[0].lower()
            chain_ag = row[4]
            chain_ag = chain_ag.replace(" | ","") 
            if row[2] == "NA":
                chain_ab = row[1]
            else:
                chain_ab = row[1]+row[2]
            print(pdb_id+"_"+chain_ab+"_"+chain_ag)
            splitter.make_pdb(pdb_id+".pdb", chain_ab, chain_ag)
    print("end")