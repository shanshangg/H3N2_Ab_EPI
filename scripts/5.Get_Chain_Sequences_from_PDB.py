import csv
from Bio import SeqIO

if __name__ == "__main__":
    chain_info = "../data/Hemagglutinin_related_sabdab_summary.tsv"
    outpath = open("chains_type_and_len.csv",'w')
    print("PDB_id,chain_id,chain_type,chain_lenth\n",end='',file=outpath)
    outpath3 = open("ChainSequence_ab_H.txt",'w')
    outpath4 = open("ChainSequence_ab_L.txt",'w')
    outpath2 = open("ChainSequence_ag.txt",'w')
    
    for fa in SeqIO.parse("ChainSequence.txt", "fasta"):
        pdb_id = fa.id[:4]
        chain_id = fa.description[-1]
        len_chain = len(fa.seq)
        seq = fa.seq
        type_chain = "ini"
        print(pdb_id+'_'+chain_id+'_'+str(len_chain))
        with open(chain_info) as f:
            chain_info_list=csv.reader(f,dialect=csv.excel_tab)
            for item in chain_info_list:
                if pdb_id == item[0]: 
                    ag_chains = item[4].replace(" | ","") 
                    
                    if chain_id == item[1]:
                        type_chain = "ab_H"
                        print('>'+pdb_id+'_'+chain_id,file=outpath3)
                        print(seq,file=outpath3)
                    elif  chain_id == item[2]:
                        type_chain ="ab_L"
                        print('>'+pdb_id+'_'+chain_id,file=outpath4)
                        print(seq,file=outpath4)
                    elif chain_id in ag_chains:
                        type_chain ="ag"
                        print('>'+pdb_id+'_'+chain_id,file=outpath2)
                        print(seq,file=outpath2)
                    else:
                        type_chain ="na"
                    break 

            print(pdb_id+'_,'+chain_id+','+type_chain+','+str(len_chain)+'\n',end='',file=outpath) 
            f.close  
    outpath.close()
    outpath2.close()
    print("end")