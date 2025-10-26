import os,re

def findAllPDBFile(base):
    for root, ds, fs in os.walk(base):
        for f in fs:
            if re.match(r'.*.pdb', f):
                # fullname = os.path.join(root, f)
                fullname = f
                yield fullname
    
if __name__ == "__main__":

    base = '/The path of the PDB files/' # change me
    for i in findAllPDBFile(base):
       for j in findAllPDBFile(base):
          print i+"_"+j
          ouputfile = '/output path/'+i+'_'+j+'.txt'
          commandline = "USalign.exe " +base+i +" " + base+j + ">" + ouputfile
          os.system(commandline)
