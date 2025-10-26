import os
import csv
from ImmuneBuilder import ABodyBuilder2
predictor = ABodyBuilder2()

if __name__ == "__main__":

    abs_seq_info = "../data/Sequences_of_antibodies_against_HA.csv"
    output_dir = "/output path/" # change me
    with open(abs_seq_info) as f:
        textfile=csv.reader(f)
        for row in textfile:
            pred_pdb = row[0]
            print(pred_pdb)
            pred_pdb = output_dir + pred_pdb + ".pdb"
            sequences = {
                "H": row[1].upper(),
                "L": row[2].upper() }
            
            antibody = predictor.predict(sequences)
            antibody.save(pred_pdb)
            print("end")