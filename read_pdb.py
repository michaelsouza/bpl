import gzip
import random
import pandas as pd

def read_pdb(fn_pdb, EPSX=1E-7):
    print('Reading %s' % fn_pdb)
    atoms = {'RECORD':[], 'serial': [], 'name': [], 'altLoc': [], 'resName': [], 'chainID': [],
             'resSeq': [], 'iCode': [], 'x': [], 'y': [], 'z': [], 'occupancy': [],
             'tempFactor': [], 'element': [], 'charge': []}
    with gzip.open(fn_pdb, 'r') as fid:
        for row in fid:
            row = row.decode('utf-8')            
            if row.startswith('ATOM') or row.startswith('HETATM'):
                atoms['RECORD'].append(row[:6].strip())
                atoms['serial'].append(int(row[6:11]))
                atoms['name'].append(row[12:16].strip())
                atoms['altLoc'].append(row[16])
                atoms['resName'].append(row[17:20].strip())
                atoms['chainID'].append(row[21])
                atoms['resSeq'].append(int(row[22:26]))
                atoms['iCode'].append(row[26])
                # add random noise to avoid colinearity
                atoms['x'].append(float(row[30:38]) + random.uniform(-EPSX, EPSX))
                atoms['y'].append(float(row[38:46]) + random.uniform(-EPSX, EPSX))
                atoms['z'].append(float(row[46:54]) + random.uniform(-EPSX, EPSX))
                atoms['occupancy'].append(float(row[54:60]))
                atoms['tempFactor'].append(float(row[60:66]))
                atoms['element'].append(row[76:78].strip())
                atoms['charge'].append(row[78:80].strip())
            if 'ENDMDL' in row:  # read only the first model
                break    
    return pd.DataFrame.from_dict(atoms)

if __name__ == "__main__":
    import sys
    df = read_pdb(sys.argv[1])
    fcsv = sys.argv[1].replace('.pdb.gz','.csv')
    print('Writing ' + fcsv)
    df.to_csv(fcsv, index=False)