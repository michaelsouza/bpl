import os
import pandas as pd
from read_pdb import read_pdb

loopDefs = pd.read_excel('table_loopDefs.xlsx')

PDB = {}
ATOMS = ['C', 'CA', 'N']

# reading pdb data
for fpdb in os.listdir('pdb'):
    if not fpdb.endswith('.pdb.gz'):
        continue
    pdbCode = fpdb.split('.')[0]
    fpdb = os.path.join('pdb', fpdb)
    pdbDataFrame = read_pdb(fpdb)
    PDB[pdbCode] = pdbDataFrame

# creating instances
for i, loop in loopDefs.iterrows():
    if loop.pdbCode not in PDB:
        print('PDB not found (%s)' % loop.pdbCode)
        continue
    df = PDB[loop.pdbCode]
    df = df[loop.resSeqFirst <= df['resSeq']]
    df = df[loop.resSeqLast  >= df['resSeq']]
    df = df.query("chainID=='%s'" % loop.chainID)
    df = df.query("name=='C' | name=='CA' | name=='N'")   
    df = df.query("altLoc==' ' | altLoc=='A'")
    wdir = 'DATA_LOOP_%02d' % loop.loopSize
    if not os.path.isdir(wdir):
        os.mkdir(wdir)
    fcsv = os.path.join(wdir, loop.pdbCode + '.csv')
    df.to_csv(fcsv,index=False)
    
