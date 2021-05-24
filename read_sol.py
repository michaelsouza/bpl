import pandas as pd
import numpy as np

def read_sol(fsol):
    df = {}
    colname = {}
    with open(fsol, 'r') as fid:
        for k, row in enumerate(fid.readlines()):
            if k == 0: # header
                for i, col in enumerate(row.split()):
                    df[col] = []
                    colname[i] = col
                continue
            for i, s in enumerate(row.split()):
                df[colname[i]].append(s)
    V = []
    for k in range(len(df['X'])):
        x = [float(df['X'][k]),float(df['Y'][k]),float(df['Z'][k])]
        V.append(x)
    df = pd.DataFrame(df)    
    V = np.array(V)
    return V, df

if __name__ == "__main__":
    fsol = 'DATA_LOOP_04/1egu.sol'
    V,_ = read_sol(fsol)    
