import pandas as pd
import os
import sys
import numpy as np


def distance(xi, xj):
    return np.sqrt((xi[0] - xj[0])**2 + (xi[1] - xj[1])**2 + (xi[2] - xj[2])**2)


def distance_range(i, a, b, c, V):
    # range of distance between i and c
    X = V[i]['XYZ']
    A = V[a]['XYZ']
    B = V[b]['XYZ']
    C = V[c]['XYZ']
    dab = distance(A, B)
    dac = distance(A, C)
    dbc = distance(B, C)
    dax = distance(A, X)
    dbx = distance(B, X)
    Cx = (dab**2 + dac**2 - dbc**2) / (2 * dab)
    Cy = np.sqrt(dac**2 - Cx**2)
    A = np.array([0, 0, 0], dtype=float)
    B = np.array([dab, 0,  0], dtype=float)
    C = np.array([Cx, Cy, 0], dtype=float)
    dAB = distance(A, B)
    dAC = distance(A, C)
    dBC = distance(B, C)
    if np.max([np.abs(dAB - dab), np.abs(dAC - dac), np.abs(dBC - dbc)]) > 1e-8:
        raise Exception('A,B,C are not correct located')
    Xx = (dab**2 + dax**2 - dbx**2) / (2*dab)
    Xy = np.sqrt(dax**2 - Xx**2)
    Xmin = np.array([Xx, Xy, 0], dtype=float)
    dAX = distance(A, Xmin)
    dBX = distance(B, Xmin)
    if np.max([np.abs(dAX - dax), np.abs(dBX - dbx)]) > 1e-8:
        raise Exception('Xmin is not correct located')
    Xmax = np.array([Xx, -Xy, 0], dtype=float)

    lij = distance(C, Xmin) + 1E-3
    uij = distance(C, Xmax) - 1E-3
    return lij, uij


def checkResidues(RESIDUES):
    for resSeq in RESIDUES:
        if len(RESIDUES[resSeq]) < 3:
            print(RESIDUES[resSeq])
            raise Exception('The %d-th residue is not well defined.' % resSeq)
    return True


def readResidues(df):
    RESIDUES = {}

    # group N, C, CA of each residue
    for _, row in df.iterrows():
        resSeq = row['resSeq']
        if resSeq not in RESIDUES:
            RESIDUES[resSeq] = {}
        atom = {'nid': None, 'XYZ': (
            row.x, row.y, row.z), 'resSeq': resSeq, 'name': row['name']}
        RESIDUES[resSeq][row['name']] = atom

    checkResidues(RESIDUES)

    # convert from dict to list of RESIDUES
    nid = 0  # node id
    for resSeq in sorted(RESIDUES):
        RESIDUES[resSeq]['resSeq'] = resSeq
        RESIDUES[resSeq]['N']['nid'] = nid
        RESIDUES[resSeq]['CA']['nid'] = nid + 1
        RESIDUES[resSeq]['C']['nid'] = nid + 2
        nid += 3
    RESIDUES = [RESIDUES[resSeq] for resSeq in sorted(RESIDUES)]

    return RESIDUES


def writeRIGID(fid, k, R, V, XYZ):
    row = '\nINFO: COLUMN RES_I RES_J  ATM_I ATM_J IDX_I IDX_J         LIJ             UIJ'
    print(row)
    fid.write(row + '\n')
    for i in range(len(R)):
        for j in range(i+1, len(R)):
            xi = XYZ[R[i]]
            xj = XYZ[R[j]]
            vi = V[R[i]]
            vj = V[R[j]]
            dij = distance(xi, xj)
            row = 'DATA: RIGID%d %5d %5d' % (k, vi['resSeq'], vj['resSeq'])
            row += '%7s %5s' % (vi['name'], vj['name'])
            row += '%6d %5d' % (vi['nid'] + 1, vj['nid'] + 1)
            row += '%16.8f' % dij
            row += '%16.8f' % dij
            print(row)
            fid.write(row + '\n')


def writeLNK(fid, R, V, i, j):
    xi = XYZ[i]
    xj = XYZ[j]
    vi = V[i]
    vj = V[j]
    row = 'DATA: LINKS  %5d %5d' % (vi['resSeq'], vj['resSeq'])
    row += '%7s %5s' % (vi['name'], vj['name'])
    row += '%6d %5d' % (vi['nid'] + 1, vj['nid'] + 1)
    dij = distance(xi, xj)
    row += '%8.5f' % dij
    row += '%8.5f' % dij
    print(row)    
    fid.write(row + '\n')


if __name__ == "__main__":
    folder = 'DATA_LOOP_08'
    if len(sys.argv) > 1:
        folder = sys.argv[1]

    print('Creating loops from folder ' + folder)

    for fcsv in os.listdir(folder):
        if not fcsv.endswith('.csv') or '_' in fcsv:
            continue
        print('\n=================================================================')
        print('Reading ' + os.path.join(folder, fcsv))
        df = pd.read_csv(os.path.join(folder, fcsv))
        RESIDUES = readResidues(df)
        print('fcsv ' + os.path.join(folder, fcsv))
        fdat = os.path.join(folder, fcsv).replace('.csv','.dat')
        fid = open(fdat, 'w')

        V = []
        XYZ = [None for _ in range(3 * len(RESIDUES))]
        for residue in RESIDUES:
            XYZ[residue['N']['nid']] = residue['N']['XYZ']
            XYZ[residue['CA']['nid']] = residue['CA']['XYZ']
            XYZ[residue['C']['nid']] = residue['C']['XYZ']
            V.append(residue['N'])
            V.append(residue['CA'])
            V.append(residue['C'])
        
        L = int(int(len(RESIDUES) / 4.0))
        FIXED_RID = [0, L, 2*L]

        # anchors
        row = 'INFO: COLUMN RESSEQ ATOM  INDEX    X        Y         Z'
        print(row)
        fid.write(row + '\n')
        for rid in FIXED_RID:
            residue = RESIDUES[rid]
            row = 'DATA: ANCHOR %6d   CA %6d' % (
                residue['resSeq'], residue['CA']['nid'] + 1)
            row = row + (' %8.5f %8.5f %8.5f' % residue['CA']['XYZ'])
            print(row)
            fid.write(row + '\n')

        # INTERVAL CONSTRAINT
        i = RESIDUES[0]['N']['nid']
        j = RESIDUES[FIXED_RID[1]]['CA']['nid']
        a = RESIDUES[FIXED_RID[0]]['CA']['nid']
        b = RESIDUES[FIXED_RID[2]]['CA']['nid']
        lij, uij = distance_range(i, a, b, j, V)
        row = '\nINFO: COLUMN RES_I RES_J  ATM_I ATM_J IDX_I IDX_J         LIJ             UIJ'
        print(row)
        fid.write(row + '\n')
        vi, vj = V[i], V[j]
        row = 'DATA: INTERV %5d %5d' % (vi['resSeq'], vj['resSeq'])
        row += '%7s %5s' % (vi['name'], vj['name'])
        row += '%6d %5d' % (vi['nid'] + 1, vj['nid'] + 1)
        row += '%16.8f' % lij
        row += '%16.8f' % uij        
        print(row)
        fid.write(row + '\n')

        # first rigid body
        R = [RESIDUES[0]['CA']['nid']]
        R.append(RESIDUES[0]['C']['nid'])
        for rid in range(FIXED_RID[0]+1, FIXED_RID[1]):
            R.append(RESIDUES[rid]['N']['nid'])
            R.append(RESIDUES[rid]['CA']['nid'])
            R.append(RESIDUES[rid]['C']['nid'])
        R.append(RESIDUES[FIXED_RID[1]]['N']['nid'])
        R.append(RESIDUES[FIXED_RID[1]]['CA']['nid'])
        writeRIGID(fid, 0, R, V, XYZ)

        # second rigid body
        R = [RESIDUES[FIXED_RID[1]]['CA']['nid']]
        R.append(RESIDUES[FIXED_RID[1]]['C']['nid'])
        for rid in range(FIXED_RID[1]+1, FIXED_RID[2]):
            R.append(RESIDUES[rid]['N']['nid'])
            R.append(RESIDUES[rid]['CA']['nid'])
            R.append(RESIDUES[rid]['C']['nid'])
        R.append(RESIDUES[FIXED_RID[2]]['N']['nid'])
        R.append(RESIDUES[FIXED_RID[2]]['CA']['nid'])
        writeRIGID(fid, 1, R, V, XYZ)

        # third rigid body
        R = [RESIDUES[FIXED_RID[2]]['CA']['nid']]
        R.append(RESIDUES[FIXED_RID[2]]['C']['nid'])
        for rid in range(FIXED_RID[2]+1, len(RESIDUES)):
            R.append(RESIDUES[rid]['N']['nid'])
            R.append(RESIDUES[rid]['CA']['nid'])
            R.append(RESIDUES[rid]['C']['nid'])
        R.append(RESIDUES[0]['N']['nid'])
        R.append(RESIDUES[0]['CA']['nid'])
        writeRIGID(fid, 2, R, V, XYZ)

        # additional constraints N-C
        row = '\nINFO: COLUMN RES_I RES_J  ATM_I ATM_J IDX_I IDX_J   LIJ     UIJ'        
        print(row)
        fid.write(row + '\n')
        i = RESIDUES[0]['N']['nid']  # last N of the first rigid body
        j = RESIDUES[0]['C']['nid']  # first C of the second rigid body
        writeLNK(fid, R, V, i, j)
        # last N of the second rigid body
        i = RESIDUES[FIXED_RID[1]]['N']['nid']
        # first C of the second rigid body
        j = RESIDUES[FIXED_RID[1]]['C']['nid']
        writeLNK(fid, R, V, i, j)
        # last N of the second rigid body
        i = RESIDUES[FIXED_RID[2]]['N']['nid']
        # first C of the second rigid body
        j = RESIDUES[FIXED_RID[2]]['C']['nid']
        writeLNK(fid, R, V, i, j)

        fid.close()

        fsol = fdat.replace('.dat','.sol')
        print('\nWriting ' + fsol)
        with open(fsol, 'w') as fid:
            row = 'RESID    NAME      NID     X        Y        Z'        
            print(row)
            fid.write(row + '\n')
            for v in V:
                row = '%5d %7s %8d ' % (v['resSeq'], v['name'],  v['nid'] + 1)
                row += '%8.5f %8.5f %8.5f' % v['XYZ']
                print(row)
                fid.write(row + '\n')