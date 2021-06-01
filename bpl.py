import sys
import numpy as np
from numpy.linalg import norm, solve
import time
from rmsd import rmsd
from read_sol import read_sol

def read_dat(fdat):
    print('Reading ' + fdat)
    df = {'ANCHOR':[], 'CONSTR':[]}
    with open(fdat, 'r') as fid:
        for row in fid:
            if 'ANCHOR' in row:
                row = row.split()
                x = np.array([float(row[-3]), float(row[-2]), float(row[-1])], dtype=float)
                i = int(row[-4]) - 1
                df['ANCHOR'].append({'i': i, 'x': x})
            elif 'INTERV' in row:
                row = row.split()
                i = int(row[-4]) - 1
                j = int(row[-3]) - 1
                lij = float(row[-2])
                uij = float(row[-1])
                df['INTERV'] = {'i': i, 'j': j, 'lij': lij, 'uij': uij}
            elif 'DATA' in row:
                row = row.split()
                ires = row[-8]
                jres = row[-7]
                iatm = row[-6]
                jatm = row[-5]
                i = int(row[-4]) - 1
                j = int(row[-3]) - 1
                lij = float(row[-2])
                uij = float(row[-1])
                df['CONSTR'].append({'i': i, 'j': j, 'lij': lij, 'uij': uij,
                                'ires': ires, 'jres': jres, 'iatm': iatm, 'jatm': jatm})
        df['ANCHOR'] = sorted(df['ANCHOR'], key=lambda u: u['i'])
        df['CONSTR'] = sorted(df['CONSTR'], key=lambda u: (u['i'], u['j']))
    return df 


def cleanD(p, D):
    S = {}
    q = np.zeros(len(p), dtype=int)
    for k, i in enumerate(p):
        q[i] = k
        S[i] = {}

    for i in range(len(p)):
        p_i = p[i]
        for j in D[p_i]:
            # q[j] < q
            if q[j] < i:
                S[p_i][j] = D[p_i][j]    
    return S


def order(df):
    # map atom type to each residue index
    RESIDUES = {}
    for c in df['CONSTR']:
        if c['ires'] not in RESIDUES:
            RESIDUES[c['ires']] = {}
        if c['jres'] not in RESIDUES:
            RESIDUES[c['jres']] = {}
        RESIDUES[c['ires']][c['iatm']] = c['i']
        RESIDUES[c['jres']][c['jatm']] = c['j']
    
    # anchors and first N
    p = [u['i'] for u in df['ANCHOR']]
    p.append(df['INTERV']['i'])

    # append remaining in N, CA, C order
    RESIDUES = [RESIDUES[k] for k in RESIDUES]
    for residue in RESIDUES:
        if residue['N'] not in p:
            p.append(residue['N'])
        if residue['CA'] not in p:
            p.append(residue['CA'])
        if residue['C'] not in p:
            p.append(residue['C'])
    return p


def constraints(df):
    lbnd = df['INTERV']['lij']
    ubnd = df['INTERV']['uij']
    D = {}
    for c in df['CONSTR']:
        i = c['i']
        j = c['j']
        dij = c['lij']
        if i not in D:
            D[i] = {}
        if j not in D:
            D[j] = {}
        D[i][j] = dij
        D[j][i] = dij
    return lbnd, ubnd, D


def solveEQ3(a, b, c, da, db, dc, dtol=1e-2, stol=1e-4):
    u = b - a
    A11 = norm(u)
    v = c - a
    A22 = norm(v)
    u = u / A11
    v = v / A22
    w = np.cross(u, v) # w perp u, v
    w = w / norm(w)
    uv = np.inner(u, v)
    A12 = A11 * uv
    A21 = A22 * uv
    # Let y = x - a, then x = y + a and y = y0*u + y1*v + y2*w
    # Using the constraints, we get
    # ||x - a|| = ||y|| = da
    # ||x - b|| = ||y - (b - a)|| = db
    # ||x - c|| = ||y - (c - a)|| = dc
    # Subtrating the square of the first from the one of two last equations, we have    
    A = [[A11, A12], [A21, A22]]
    B = [(da**2 - db**2 + A11**2)/2.0, (da**2 - dc**2 + A22**2)/2.0]
    y0, y1 = solve(A, B)
    s = da**2 - y0**2 - y1**2 - 2.0 * y0 * y1 * uv
    if s < 0 and np.abs(s) < stol:
        #print('Warning: base is almost plane (s=%g)' % s)
        s = 0
    if s < 0: # there is no solution
        # print('solveEQ3:: there is no solution (s=%g)' % s)
        return False, None, None
    
    proj_x = a + y0 * u + y1 * v # proj on the plane(a,b,c)
    
    y2 = np.sqrt(s)

    xpos = proj_x + y2 * w
    DA = norm(a - xpos)
    DB = norm(b - xpos)
    DC = norm(c - xpos)
    eij = np.max([np.abs(DA - da), np.abs(DB - db), np.abs(DC - dc)])
    if eij > dtol:
        raise Exception('xpos is not correct located')

    xneg = proj_x - y2 * w
    DA = norm(a - xneg)
    DB = norm(b - xneg)
    DC = norm(c - xneg)
    eij = np.max([np.abs(DA - da), np.abs(DB - db), np.abs(DC - dc)])
    if eij > dtol:
        raise Exception('xneg is not correct located')    
    
    return True, xpos, xneg


def viable(i, u, x, D, dtol=1e-2):
    for j in D[i]:
        dij = D[i][j]
        xj = x[j]
        DIJ = norm(u - xj)
        eij = np.abs(dij - DIJ)
        if eij > dtol:
            # print('unfeasible eij=%g' % eij)
            return False
    return True


def checkBinarySolution(s, p, x, D, dtol=1e-2):
    for i in range(4, len(s)):
        V = [j for j in D[p[i]]]
        a = x[V[0]]
        b = x[V[1]]
        c = x[V[2]]
        da = D[p[i]][V[0]]
        db = D[p[i]][V[1]]
        dc = D[p[i]][V[2]]
        _, xpos, xneg = solveEQ3(a, b, c, da, db, dc)
        if s[i] == 0:
            x[p[i]] = xpos
        else:
            x[p[i]] = xneg
    return errmax(p, x, D) < dtol
        

def errmax(p, x, D):
    eijMax = 0
    for k in range(4, len(p)):
        i = p[k]
        xi = x[i]
        for j in D[i]:
            dij = D[i][j]
            xj = x[j]
            DIJ = norm(xi - xj)
            eij = np.abs(dij - DIJ)
            if eij > eijMax:
                eijMax = eij                
    return eijMax

def bp(p, i, x, D, s, df, fid):
    if i == len(p):
        df['nsols'] += 1
        d13 = norm(x[p[1]]-x[p[3]])
        eijMax = errmax(p, x, D)
        rval, _, _ = rmsd(x, df['xsol'])
        fid.write('d13=%.12f eijMax=%.3E rmsd=%.3E s=[ ' % (d13, eijMax, rval))
        for k in range(len(s)):
            fid.write('%d ' % s[k])
        fid.write(']\n')
        return
    V = [j for j in D[p[i]]] # adjacent antecessors
    a = x[V[0]]
    b = x[V[1]]
    c = x[V[2]]
    da = D[p[i]][V[0]]
    db = D[p[i]][V[1]]
    dc = D[p[i]][V[2]]
    solved, xpos, xneg = solveEQ3(a, b, c, da, db, dc)
    if not solved:
        return
    # try xpos solution
    if viable(p[i], xpos, x, D):
        s[i] = 0
        x[p[i]] = xpos
        bp(p, i+1, x, D, s, df, fid)
    # try xneg solution
    if viable(p[i], xneg, x, D):
        s[i] = 1
        x[p[i]] = xneg
        bp(p, i+1, x, D, s, df, fid)


def bpl(p, lbnd, ubnd, D, df, fid, num=100):    
    df['nsols'] = 0

    x = np.zeros((len(p), 3), dtype=float)

    # set base    
    x[p[0]] = df['ANCHOR'][0]['x']
    x[p[1]] = df['ANCHOR'][1]['x']
    x[p[2]] = df['ANCHOR'][2]['x']

    # solution binary representation
    s = np.zeros(len(p), dtype=int) 

    d03 = D[p[3]][p[0]]
    d23 = D[p[3]][p[2]]
    D13 = [d13 for d13 in np.linspace(lbnd, ubnd, num=num, endpoint=True, dtype=float)]

    for d13 in D13:
        # only positive oriented branch
        solved, xpos, xneg = solveEQ3(x[p[0]], x[p[1]], x[p[2]], d03, d13, d23)
        if not solved:
            raise Exception('ibpl::the first system could not be solved')
        
        x[p[3]] = xpos
        s[3] = 0
        bp(p, 4, x, D, s, df, fid)

        # TODO Consider to use only one root (first branch)
        x[p[3]] = xneg
        s[3] = 1
        bp(p, 4, x, D, s, df, fid)


if __name__ == '__main__':    
    # default parameters
    num = 100
    fdat = 'DATA_LOOP_08/1i0h.dat'

    # read input
    if len(sys.argv) > 1:
        fdat = sys.argv[1]
    if len(sys.argv) > 2:
        num = int(sys.argv[2])

    df = read_dat(fdat)
    df['xsol'], _ = read_sol(fdat.replace('.dat', '.sol'))
    
    lbnd, ubnd, D = constraints(df)
    p = order(df)
    D = cleanD(p, D)
    
    flog = fdat.replace('.dat', '_%d.log' % num)
    fid = open(flog, 'w')

    tic = time.time()
    bpl(p, lbnd, ubnd, D, df, fid, num)
    toc = time.time() - tic

    nedges = 0
    for i in D:
        nedges += len(D[i])

    fid.write('num ..... %d\n' % num)
    fid.write('nnodes .. %d\n' % len(p))
    fid.write('nedges .. %d\n' % nedges)
    fid.write('nsols ... %d\n' % df['nsols'])    
    fid.write('tsecs ... %f\n' % toc)
    
    fid.close()
