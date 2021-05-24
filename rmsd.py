import numpy as np

def center(x):
    c = [0.0,0.0,0.0]
    for i in range(len(x)):
        c[0] += x[i,0]
        c[1] += x[i,1]
        c[2] += x[i,2]
    n = len(x)
    c[0] /= n
    c[1] /= n
    c[2] /= n
    return c

def translation(x, c):
    for i in range(len(x)):
        x[i] = x[i] - c
    return x

def rmsd(x, y):
    x = np.copy(x)
    y = np.copy(y)

    if x.shape != y.shape:
        raise Exception('x and y must have the same shape')
    
    # translation to origin
    x = translation(x, center(x))
    y = translation(y, center(y))

    M = np.transpose(y) @ x

    U,_,V = np.linalg.svd(M)

    Q = U @ V

    y = y @ Q

    return np.linalg.norm(x -y, ord='fro'), x, y


if __name__ == "__main__":
    natoms = 5
    x = [[2.523,  12.637,  -4.449],
         [2.074,  12.786,  -5.822],
         [3.231,  13.19,  -6.737],
         [7.698,  18.571,  -9.357],
         [8.726,  18.865, -10.341],
         [9.905,  19.599,  -9.7],
         [13.695,  26.437,  -9.645],
         [15.02,  26.992,  -9.428],
         [14.963,  28.2,  -8.492],
         [16.687,  32.163,  -9.003],
         [17.264,  33.262,  -9.758],
         [18.27,  32.75, -10.79]]
    x = np.array(x)

    # add some noise and reflection
    y = x + 5
    y = -y
    print('|x-y|=%g' % np.linalg.norm(x-y,'fro'))

    r, x, y = rmsd(x, y)

    print('|x-y|=%g' % np.linalg.norm(x-y,'fro'))
    
