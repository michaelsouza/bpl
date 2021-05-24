import sys
from read_sol import read_sol
from numpy.linalg import norm

V, df = read_sol(sys.argv[1])
print(df)
print('|V[0]-V[4]| = %.12g' % norm(V[0] - V[4]))
for i in range(len(V)):
    row = '% 8.3f % 8.3f % 8.3f' % (V[i][0],V[i][1],V[i][2],)
    print('%2d' % i, row)