# USAGE: python call_ibpl.py | tee call_ibpl.log
import os

WDIR = ['DATA_LOOP_04', 'DATA_LOOP_08', 'DATA_LOOP_12']
# NUM = [10, 100, 1000]
NUM = [10]
if __name__ == "__main__":
    njobs = 4
    cmd = 'parallel -j %d -- ' % njobs
    for wdir in WDIR:
        for num in NUM:
            for fdat in os.listdir(wdir):
                if not fdat.endswith('.dat'):
                    continue
                fdat = os.path.join(wdir, fdat)
                cmd += '"python ibpl.py %s %d" ' % (fdat, num)
    print(cmd)
    os.system(cmd)
