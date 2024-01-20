import numpy as np
import os
import sys
from tqdm import tqdm
from multiprocessing import Pool


def txt2npy(inputs):
    t, i, L, model_name, Q = inputs
    if Q == None:
        import_path_name = f'txtfile/{model_name}/L{L}T{t}_{i}.txt'
        export_path_name = f'dataset/{model_name}/L{L}/L{L}T{t}_{i}.npy'
    else:
        import_path_name = f'txtfile/{model_name}/q={Q}/L{L}T{t}_{i}.txt'
        export_path_name = f'dataset/{model_name}/L{L}_q={Q}/L{L}T{t}_{i}.npy'
    if (os.path.exists(import_path_name)):
        spin = np.empty((L, L))
        for read in open(import_path_name).readlines():
            read = read[:-2].split(' ')
            ix = int(read[0])
            iy = int(read[1])
            spin[ix, iy] = read[2]
            np.save(export_path_name, spin)
    else:
        print('no input configuration')
        sys.exit()


def convert_txt2npy(L, nconf, ndata, model_name, Q=None):
    values = [(t, i, L, model_name, Q)
              for t in range(nconf) for i in range(ndata)]
    list(tqdm(Pool().imap(txt2npy, values), total=nconf*ndata))
    Pool().close()
    Pool().join()


if __name__ == '__main__':
    convert_txt2npy(L=64, nconf=31, ndata=50, model_name="2d_Ising")
    # convert_txt2npy(L=64, nconf=31, ndata=50, model_name="2d_Potts")
    # convert_txt2npy(L=64, nconf=31, ndata=25, model_name="2d_Clock", Q=2)
    # convert_txt2npy(L=64, nconf=31, ndata=25, model_name="2d_Clock", Q=4)
    # convert_txt2npy(L=64, nconf=31, ndata=25, model_name="2d_Clock", Q=5)
