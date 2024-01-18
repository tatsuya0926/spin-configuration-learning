import numpy as np
import os
import sys
from tqdm import tqdm


def convert_2d_Ising_txt2npy():
    NX = 64  # number of sites along x-direction
    NY = 64  # number of sites along y-direction
    nconf = 31
    ndata = 50
    for t in tqdm(range(nconf)):
        for i in range(ndata):
            if (os.path.exists(f'txtfile/2d_Ising/L{NX}T{t}_{i}.txt')):
                spin = np.empty((NX, NY))
                for read in open(f'txtfile/2d_Ising/L{NX}T{t}_{i}.txt').readlines():
                    read = read[:-2].split(' ')
                    ix = int(read[0])
                    iy = int(read[1])
                    spin[ix, iy] = read[2]
                    np.save(f'dataset/2d_Ising/L{NX}/L{NX}T{t}_{i}.npy', spin)
            else:
                print('no input configuration')
                sys.exit()


def convert_2d_Potts_txt2npy():
    NX = 64  # number of sites along x-direction
    NY = 64  # number of sites along y-direction
    nconf = 31
    ndata = 25
    for t in tqdm(range(nconf)):
        for i in range(ndata):
            if (os.path.exists(f'txtfile/2d_Potts/L{NX}T{t}_{i}.txt')):
                spin = np.empty((NX, NY))
                for read in open(f'txtfile/2d_Potts/L{NX}T{t}_{i}.txt').readlines():
                    read = read[:-2].split(' ')
                    ix = int(read[0])
                    iy = int(read[1])
                    spin[ix, iy] = read[2]
                    np.save(f'dataset/2d_Potts/L{NX}/L{NX}T{t}_{i}.npy', spin)
            else:
                print('no input configuration')
                sys.exit()


def convert_2d_Clock_txt2npy():
    NX = 64  # number of sites along x-direction
    NY = 64  # number of sites along y-direction
    nconf = 51
    ndata = 25
    Q = 4
    for t in tqdm(range(nconf)):
        for i in range(ndata):
            if (os.path.exists(f'txtfile/2d_Clock/L{NX}T{t}_{i}.txt')):
                spin = np.empty((NX, NY))
                for read in open(f'txtfile/2d_Clock/L{NX}T{t}_{i}.txt').readlines():
                    read = read[:-2].split(' ')
                    ix = int(read[0])
                    iy = int(read[1])
                    spin[ix, iy] = read[2]
                    np.save(
                        f'dataset/2d_Clock/L{NX}_q={Q}/L{NX}T{t}_{i}.npy', spin)
            else:
                print('no input configuration')
                sys.exit()


if __name__ == '__main__':
    # convert_2d_Potts_txt2npy()
    convert_2d_Clock_txt2npy()
