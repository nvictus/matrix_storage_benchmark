from __future__ import division, print_function
from multiprocessing import Pool
import os.path as op
import random
import sys
import time
import argparse

import numpy as np
import h5py
import cooler
from cooler import ice


def main():
    parser = argparse.ArgumentParser(description="""

    python matrix_storage_benchmark.py matrix_tsv
""")

    parser.add_argument('matrix_tsv', nargs=1)
    parser.add_argument('-s', '--square', default=256,
                                        help="The size of the square within which to return values", type=int)
    parser.add_argument('-i', '--iterations', default=100,
                                        help="The number of times to run the range query", type=int)
    args = parser.parse_args()
    infilepath = args.matrix_tsv[0]
    outfilepath = op.join(op.dirname(infilepath), 'test.cool')


    # Fetch hg19 chromsizes from UCSC
    chromsizes = cooler.read_chromsizes('http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz')


    # Build "index"
    t1 = time.time()
    chromsizes = chromsizes.loc['chrX':'chrX']
    chroms = ['chrX']
    lengths = [chromsizes['chrX']]
    binsize = 5000
    bins = cooler.binnify(chromsizes, binsize)
    chunksize = 10000000
    reader = cooler.io.SparseLoader(infilepath, chunksize)

    with h5py.File(outfilepath, 'w') as h5:
        h5opts = dict(compression='gzip', compression_opts=6)
        cooler.io.create(h5, chroms, lengths, bins, reader, binsize, h5opts=h5opts)

    c = cooler.Cooler(outfilepath)
    matrix = c.matrix()
    print("Time creating index: {:.2f} seconds".format(time.time() - t1))


    # Normalization
    t15 = time.time()
    N_CPUS = 5
    chunksize = int(100e6)

    with h5py.File(outfilepath, 'r+') as h5, Pool(N_CPUS) as pool:
        bias = ice.iterative_correction(
            h5, chunksize=chunksize, tol=1e-05, min_nnz=10,
            cis_only=False, ignore_diags=3, map=pool.map)

        h5opts = dict(compression='gzip', compression_opts=6)
        h5['bins'].create_dataset('weight', data=bias, **h5opts)
    print("Time for normalization: {:.2f} seconds".format(time.time() - t15))


    # The bounds of the contact coordinates
    min_x = 0
    min_y = 0
    max_x = c.shape[0]
    max_y = c.shape[1]
    print("max_x:", max_x)
    print("max_y:", max_y)


    # Range queries
    square_size = args.square
    t2 = time.time()

    for i in range(args.iterations):
        point1 = random.randint(min_x, max_x - square_size)
        point2 = random.randint(min_y, max_y - square_size)
        mat = matrix[point1 : point1+square_size, point2 : point2+square_size]
        selected_points = list(zip(mat.row, mat.col, mat.data))

    t25 = time.time()
    print("Time performing range queries (256x256): {:.2f} seconds (per query): {:.2f} seconds".format(t25 - t2, (t25 - t2) / args.iterations))

    for i in range(args.iterations):
        point1 = random.randint(min_x, max_x - square_size * 8)
        point2 = random.randint(min_y, max_y - square_size * 8)
        mat = matrix[point1 : point1+square_size*8, point2 : point2+square_size*8]
        selected_points = list(zip(mat.row, mat.col, mat.data))

    t3 = time.time()
    print("Time performing range queries (2048 x 2048): {:.2f} seconds (per query): {:.2f} seconds".format(t3 - t25, (t3 - t25) / args.iterations))

    for i in range(args.iterations):
        point1 = random.randint(min_x, max_x - square_size)
        mat = matrix[point1, :]
        selected_points = list(zip(mat.row, mat.col, mat.data))

    t4 = time.time()
    print("Time slicing across first dimension: {:.2f} seconds (per query): {:.2f} seconds".format(t4 - t3, (t4 - t3) / args.iterations))

    for i in range(args.iterations):
        point2 = random.randint(min_y, max_y - square_size)
        mat = matrix[:, point2]
        selected_points = list(zip(mat.row, mat.col, mat.data))

    t5 = time.time()
    print("Time slicing across second dimension: {:.2f} seconds (per query): {:.2f} seconds".format(t5 - t4, (t5 - t4) / args.iterations))

    selected_points = []
    for i in range(args.iterations):
        for pix in c.pixels().iterchunks(size=100000):
            diag = pix[pix.bin1_id == pix.bin2_id]
            selected_points.extend( list(zip(diag['bin1_id'], diag['bin2_id'], diag['count'])) )

    t6 = time.time()
    print("Time slicing across the diagonal: {:.2f} seconds (per query): {:.2f} seconds".format(t6 - t5, (t6 - t5) / args.iterations))


    # Dump
    print("Size of index: {} bytes".format(op.getsize(outfilepath)))
    with open('/tmp/tmp.tsv', 'wt') as f:
        for pix in c.pixels().iterchunks(size=100000):
            pix.to_csv(f, sep='\t', index=False, header=False)
    print("Time outputting the index: {:.2f}".format(time.time() - t6))
    print("Size of output: {} bytes".format(op.getsize('/tmp/tmp.tsv')))


if __name__ == '__main__':
    main()
