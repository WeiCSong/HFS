"""
Description:
    This script is run after `fasta_cli.py`. It computes the sequence-class
    level variant effect scores based on Sei chromatin profile predictions

Usage:
    sequence_class.py <ref>
    sequence_class.py -h | --help

Options:
    <ref>    Ref HDF5 filepath.

"""
import os
from os.path import exists
from docopt import docopt
import h5py
import numpy as np
import pandas as pd



def get_targets(filename):
    targets = []
    with open(filename, 'r') as file_handle:
        for line in file_handle:
            targets.append(line.strip())
    return targets


def get_data(filename):
    fh = h5py.File(filename, 'r')
    data = fh["data"][()]
    fh.close()
    return data



if __name__ == "__main__":
    arguments = docopt(
        __doc__,
        version='1.0.0')
    ref_fp = arguments['<ref>']

    results_dir, filename = os.path.split(ref_fp)
    filename_prefix = '.'.join(filename.split('.')[:-1])
    print(filename_prefix)

    sei_dir = "/dssg/home/acct-bmelgn/bmelgn-3/FLAT/seimodel/"
    chromatin_profiles = get_targets(os.path.join(sei_dir, "target.names"))
    seqclass_names = get_targets(os.path.join(sei_dir, "seqclass.names"))


    chromatin_profile_ref = get_data(ref_fp)

    clustervfeat = np.load(os.path.join(sei_dir, 'projvec_targets.npy'))
    histone_inds = np.load(os.path.join(sei_dir, 'histone_inds.npy'))
    
        
    refproj = (np.dot(chromatin_profile_ref,clustervfeat.T) /
               np.linalg.norm(clustervfeat,axis=1))
    dat=pd.DataFrame(refproj)
    dat=dat.iloc[:,0:39]
    dat.to_csv(ref_fp.replace("h5", "txt"),sep="\t")
