

"""
Description:
    This CLI script takes in a FASTA file of 4096bp sequences
    (note each description line must have unique identifiers)
    and writes the predictions to the specified output directory.

Output:
    Sei chromatin profile model predictions

Usage:
    fasta_cli.py <fasta> <output-dir>
    fasta_cli.py -h | --help

Options:
    -h --help               Show this screen.

    <fasta>                 Input FASTA
    <output-dir>            Output directory
"""

import os

from docopt import docopt

from selene_sdk.sequences import Genome
from selene_sdk.utils import load_path
from selene_sdk.utils import parse_configs_and_run
from selene_sdk import __version__


if __name__ == "__main__":
    arguments = docopt(
        __doc__,
        version=__version__)

    def run_config(config_yml):
        configs = load_path(config_yml, instantiate=False)
        configs["prediction"]["input_path"] = arguments["<fasta>"]
        configs["prediction"]["output_dir"] = arguments["<output-dir>"]
        parse_configs_and_run(configs)
    print(arguments)
    run_config("/dssg/home/acct-bmelgn/bmelgn-3/FLAT/code/seqs_fasta.yml")




