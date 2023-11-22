"""
A CLI for searching the CADC CFHT MegaPrime archive for observations of a moving object in a time interval.
"""
import argparse
import logging
import os.path
import warnings

import vos
from astropy.time import Time
from .core import search


def main():
    """
    The CLI entry point.
    """
    parser = argparse.ArgumentParser(description='Search the CADC CFHT MegaPrime archive for observations of '
                                                 'a moving object in a time interval.')
    parser.add_argument('ast_file', type=str, help='The path to the ast file containing '
                                                   'the observations of the ')
    parser.add_argument('start_date', type=str, help='The start date of the time interval to search')
    parser.add_argument('end_date', type=str, help='The end date of the time interval to search')
    parser.add_argument('--output', type=str,
                        help="Name of output file defaults to stdout",
                        default=None)
    parser.add_argument('--format', type=str,
                        help=('The output format, an astropy write format string.  '
                              'See https://docs.astropy.org/en/stable/io/unified.html#built-in-readers-writers'),
                        default='ascii.csv')
    parser.add_argument('--delimiter', type=str, help='The output column delimiter', default='\t')
    parser.add_argument('--dblist', help="Text file containing list of observation ids of interest.",
                        default=None)
    parser.add_argument('--ssois', help="Use the SSOIS service to find observations of the moving object.",
                        action='store_true', default=False)
    parser.add_argument('--log-level', type=str, help='The log level to use',
                        default='ERROR',
                        choices=['DEBUG',
                                 'INFO',
                                 'WARNING',
                                 'ERROR',
                                 'CRITICAL'])
    args = parser.parse_args()
    warnings.simplefilter('ignore', category=UserWarning)
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s %(levelname)s %(message)s')
    start_date = Time(args.start_date)
    end_date = Time(args.end_date)
    ast_file = args.ast_file
    observation_id_list = None
    if args.dblist is not None:
        observation_id_list = [observation_id.strip() for observation_id in open(args.dblist).readlines()]

    result_table = search(ast_file, start_date, end_date, observation_id_list, use_ssois=args.ssois)
    result_table.sort('MJD')
    result_table.write(args.output,
                       format='ascii',
                       delimiter='\t',
                       overwrite=True)
