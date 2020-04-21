"""Main program for FMI-PPN.

FMI-PPN (Finnish Meteorological Institute Probabilistic Precipitation
Nowcaster) is a script for weather radar-based nowcasting. It is
heavily based on pySTEPS initiative.

For more information about pySTEPS, see https://pysteps.github.io .

Author: Petteri Karsisto
Year: 2019
"""
import argparse

import ppn


def get_input_arguments():
    """Wrapper for reading input arguments from command line (via argparse).

    Return dictionary object with parsed arguments.
    """
    # Add more options if necessary
    parser = argparse.ArgumentParser(description="Command line interface for FMI-PPN")
    parser.add_argument("-c", "--config", help="Select configuration settings")
    parser.add_argument("-t", "--timestamp", help="Nowcast initialization time",
                        metavar="YYYYMMDDHHMM")

    return vars(parser.parse_args())


def main():
    """Pass commandline arguments to ppn.run() method."""
    # Read command line arguments
    args = get_input_arguments()
    # Unpack dictionary into keyword arguments
    # Unused arguments should be ignored silently.
    ppn.run(**args)


if __name__ == "__main__":
    main()
