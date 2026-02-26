"""spec2nii module containing utility functions 
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2026 University of Oxford
"""
from chardet.universaldetector import UniversalDetector
from pathlib import Path


def detect_file_encoding(filepath: Path | str) -> str:
    """Use chardet to estimate the encoding of the file

    :param filepath: Path to file
    :type filepath: Path | str
    :return: Most likely encoding string
    :rtype: str
    """

    detector = UniversalDetector()
    with open(filepath, 'rb') as f:
        for line in f:
            detector.feed(line)
            if detector.done:
                break
    detector.close()

    return detector.result['encoding']
