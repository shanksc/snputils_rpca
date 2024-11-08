from pathlib import Path
from typing import Union

from snputils.ancestry.genobj.local import LocalAncestryObject


def read_lai(filename: Union[str, Path], **kwargs) -> LocalAncestryObject:
    """
    Automatically detect the file format and read it into a LocalAncestryObject.

    Args:
        filename: Filename of the file to read.
        **kwargs: Additional arguments passed to the reader method.

    Raises:
        ValueError: If the filename does not have an extension or the extension is not supported.
    """
    from snputils.ancestry.io.local.read.auto import LAIReader

    return LAIReader(filename).read(**kwargs)


def read_msp(file: Union[str, Path]) -> 'LocalAncestryObject':
    """
    Read data from an `.msp` file and construct a `LocalAncestryObject`.

    This function processes the input file to extract the necessary information including 
    the Q and P matrices, sample identifiers, SNP identifiers, and ancestry map.

    Args:
        file (str or pathlib.Path): 
            Path to the `.msp` file containing LAI info.

    Returns:
            LocalAncestryObject:
                A local ancestry object instance.
    """
    from snputils.ancestry.io.local.read.msp import MSPReader

    return MSPReader(file).read()
