from __future__ import annotations

import pathlib
from typing import Union


class LAIReader:
    def __new__(cls,
                filename: Union[str, pathlib.Path]) -> LAIReader:
        """
        Automatically detect the local ancestry data file format from the file extension, and return its corresponding reader.

        Args:
            filename: Filename of the file to read.

        Raises:
            ValueError: If the filename does not have an extension or the extension is not supported.
        """
        filename = pathlib.Path(filename)
        suffixes = filename.suffixes
        if not suffixes:
            raise ValueError("The filename should have an extension when using LAIReader.")

        extension = suffixes[-1].lower()

        if extension == ".msp":
            from snputils.ancestry.io.local.read.msp import MSPReader

            return MSPReader(filename)
        else:
            raise ValueError(f"File format not supported: {filename}")
