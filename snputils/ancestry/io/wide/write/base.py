import abc
from pathlib import Path
from typing import Union

from snputils.ancestry.genobj.wide import GlobalAncestryObject

class WideBaseWriter(abc.ABC):
    """
    Abstract class for global ancestry writers.
    """ 
    def __init__(
        self, 
        wideobj: GlobalAncestryObject, 
        file_prefix: Union[str, Path]
    ) -> None:
        """
        Args:
            wideobj (GlobalAncestryObject): 
                A wide ancestry object instance.
            file_prefix (str or pathlib.Path): 
                The prefix for the output file names, including any parent directories. 
                This prefix is used to generate the output file names and should not include 
                the file extensions.
        """
        self.__wideobj = wideobj
        self.__file_prefix = Path(file_prefix)

    @abc.abstractmethod
    def write(self) -> None:
        """
        Abstract method to write wide ancestry files.

        Subclasses must implement this method to write the data from
        the Q and P files, as well as sample identifiers, SNPs, and ancestries if available.
        """
        pass

    @property
    def wideobj(self) -> GlobalAncestryObject:
        """
        Retrieve `wideobj`.

        Returns:
            GlobalAncestryObject: A wide ancestry object instance.
        """
        return self.__wideobj

    @property
    def file_prefix(self) -> Path:
        """
        Retrieve `file_prefix`.

        Returns:
            pathlib.Path: 
                The prefix for the output file names, including any parent directories. 
                This prefix is used to generate the output file names and should not include 
                the file extensions.
        """
        return self.__file_prefix
