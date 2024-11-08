import abc
from pathlib import Path
from typing import Union

from snputils.ancestry.genobj.local import LocalAncestryObject


class LAIBaseReader(abc.ABC):
    """
    Abstract class for local ancestry readers.
    """
    def __init__(self, file: Union[str, Path]) -> None:
        """
        Args:
            file (str or pathlib.Path): 
                Path to the `.msp` file containing LAI info.
        """
        self.__file = Path(file)

    @property
    def file(self) -> Path:
        """
        Retrieve `file`.

        Returns:
            pathlib.Path: Path to the file containing LAI info.
        """
        return self.__file

    @abc.abstractmethod
    def read(self) -> 'LocalAncestryObject':
        """
        Abstract method to read data from the provided `file`.

        Subclasses must implement this method to read and parse the data. 
        The implementation should construct an instance of `snputils.ancestry.genobj.local.LocalAncestryObject` 
        based on the read data.
        """
        pass
