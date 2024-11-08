import abc
from pathlib import Path
from typing import Union

from snputils.ancestry.genobj.local import LocalAncestryObject


class LAIBaseWriter(abc.ABC):
    """
    Abstract class for local ancestry writers.
    """
    def __init__(self, laiobj: LocalAncestryObject, file=Union[str, Path]) -> None:
        """
        Args:
            laiobj (LocalAncestryObject):
                A local ancestry object instance.
            file (str or pathlib.Path): 
                Path to the output `.msp` file containing LAI info.
        """
        self.__laiobj = laiobj
        self.__file = Path(file)

    @property
    def laiobj(self) -> LocalAncestryObject:
        """
        Retrieve `laiobj`. 

        Returns:
            laiobj (LocalAncestryObject):
                A local ancestry object instance.
        """
        return self.__laiobj

    @property
    def file(self) -> str:
        """
        Retrieve `file`.

        Returns:
            pathlib.Path: Path to the output `.msp` file containing LAI info.
        """
        return self.__file

    @abc.abstractmethod
    def write(self) -> None:
        """
        Abstract method to write the data from `laiobj` to the specified output `file`.

        Subclasses must implement this method to perform the actual writing of data.
        """
        pass
