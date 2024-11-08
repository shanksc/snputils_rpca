import abc
from pathlib import Path
from typing import Union, Optional
import logging
import numpy as np

log = logging.getLogger(__name__)


class WideBaseReader(abc.ABC):
    """
    Abstract class for global ancestry readers.
    """
    def __init__(
        self,
        Q_file: Union[str, Path],
        P_file: Union[str, Path],
        sample_file: Optional[Union[str, Path]] = None,
        snp_file: Optional[Union[str, Path]] = None,
        ancestry_file: Optional[Union[str, Path]] = None,
    ) -> None:
        """
        Args:
            Q_file (str or pathlib.Path):
                Path to the file containing the Q matrix (per-sample ancestry proportions).
            P_file (str or pathlib.Path):
                Path to the file containing the P/F matrix (per-ancestry SNP frequencies).
            sample_file (str or pathlib.Path, optional):
                Path to the file containing sample identifiers. If None, sample identifiers are not loaded.
            snp_file (str or pathlib.Path, optional):
                Path to the file containing SNP identifiers. If None, SNP identifiers are not loaded.
            ancestry_file (str or pathlib.Path, optional):
                Path to the file containing ancestry labels for each sample. If None, ancestries are not loaded.
        """
        self.__Q_file = Path(Q_file)
        self.__P_file = Path(P_file)
        self.__sample_file = Path(sample_file) if sample_file is not None else None
        self.__snp_file = Path(snp_file) if snp_file is not None else None
        self.__ancestry_file = Path(ancestry_file) if ancestry_file is not None else None

    @property
    def Q_file(self) -> Path:
        """
        Retrieve `Q_file`.

        Returns:
            pathlib.Path: Path to the file containing the Q matrix (per-sample ancestry proportions).
        """
        return self.__Q_file

    @property
    def P_file(self) -> Path:
        """
        Retrieve `P_file`.

        Returns:
            pathlib.Path: Path to the file containing the P/F matrix (per-ancestry SNP frequencies).
        """
        return self.__P_file

    @property
    def sample_file(self) -> Optional[Path]:
        """
        Retrieve `sample_file`.

        Returns:
            pathlib.Path: Path to the file containing sample identifiers. If None, sample identifiers are not loaded.
        """
        return self.__sample_file
    
    @property
    def snp_file(self) -> Optional[Path]:
        """
        Retrieve `snp_file`.

        Returns:
            pathlib.Path: Path to single-column text file storing SNP ID in order.
        """
        return self.__snp_file

    @property
    def ancestry_file(self) -> Optional[Path]:
        """
        Retrieve `ancestry_file`.

        Returns:
            pathlib.Path: Path to the file containing ancestry labels for each sample. If None, ancestries are not loaded.
        """
        return self.__ancestry_file

    def _read_sample_ids(self) -> Optional[np.ndarray]:
        """
        Read sample identifiers from `sample_file` file.
        """
        if self.sample_file is None:
            return None
        
        log.info(f"Reading sample identifiers from '{self.sample_file}'...")
        if self.sample_file.suffix == ".txt":
            return np.genfromtxt(self.sample_file, dtype=str)
        elif self.sample_file.suffix == ".fam":
            return np.genfromtxt(self.sample_file, dtype=str, usecols=1)
        else:
            raise ValueError("Invalid file format for sample identifiers. Should be a single-column txt file or a .fam file.")
    
    def _read_snps(self) -> Optional[np.ndarray]:
        """
        Read SNP identifiers from `snp_file` file.
        """
        if self.snp_file is None:
            return None
        
        log.info(f"Reading SNP identifiers from '{self.snp_file}'...")
        if self.snp_file.suffix == ".txt":
            return np.genfromtxt(self.snp_file, dtype=str)
        elif self.snp_file.suffix == ".bim":
            return np.genfromtxt(self.snp_file, dtype=str, usecols=1)
        else:
            raise ValueError("Invalid file format for SNP identifiers. Should be a single-column txt file or a .bim file.")

    def _read_ancestries(self) -> Optional[np.ndarray]:
        """
        Read ancestries for each sample from `ancestry_file`.
        """
        if self.ancestry_file is None:
            return None
        
        log.info(f"Reading ancestries for each sample from '{self.ancestry_file}'...")
        if self.ancestry_file.suffix == ".map":
            return np.genfromtxt(self.ancestry_file, dtype=str)
        else:
            raise ValueError("Invalid file format for ancestry labels. Should be a single-column txt file.")

    @abc.abstractmethod
    def read(self) -> None:
        """
        Abstract method to read data from the provided files.

        Subclasses must implement this method to read and parse the data from
        the Q and P files, as well as sample identifiers, SNPs, and ancestries if provided.
        The implementation should construct an instance of `snputils.ancestry.genobj.GlobalAncestryObject` 
        based on the read data.
        """
        pass
