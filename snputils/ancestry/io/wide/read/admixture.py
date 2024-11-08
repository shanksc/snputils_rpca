import logging
import numpy as np
from pathlib import Path
from typing import Union, Optional

log = logging.getLogger(__name__)

from .base import WideBaseReader
from snputils.ancestry.genobj.wide import GlobalAncestryObject


class AdmixtureReader(WideBaseReader):
    """
    Reader class for ADMIXTURE output files.
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

    def read(self) -> 'GlobalAncestryObject':
        """
        Read ADMIXTURE output files and construct a `GlobalAncestryObject` instance.

        This method processes the provided Q and P matrix files to extract per-sample ancestry proportions
        and per-ancestry SNP frequencies. Optional files for sample identifiers, SNP identifiers, and 
        ancestry labels are also read, if provided.

        Returns:
            GlobalAncestryObject: 
                A global ancestry object instance.
        """
        log.info(f"Reading Q matrix from '{self.Q_file}'...")
        Q_mat = np.genfromtxt(self.Q_file, delimiter=' ')
        log.info(f"Reading P matrix from '{self.P_file}'...")
        P_mat = np.genfromtxt(self.P_file, delimiter=' ')

        samples = self._read_sample_ids()
        snps = self._read_snps()
        ancestries = self._read_ancestries()

        return GlobalAncestryObject(
            Q_mat,
            P_mat,
            samples=samples,
            snps=snps,
            ancestries=ancestries
        )

WideBaseReader.register(AdmixtureReader)
