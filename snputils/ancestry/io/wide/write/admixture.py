import logging
from pathlib import Path
from typing import Union, Optional
import numpy as np

from .base import WideBaseWriter
from snputils.ancestry.genobj.wide import GlobalAncestryObject

log = logging.getLogger(__name__)


class AdmixtureWriter(WideBaseWriter):
    """
    A class for writing data stored in a `GlobalAncestryObject` instance into the multiple ADMIXTURE files.
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
        super(AdmixtureWriter, self).__init__(wideobj, file_prefix)
        self.__Q_file = self.file_prefix.with_suffix(f".{self.wideobj.n_ancestries}.Q")
        self.__P_file = self.file_prefix.with_suffix(f".{self.wideobj.n_ancestries}.P")

        self.__sample_file = self.file_prefix.with_suffix(".sample_ids.txt") if self.wideobj.samples is not None else None
        self.__snp_file = self.file_prefix.with_suffix(".snp_ids.txt") if self.wideobj.snps is not None else None
        self.__ancestry_file = self.file_prefix.with_suffix(".map") if self.wideobj.ancestries is not None else None

    @property
    def wideobj(self) -> GlobalAncestryObject:
        """
        Retrieve `wideobj`.

        Returns:
            GlobalAncestryObject: A wide ancestry object instance.
        """
        return self._wideobj
    
    @property
    def file_prefix(self) -> Path:
        """
        Retrieve `file_prefix`.

        Returns:
            pathlib.Path: 
                The prefix for the output file names, including any parent directories. 
                For example, if `"parent1/fname"` is provided, the Q matrix will be saved as 
                `"parent1/fname.K.Q"` and the P matrix as `"parent1/fname.K.P"`, where `K` is 
                the number of ancestries.
        """
        return self._file_prefix

    @property
    def Q_file(self) -> str:
        """Retrieve `Q_file`.

        Returns:
            str: Path to store the file containing the Q matrix (per-sample ancestry proportions).
        """
        return self.__Q_file
    
    @property
    def P_file(self) -> str:
        """Retrieve `P_file`.

        Returns:
            str: Path to store the file containing the P/F matrix (per-ancestry SNP frequencies).
        """
        return self.__P_file
    
    @property
    def sample_file(self) -> Optional[str]:
        """Retrieve `sample_file`.

        Returns:
            str: Path to store the file containing sample identifiers. If None, sample identifiers are not saved.
        """
        return self.__sample_file
    
    @property
    def snp_file(self) -> Optional[str]:
        """Retrieve `snp_file`.

        Returns:
            str: Path to store the file containing SNP identifiers. If None, SNP identifiers are not saved.
        """
        return self.__snp_file
    
    @property
    def ancestry_file(self) -> Optional[str]:
        """Retrieve `ancestry_file`.

        Returns:
            str: Path to store the file containing ancestry labels for each sample. If None, ancestries are not saved.
        """
        return self.__ancestry_file

    def _write_Q(self):
        log.info(f"Writing Q matrix to '{self.Q_file}'...")
        np.savetxt(self.Q_file, self.wideobj.Q, delimiter=" ")
        log.info(f"Finished writing Q matrix to '{self.Q_file}'.")

    def _write_P(self):
        log.info(f"Writing P matrix to '{self.P_file}'...")
        np.savetxt(self.P_file, self.wideobj.P, delimiter=" ")
        log.info(f"Finished writing P matrix to '{self.P_file}'.")

    def _write_sample_ids(self):
        if self.wideobj.samples is not None:
            log.info(f"Writing sample IDs to '{self.sample_file}'...")
            np.savetxt(self.sample_file, self.wideobj.samples, fmt="%s")
            log.info(f"Finished writing sample IDs to '{self.sample_file}'.")

    def _write_snps(self):
        if self.wideobj.snps is not None:
            log.info(f"Writing SNP IDs to '{self.snp_file}'...")
            np.savetxt(self.snp_file, self.wideobj.snps, fmt="%s")
            log.info(f"Finished writing SNP IDs to '{self.snp_file}'.")

    def _write_ancestries(self):
        if self.wideobj.ancestries is not None:
            log.info(f"Writing ancestry information to '{self.ancestry_file}'...")
            np.savetxt(self.ancestry_file, self.wideobj.ancestries, fmt="%s")
            log.info(f"Finished writing ancestry information to '{self.ancestry_file}'.")

    def write(self) -> None:
        """
        Write the data contained in the `wideobj` instance into the multiple ADMIXTURE files.
        The output filenames are based on the specified `file_prefix` and include the following:

            - Q matrix file: `<file_prefix>.K.Q`.
            - P matrix file: `<file_prefix>.K.P`.
            - Sample IDs file: `<file_prefix>.sample_ids.txt` (if sample IDs are available).
            - SNP IDs file: `<file_prefix>.snp_ids.txt` (if SNP IDs are available).
            - Ancestry file: `<file_prefix>.map` (if ancestries information is available).

        where `K` is the total number of ancestries.
        """
        log.info(f"Preparing to write ADMIXTURE files with prefix '{self.file_prefix}'...")
        
        self.file_prefix.parent.mkdir(parents=True, exist_ok=True)
        
        self._write_Q()
        self._write_P()
        self._write_sample_ids()
        self._write_snps()
        self._write_ancestries()

        log.info(f"Finished writing all ADMIXTURE files with prefix '{self.file_prefix}'.")

WideBaseWriter.register(AdmixtureWriter)
