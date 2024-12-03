from pathlib import Path
import numpy as np
import copy
import warnings
from typing import Union, List, Dict, Sequence, Optional

from .base import AncestryObject


class SNPLevelAncestryObject(AncestryObject):
    """
    A class for SNP-level Local Ancestry Inference (LAI) data, including genotype and variant metadata.
    """
    def __init__(
        self,
        haplotypes: List[str],
        lai: np.ndarray,
        calldata_gt: Optional[np.ndarray] = None,
        samples: Optional[np.ndarray] = None,
        ancestry_map: Optional[Dict[str, str]] = None,
        variants_ref: Optional[np.ndarray] = None,
        variants_alt: Optional[np.ndarray] = None,
        variants_chrom: Optional[np.ndarray] = None,
        variants_filter_pass: Optional[np.ndarray] = None,
        variants_id: Optional[np.ndarray] = None,
        variants_pos: Optional[np.ndarray] = None,
        variants_qual: Optional[np.ndarray] = None
    ) -> None:
        """
        Args:
            haplotypes (list of str of length n_haplotypes): 
                A list of unique haplotype identifiers.
            lai (array of shape (n_snps, n_haplotypes)): 
                A 2D array containing local ancestry inference values, where each row represents a SNP,
                and each column corresponds to a haplotype phase for each sample.
            calldata_gt (array, optional): 
                An array containing genotype data for each sample. This array is 3D with shape 
                `(n_snps, n_samples, 2)`.
            samples (array of shape (n_sampels,), optional): 
                An array containing unique sample identifiers.
            ancestry_map (dict of str to str, optional): 
                A dictionary mapping ancestry codes to region names.
            variants_ref (array of shape (n_snps,), optional): 
                An array containing the reference allele for each SNP.
            variants_alt (array of shape (n_snps,), optional): 
                An array containing the alternate allele for each SNP.
            variants_chrom (array of shape (n_snps,), optional): 
                An array containing the chromosome for each SNP.
            variants_filter_pass (array of shape (n_snps,), optional): 
                An array indicating whether each SNP passed control checks.
            variants_id (array of shape (n_snps,), optional): 
                An array containing unique identifiers (IDs) for each SNP.
            variants_pos (array of shape (n_snps,), optional): 
                An array containing the chromosomal positions for each SNP.
            variants_qual (array of shape (n_snps,), optional): 
                An array containing the Phred-scaled quality score for each SNP.
        """
        if lai.ndim != 2:
            raise ValueError("`lai` must be a 2D array with shape (n_windows, n_haplotypes).")
        
        # Determine the number of unique ancestries and samples from the LAI array
        n_ancestries = len(np.unique(lai))
        n_haplotypes = lai.shape[1]
        n_samples = n_haplotypes // 2

        super(SNPLevelAncestryObject, self).__init__(n_samples, n_ancestries)

        self.__haplotypes = haplotypes
        self.__lai = lai
        self.__calldata_gt = calldata_gt
        self.__samples = samples
        self.__ancestry_map = ancestry_map
        self.__variants_ref = variants_ref
        self.__variants_alt = variants_alt
        self.__variants_chrom = variants_chrom
        self.__variants_filter_pass = variants_filter_pass
        self.__variants_id = variants_id
        self.__variants_pos = variants_pos
        self.__variants_qual = variants_qual

        # Perform sanity check to ensure all unique ancestries in LAI data are represented in the ancestry map
        self._sanity_check()

    def __getitem__(self, key):
        """
        To access an attribute of the class using the square bracket notation,
        similar to a dictionary.
        """
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(f'Invalid key: {key}')

    def __setitem__(self, key, value):
        """
        To set an attribute of the class using the square bracket notation,
        similar to a dictionary.
        """
        try:
            setattr(self, key, value)
        except AttributeError:
            raise KeyError(f'Invalid key: {key}')

    @property
    def haplotypes(self) -> List[str]:
        """
        Retrieve `haplotypes`.

        Returns:
            **list of length n_haplotypes:** A list of unique haplotype identifiers.
        """
        return self.__haplotypes

    @haplotypes.setter
    def haplotypes(self, x):
        """
        Update `haplotypes`.
        """
        self.__haplotypes = x

    @property
    def lai(self) -> np.ndarray:
        """
        Retrieve `lai`.

        Returns:
            **array of shape (n_windows, n_haplotypes):** 
                A 2D array containing local ancestry inference values, where each row represents a SNP,
                and each column corresponds to a haplotype phase for each sample.
        """
        return self.__lai

    @lai.setter
    def lai(self, x):
        """
        Update `lai`.
        """
        self.__lai = x

    @property
    def calldata_gt(self) -> np.ndarray:
        """
        Retrieve `calldata_gt`.

        Returns:
            **array:** 
                An array containing genotype data for each sample. This array is 3D with shape 
                `(n_snps, n_samples, 2)`.
        """
        return self.__calldata_gt

    @calldata_gt.setter
    def calldata_gt(self, x: np.ndarray):
        """
        Update `calldata_gt`.
        """
        self.__calldata_gt = x

        @property
    def samples(self) -> Optional[np.ndarray]:
        """
        Retrieve `samples`.

        Returns:
            **array of shape (n_sampels,):** 
                An array containing unique sample identifiers.
        """
        return self.__samples

    @samples.setter
    def samples(self, x: np.ndarray):
        """
        Update `samples`.
        """
        self.__samples = x

    @property
    def ancestry_map(self) -> Optional[Dict[str, str]]:
        """
        Retrieve `ancestry_map`.

        Returns:
            **dict of str to str:** A dictionary mapping ancestry codes to region names.
        """
        return self.__ancestry_map

    @ancestry_map.setter
    def ancestry_map(self, x):
        """
        Update `ancestry_map`.
        """
        self.__ancestry_map = x

    @property
    def variants_ref(self) -> Optional[np.ndarray]:
        """
        Retrieve `variants_ref`.

        Returns:
            **array of shape (n_snps,):** An array containing the reference allele for each SNP.
        """
        return self.__variants_ref

    @variants_ref.setter
    def variants_ref(self, x: np.ndarray):
        """
        Update `variants_ref`.
        """
        self.__variants_ref = x

    @property
    def variants_alt(self) -> Optional[np.ndarray]:
        """
        Retrieve `variants_alt`.

        Returns:
            **array of shape (n_snps,):** An array containing the alternate allele for each SNP.
        """
        return self.__variants_alt

    @variants_alt.setter
    def variants_alt(self, x: np.ndarray):
        """
        Update `variants_alt`.
        """
        self.__variants_alt = x

    @property
    def variants_chrom(self) -> Optional[np.ndarray]:
        """
        Retrieve `variants_chrom`.

        Returns:
            **array of shape (n_snps,):** An array containing the chromosome for each SNP.
        """
        return self.__variants_chrom

    @variants_chrom.setter
    def variants_chrom(self, x: np.ndarray):
        """
        Update `variants_chrom`.
        """
        self.__variants_chrom = x

    @property
    def variants_filter_pass(self) -> Optional[np.ndarray]:
        """
        Retrieve `variants_filter_pass`.

        Returns:
            **array of shape (n_snps,):** An array indicating whether each SNP passed control checks.
        """
        return self.__variants_filter_pass

    @variants_filter_pass.setter
    def variants_filter_pass(self, x: np.ndarray):
        """
        Update `variants_filter_pass`.
        """
        self.__variants_filter_pass = x

    @property
    def variants_id(self) -> Optional[np.ndarray]:
        """
        Retrieve `variants_id`.

        Returns:
            **array of shape (n_snps,):** An array containing unique identifiers (IDs) for each SNP.
        """
        return self.__variants_id

    @variants_id.setter
    def variants_id(self, x: np.ndarray):
        """
        Update `variants_id`.
        """
        self.__variants_id = x

    @property
    def variants_pos(self) -> Optional[np.ndarray]:
        """
        Retrieve `variants_pos`.

        Returns:
            **array of shape (n_snps,):** An array containing the chromosomal positions for each SNP.
        """
        return self.__variants_pos

    @variants_pos.setter
    def variants_pos(self, x: np.ndarray):
        """
        Update `variants_pos`.
        """
        self.__variants_pos = x

    @property
    def variants_qual(self) -> Optional[np.ndarray]:
        """
        Retrieve `variants_qual`.

        Returns:
            **array of shape (n_snps,):** An array containing the Phred-scaled quality score for each SNP.
        """
        return self.__variants_qual

    @variants_qual.setter
    def variants_qual(self, x: np.ndarray):
        """
        Update `variants_qual`.
        """
        self.__variants_qual = x

    @property
    def n_samples(self) -> int:
        """
        Retrieve `n_samples`.

        Returns:
            **int:** 
                The total number of samples.
        """
        if self.__samples is not None:
            return len(self.__samples)
        elif self.__haplotypes is not None:
            # Divide by 2 because each sample has two associated haplotypes
            return len(self.__haplotypes) // 2
        else:
            # Divide by 2 because columns represent haplotypes
            return self.__lai.shape[1] // 2

    @property
    def n_ancestries(self) -> int:
        """
        Retrieve `n_ancestries`.

        Returns:
            **int:** The total number of unique ancestries.
        """
        return len(np.unique(self.__lai))
    
    @property
    def n_haplotypes(self) -> int:
        """
        Retrieve `n_haplotypes`.

        Returns:
            **int:** The total number of haplotypes.
        """
        if self.__haplotypes is not None:
            return len(self.__haplotypes)
        else:
            return self.__lai.shape[1]

    @property
    def n_windows(self) -> int:
        """
        Retrieve `n_windows`.

        Returns:
            **int:** The total number of genomic windows.
        """
        return self.__lai.shape[0]

    def copy(self) -> 'SNPLevelAncestryObject':
        """
        Create and return a copy of `self`.

        Returns:
            **SNPLevelAncestryObject:** 
                A new instance of the current object.
        """
        return copy.copy(self)

    def keys(self) -> List[str]:
        """
        Retrieve a list of public attribute names for `self`.

        Returns:
            **list of str:** 
                A list of attribute names, with internal name-mangling removed, 
                for easier reference to public attributes in the instance.
        """
        return [attr.replace('_SNPLevelAncestryObject__', '').replace('_AncestryObject__', '') for attr in vars(self)]

    def save_pickle(self, file: Union[str, Path]) -> None:
        """
        Save `self` in serialized form to a `.pkl` file.
        If the file already exists, it will be overwritten.

        Args:
            file (str or pathlib.Path): 
                Path to the file where the data will be saved. It should end with `.pkl`. 
                If the provided path does not have this extension, it will be appended.
        """
        import pickle
        with open(file, 'wb') as file:
            pickle.dump(self, file)
