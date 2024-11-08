import pathlib
import numpy as np
import copy
import warnings
from typing import Union, List, Dict, Sequence, Optional

from .base import AncestryObject


class LocalAncestryObject(AncestryObject):
    """
    A class for Local Ancestry Inference (LAI) data.
    """
    def __init__(
        self,
        haplotypes: List,
        lai: np.ndarray,
        samples: Optional[List] = None,
        ancestry_map: Optional[Dict] = None,
        window_sizes: Optional[np.ndarray] = None,
        centimorgan_pos: Optional[np.ndarray] = None,
        chromosomes: Optional[np.ndarray] = None,
        physical_pos: Optional[np.ndarray] = None
    ) -> None:
        """
        Args:
            haplotypes (list): 
                A list of unique haplotype identifiers with a length of `n_haplotypes`.
            lai (array of shape (n_windows, n_haplotypes)): 
                A 2D array containing local ancestry inference values, where each row represents a 
                genomic window, and each column corresponds to a haplotype phase for each sample.
            samples (list, optional): 
                A list of unique sample identifiers with a length of `n_samples`.
            ancestry_map (dict, optional): 
                A dictionary mapping ancestry codes to region names.
            window_sizes (array of shape (n_windows,), optional): 
                An array specifying the number of SNPs in each genomic window.
            centimorgan_pos (array of shape (n_windows, 2), optional): 
                A 2D array containing the start and end centimorgan positions for each window.
            chromosomes (array of shape (n_windows,), optional): 
                An array with chromosome numbers corresponding to each genomic window.
            physical_pos (array of shape (n_windows, 2), optional): 
                A 2D array containing the start and end physical positions for each window.
        """
        if lai.ndim != 2:
            raise ValueError("`lai` must be a 2D array with shape (n_windows, n_haplotypes).")
        
        # Determine the number of unique ancestries and samples from the LAI array
        n_ancestries = len(np.unique(lai))
        n_haplotypes = lai.shape[1]
        n_samples = n_haplotypes // 2

        super(LocalAncestryObject, self).__init__(n_samples, n_ancestries)

        self.__haplotypes = haplotypes
        self.__lai = lai
        self.__window_sizes = window_sizes
        self.__centimorgan_pos = centimorgan_pos
        self.__samples = samples
        self.__chromosomes = chromosomes
        self.__physical_pos = physical_pos
        self.__ancestry_map = ancestry_map

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
    def haplotypes(self) -> List:
        """
        Retrieve `haplotypes`.

        Returns:
            List: A list of unique haplotype identifiers with a length of `n_samples*2`.
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
            numpy.ndarray: A 2D array containing local ancestry inference values, where each row represents a 
                           genomic window, and each column corresponds to a haplotype phase for each sample.
        """
        return self.__lai

    @lai.setter
    def lai(self, x):
        """
        Update `lai`.
        """
        self.__lai = x

    @property
    def samples(self) -> Optional[List]:
        """
        Retrieve `samples`.

        Returns:
            List: A list of unique sample identifiers with a length of `n_samples`.
        """
        return self.__samples
    
    @samples.setter
    def samples(self, x):
        """
        Update `samples`.
        """
        self.__samples = x

    @property
    def ancestry_map(self) -> Optional[Dict]:
        """
        Retrieve `ancestry_map`.

        Returns:
            Dict: A dictionary mapping ancestry codes to region names.
        """
        return self.__ancestry_map

    @ancestry_map.setter
    def ancestry_map(self, x):
        """
        Update `ancestry_map`.
        """
        self.__ancestry_map = x

    @property
    def window_sizes(self) -> Optional[np.ndarray]:
        """
        Retrieve `window_sizes`.

        Returns:
            numpy.ndarray: An array specifying the number of SNPs in each genomic window.
        """
        return self.__window_sizes
        
    @window_sizes.setter
    def window_sizes(self, x):
        """
        Update `window_sizes`.
        """
        self.__window_sizes = x

    @property
    def centimorgan_pos(self) -> Optional[np.ndarray]:
        """
        Retrieve `centimorgan_pos`.

        Returns:
            numpy.ndarray: A 2D array containing the start and end centimorgan positions for each window.
        """
        return self.__centimorgan_pos

    @centimorgan_pos.setter
    def centimorgan_pos(self, x):
        """
        Update `centimorgan_pos`.
        """
        self.__centimorgan_pos = x

    @property
    def chromosomes(self) -> Optional[np.ndarray]:
        """
        Retrieve `chromosomes`.

        Returns:
            numpy.ndarray: An array with chromosome numbers corresponding to each genomic window.
        """
        return self.__chromosomes
        
    @chromosomes.setter
    def chromosomes(self, x):
        """
        Update `chromosomes`.
        """
        self.__chromosomes = x

    @property
    def physical_pos(self) -> Optional[np.ndarray]:
        """
        Retrieve `physical_pos`.

        Returns:
            numpy.ndarray: A 2D array containing the start and end physical positions for each window.
        """
        return self.__physical_pos

    @physical_pos.setter
    def physical_pos(self, x):
        """
        Update `physical_pos`.
        """
        self.__physical_pos = x

    @property
    def n_samples(self) -> int:
        """
        Retrieve `n_samples`.

        Returns:
            int: The total number of samples. If `samples` is available, returns its length; 
                 otherwise, calculates based on the number of `haplotypes` or `lai` array dimensions.
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
            int: The total number of unique ancestries.
        """
        return len(np.unique(self.__lai))
    
    @property
    def n_haplotypes(self) -> int:
        """
        Retrieve `n_haplotypes`.

        Returns:
            int: The total number of haplotypes.
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
            int: The total number of genomic windows.
        """
        return self.__lai.shape[0]

    def copy(self) -> 'LocalAncestryObject':
        """
        Create and return a copy of the current `LocalAncestryObject` instance.

        Returns:
            LocalAncestryObject: 
                A new instance of the current object.
        """
        return copy.copy(self)

    def keys(self) -> List:
        """
        Retrieve a list of public attribute names for this `LocalAncestryObject` instance.

        Returns:
            List: A list of attribute names, with internal name-mangling removed, 
                  for easier reference to public attributes in the instance.
        """
        return [attr.replace('_LocalAncestryObject__', '').replace('_AncestryObject__', '') for attr in vars(self)]

    def filter_windows(
            self,
            indexes: Union[int, Sequence[int], np.ndarray],
            include: bool = True,
            inplace: bool = False
        ) -> Optional['LocalAncestryObject']:
        """
        Filter genomic windows in the `LocalAncestryObject` based on their indexes.

        This method allows inclusion or exclusion of specific genomic windows from the
        `LocalAncestryObject` by specifying their indexes. Negative indexes are supported
        and follow NumPy's indexing conventions. It updates the `lai`, `chromosomes`, 
        `centimorgan_pos`, and `physical_pos` attributes accordingly.

        Args:
            indexes (int or array-like of int): 
                Indexes of the windows to include or exclude. Can be a single integer or a
                sequence of integers. Negative indexes are supported.
            include (bool, default=True): 
                If True, includes only the specified windows. If False, excludes the specified
                windows. Default is True.
            inplace (bool, default=False): 
                If True, modifies `self` in place. If False, returns a new `LocalAncestryObject` with 
                the windows filtered. Default is False.

        Returns:
            Optional[LocalAncestryObject]: Returns a new `LocalAncestryObject` with the specified 
            windows filtered if `inplace=False`. If `inplace=True`, modifies the object in place and 
            returns None.
        """
        # Convert indexes to a NumPy array
        indexes = np.atleast_1d(indexes)

        # Get total number of windows
        n_windows = self.n_windows

        # Validate indexes, allowing negative indexes
        if np.any((indexes < -n_windows) | (indexes >= n_windows)):
            raise IndexError("One or more indexes are out of bounds.")

        # Create boolean mask
        mask = np.zeros(n_windows, dtype=bool)
        mask[indexes] = True

        # Invert mask if `include=False`
        if not include:
            mask = ~mask
        
        # Filter `lai`
        filtered_lai = self['lai'][mask, :] 
        
        # Filter `chromosomes`, `centimorgan_pos`, and `physical_pos`, checking if they are None before filtering
        filtered_chromosomes = self['chromosomes'][mask] if self['chromosomes'] is not None else None
        filtered_centimorgan_pos = self['centimorgan_pos'][mask, :] if self['centimorgan_pos'] is not None else None
        filtered_physical_pos = self['physical_pos'][mask, :] if self['physical_pos'] is not None else None

        # Modify the original object if `inplace=True`, otherwise create and return a copy
        if inplace:
            self['lai'] = filtered_lai
            if filtered_chromosomes is not None:
                self['chromosomes'] = filtered_chromosomes
            if filtered_centimorgan_pos is not None:
                self['centimorgan_pos'] = filtered_centimorgan_pos
            if filtered_physical_pos is not None:
                self['physical_pos'] = filtered_physical_pos
            return None
        else:
            laiobj = self.copy()
            laiobj['lai'] = filtered_lai
            if filtered_chromosomes is not None:
                laiobj['chromosomes'] = filtered_chromosomes
            if filtered_centimorgan_pos is not None:
                laiobj['centimorgan_pos'] = filtered_centimorgan_pos
            if filtered_physical_pos is not None:
                laiobj['physical_pos'] = filtered_physical_pos
            return laiobj

    def filter_samples(
        self,
        samples: Union[str, Sequence[str], np.ndarray, None] = None,
        indexes: Union[int, Sequence[int], np.ndarray, None] = None,
        include: bool = True,
        inplace: bool = False
    ) -> Optional['LocalAncestryObject']:
        """
        Filter samples in the `LocalAncestryObject` based on sample names or indexes.

        This method allows inclusion or exclusion of specific samples by their names,
        indexes, or both. When both samples and indexes are provided, the union of
        the specified samples is used. Negative indexes are supported and follow NumPy's indexing 
        conventions. It updates the `lai`, `samples`, and `haplotypes` attributes accordingly.

        Args:
            samples (str or array_like of str, optional): 
                 Names of the samples to include or exclude. Can be a single sample name or a
                 sequence of sample names. Default is None.
            indexes (int or array_like of int, optional):
                Indexes of the samples to include or exclude. Can be a single index or a sequence
                of indexes. Negative indexes are supported. Default is None.
            include (bool, default=True): 
                If True, includes only the specified samples. If False, excludes the specified
                samples. Default is True.
            inplace (bool, default=False): 
                If True, modifies `self` in place. If False, returns a new `LocalAncestryObject` with the 
                samples filtered. Default is False.

        Returns:
            Optional[LocalAncestryObject]: A new LocalAncestryObject with the specified samples 
            filtered if `inplace=False`. If inplace=True, modifies `self` in place and returns None.
        """
        if samples is None and indexes is None:
            raise UserWarning("At least one of 'samples' or 'indexes' must be provided.")

        n_haplotypes = self.n_haplotypes
        n_samples = self.n_samples

        # Create mask based on sample names
        if samples is not None:
            samples = np.atleast_1d(samples)
            # Extract sample names from haplotype identifiers
            haplotype_ids = np.array(self['haplotypes'])
            sample_names = np.array([hap.split('.')[0] for hap in haplotype_ids])
            # Create mask for haplotypes belonging to specified samples
            mask_samples = np.isin(sample_names, samples)
        else:
            mask_samples = np.zeros(n_haplotypes, dtype=bool)

        # Create mask based on sample indexes
        if indexes is not None:
            indexes = np.atleast_1d(indexes)

            # Validate indexes, allowing negative indexes
            out_of_bounds_indexes = indexes[(indexes < -n_samples) | (indexes >= n_samples)]
            if out_of_bounds_indexes.size > 0:
                raise ValueError(f"One or more sample indexes are out of bounds.")

            # Adjust negative indexes
            indexes = np.mod(indexes, n_samples)
            
            # Get haplotype indexes for the specified sample indexes
            haplotype_indexes = np.concatenate([2*indexes, 2*indexes+1])
            # Create mask for haplotypes
            mask_indexes = np.zeros(n_haplotypes, dtype=bool)
            mask_indexes[haplotype_indexes] = True
        else:
            mask_indexes = np.zeros(n_haplotypes, dtype=bool)

        # Combine masks using logical OR (union of samples)
        mask_combined = mask_samples | mask_indexes

        if not include:
            mask_combined = ~mask_combined

        # Filter `lai`
        filtered_lai = self['lai'][:, mask_combined]

        # Filter `haplotypes`
        filtered_haplotypes = np.array(self['haplotypes'])[mask_combined].tolist()

        # Filter `samples`, checking if they are None before filtering
        sample_mask = mask_combined.reshape(-1, 2).any(axis=1)
        filtered_samples = np.array(self['samples'])[sample_mask].tolist() if self['samples'] is not None else None

        if inplace:
            self['haplotypes'] = filtered_haplotypes
            self['samples'] = filtered_samples
            self['lai'] = filtered_lai
            return None
        else:
            laiobj = self.copy()
            laiobj['haplotypes'] = filtered_haplotypes
            laiobj['samples'] = filtered_samples
            laiobj['lai'] = filtered_lai
            return laiobj

    def _sanity_check(self) -> None:
        """
        Perform sanity checks on the parsed data to ensure data integrity.

        This method checks that all unique ancestries in LAI are represented 
        in the ancestry map.

        Args:
            lai (np.ndarray): The LAI data array.
            ancestry_map (dict, optional): A dictionary mapping ancestry codes to region names, if available.
        """
        # Get unique ancestries from LAI data
        unique_ancestries = np.unique(self.lai)

        if self.ancestry_map is not None:
            # Check if all unique ancestries in the LAI are present in the ancestry map
            for ancestry in unique_ancestries:
                if str(ancestry) not in self.ancestry_map:
                    warnings.warn(
                        f"Ancestry '{ancestry}' found in LAI data is not represented in the ancestry map."
                    )

    def save(self, file: Union[str, pathlib.Path]) -> None:
        """
        Save the data stored in the `LocalAncestryObject` instance to a `.msp` file.

        Args:
            file (str or pathlib.Path): 
                The path to the file where the data will be saved. It should end with `.msp`. 
                If the provided path does not have this extension, it will be appended.
        """
        from snputils.ancestry.io.local.write import MSPWriter

        MSPWriter(self, file).write()
