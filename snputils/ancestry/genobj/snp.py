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
    def n_snps(self) -> int:
        """
        Retrieve `n_snps`.

        Returns:
            **int:** The total number of SNPs.
        """
        if self.lai is not None:
            return self.lai.shape[0]
        elif self.__calldata_gt is not None:
            return self.__calldata_gt.shape[0]
        elif self.variants_ref is not None:
            return len(self.variants_ref)
        elif self.variants_pos is not None:
            return len(self.variants_pos)
        else:
            raise ValueError("Unable to determine the total number of SNPs: no relevant data is available.")

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

    def filter_variants(
        self, 
        chrom: Optional[Union[str, Sequence[str], np.ndarray, None]] = None, 
        pos: Optional[Union[int, Sequence[int], np.ndarray, None]] = None, 
        indexes: Optional[Union[int, Sequence[int], np.ndarray, None]] = None, 
        include: bool = True, 
        inplace: bool = False
    ) -> Optional['SNPLevelAncestryObject']:
        """
        Filter variants based on specified chromosome names, variant positions, or variant indexes.

        This method updates the `lai`, `calldata_gt`, `variants_ref`, `variants_alt`, `variants_chrom`, 
        `variants_filter_pass`, `variants_id`, `variants_pos`, and `variants_qual` attributes to include or 
        exclude the specified variants. The filtering criteria can be based on chromosome names, variant 
        positions, or indexes. If multiple criteria are provided, their union is used for filtering. 
        The order of the variants is preserved.

        Negative indexes are supported and follow 
        [NumPy's indexing conventions](https://numpy.org/doc/stable/user/basics.indexing.html).

        Args:
            chrom (str or array_like of str, optional): 
                Chromosome(s) to filter variants by. Can be a single chromosome as a string or a sequence 
                of chromosomes. If both `chrom` and `pos` are provided, they must either have matching lengths 
                (pairing each chromosome with a position) or `chrom` should be a single value that applies to 
                all positions in `pos`. Default is None. 
            pos (int or array_like of int, optional): 
                Position(s) to filter variants by. Can be a single position as an integer or a sequence of positions. 
                If `chrom` is also provided, `pos` should either match `chrom` in length or `chrom` should be a 
                single value. Default is None.
            indexes (int or array_like of int, optional): 
                Index(es) of the variants to include or exclude. Can be a single index or a sequence
                of indexes. Negative indexes are supported. Default is None.
            include (bool, default=True): 
                If True, includes only the specified variants. If False, excludes the specified
                variants. Default is True.
            inplace (bool, default=False): 
                If True, modifies `self` in place. If False, returns a new `SNPLevelAncestryObject` with the variants 
                filtered. Default is False.

        Returns:
            **Optional[SNPLevelAncestryObject]:** 
                A new `SNPLevelAncestryObject` with the specified variants filtered if `inplace=False`. 
                If `inplace=True`, modifies `self` in place and returns None.
        """
        if chrom is None and pos is None and indexes is None:
            raise ValueError("At least one of 'chrom', 'pos', or 'indexes' must be provided.")

        n_snps = self.n_snps

        # Convert inputs to arrays for consistency
        chrom = np.atleast_1d(chrom) if chrom is not None else None
        pos = np.atleast_1d(pos) if pos is not None else None
        indexes = np.atleast_1d(indexes) if indexes is not None else None

        # Validate chrom and pos lengths if both are provided
        if chrom is not None and pos is not None:
            if len(chrom) != len(pos) and len(chrom) > 1:
                raise ValueError(
                    "When both 'chrom' and 'pos' are provided, they must either be of the same length "
                    "or 'chrom' must be a single value."
                )

        # Create a mask for chromosome and position filtering
        mask_combined = np.zeros(n_snps, dtype=bool)
        if chrom is not None and pos is not None:
            if len(chrom) == 1:
                # Apply single chromosome to all positions in `pos`
                mask_combined = (self['variants_chrom'] == chrom[0]) & np.isin(self['variants_pos'], pos)
            else:
                # Vectorized pair matching for chrom and pos
                query_pairs = np.array(
                    list(zip(chrom, pos)),
                    dtype=[
                        ('chrom', self['variants_chrom'].dtype),
                        ('pos', self['variants_pos'].dtype)
                    ]
                )
                data_pairs = np.array(
                    list(zip(self['variants_chrom'], self['variants_pos'])),
                    dtype=[
                        ('chrom', self['variants_chrom'].dtype),
                        ('pos', self['variants_pos'].dtype)
                    ]
                )
                mask_combined = np.isin(data_pairs, query_pairs)

        elif chrom is not None:
            # Only chromosome filtering
            mask_combined = np.isin(self['variants_chrom'], chrom)
        elif pos is not None:
            # Only position filtering
            mask_combined = np.isin(self['variants_pos'], pos)

        # Create mask based on indexes if provided
        if indexes is not None:
            # Validate indexes, allowing negative indexes
            out_of_bounds_indexes = indexes[(indexes < -n_snps) | (indexes >= n_snps)]
            if out_of_bounds_indexes.size > 0:
                raise ValueError(f"One or more sample indexes are out of bounds.")

            # Handle negative indexes and check for out-of-bounds indexes
            adjusted_indexes = np.mod(indexes, n_snps)

            # Create mask for specified indexes
            mask_indexes = np.zeros(n_snps, dtype=bool)
            mask_indexes[adjusted_indexes] = True

            # Combine with `chrom` and `pos` mask using logical OR (union of all specified criteria)
            mask_combined = mask_combined | mask_indexes

        # Invert mask if `include` is False
        if not include:
            mask_combined = ~mask_combined

        # Define keys to filter
        keys = [
            'lai', 'calldata_gt', 'variants_ref', 'variants_alt', 
            'variants_chrom', 'variants_filter_pass', 'variants_id', 
            'variants_pos', 'variants_qual'
        ]

        # Apply filtering based on inplace parameter
        if inplace:
            for key in keys:
                if self[key] is not None:
                    self[key] = np.asarray(self[key])[mask_combined]
            return None
        else:
            snpobj = self.copy()
            for key in keys:
                if snpobj[key] is not None:
                    snpobj[key] = np.asarray(snpobj[key])[mask_combined]
            return snpobj

    def filter_samples(
        self, 
        samples: Optional[Union[str, Sequence[str], np.ndarray, None]] = None,
        indexes: Optional[Union[int, Sequence[int], np.ndarray, None]] = None,
        include: bool = True,
        inplace: bool = False
    ) -> Optional['SNPLevelAncestryObject']:
        """
        Filter samples based on specified names or indexes.

        This method updates the `samples`, `calldata_gt`, `haplotypes`, and `lai` attributes 
        to include or exclude the specified samples. The order of the samples is preserved.

        If both samples and indexes are provided, any sample matching either a name in samples 
        or an index in indexes will be included or excluded.

        This method allows inclusion or exclusion of specific samples by their names or indexes. 
        When both sample names and indexes are provided, the union of the specified samples is used. 
        Negative indexes are supported and follow [NumPy's indexing conventions](https://numpy.org/doc/stable/user/basics.indexing.html).

        Args:
            samples (str or array_like of str, optional): 
                Name(s) of the samples to include or exclude. Can be a single sample name or a
                sequence of sample names. Default is None.
            indexes (int or array_like of int, optional):
                Index(es) of the samples to include or exclude. Can be a single index or a sequence
                of indexes. Negative indexes are supported. Default is None.
            include (bool, default=True): 
                If True, includes only the specified samples. If False, excludes the specified samples. Default is True.
            inplace (bool, default=False): 
                If True, modifies `self` in place. If False, returns a new `SNPLevelAncestryObject` 
                with the samples filtered. Default is False.

        Returns:
            **Optional[SNPLevelAncestryObject]:** 
                A new `SNPLevelAncestryObject` with the specified samples filtered if `inplace=False`. 
                If `inplace=True`, modifies `self` in place and returns None.
        """
        if samples is None and indexes is None:
            raise ValueError("At least one of 'samples' or 'indexes' must be provided.")

        n_samples = self.n_samples
        sample_names = np.array(self['samples'])

        # Create mask based on sample names
        if samples is not None:
            samples = np.atleast_1d(samples)
            mask_samples = np.isin(sample_names, samples)
            missing_samples = samples[~np.isin(samples, sample_names)]
            if missing_samples.size > 0:
                raise ValueError(f"The following specified samples were not found: {missing_samples.tolist()}")
        else:
            mask_samples = np.zeros(n_samples, dtype=bool)

        # Create mask based on sample indexes
        if indexes is not None:
            indexes = np.atleast_1d(indexes)

            # Validate indexes, allowing negative indexes
            out_of_bounds_indexes = indexes[(indexes < -n_samples) | (indexes >= n_samples)]
            if out_of_bounds_indexes.size > 0:
                raise ValueError(f"One or more sample indexes are out of bounds.")
            
            # Handle negative indexes
            adjusted_indexes = np.mod(indexes, n_samples)

            mask_indexes = np.zeros(n_samples, dtype=bool)
            mask_indexes[adjusted_indexes] = True
        else:
            mask_indexes = np.zeros(n_samples, dtype=bool)

        # Combine masks using logical OR (union of samples)
        mask_combined = mask_samples | mask_indexes

        if not include:
            mask_combined = ~mask_combined

        # Filter `samples`
        filtered_samples = sample_names[mask_combined].tolist()

        # Filter `calldata_gt`
        filtered_calldata_gt = np.array(self['calldata_gt'])[:, mask_combined, ...]

        # Filter `haplotypes` (2 haplotypes per sample)
        haplotype_mask = np.repeat(mask_combined, 2)
        filtered_haplotypes = np.array(self['haplotypes'])[haplotype_mask].tolist()

        # Filter `lai`
        filtered_lai = np.array(self['lai'])[:, haplotype_mask]

        if inplace:
            self['samples'] = filtered_samples
            self['calldata_gt'] = filtered_calldata_gt
            self['haplotypes'] = filtered_haplotypes
            self['lai'] = filtered_lai
            return None
        else:
            snpobj = self.copy()
            snpobj['samples'] = filtered_samples
            snpobj['calldata_gt'] = filtered_calldata_gt
            snpobj['haplotypes'] = filtered_haplotypes
            snpobj['lai'] = filtered_lai
            return snpobj

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
