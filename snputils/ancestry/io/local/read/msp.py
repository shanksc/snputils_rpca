import pathlib
from typing import List, Dict, Optional
import logging
import warnings
import numpy as np
import pandas as pd
from typing import Union

from .base import LAIBaseReader
from snputils.ancestry.genobj.local import LocalAncestryObject

log = logging.getLogger(__name__)


class MSPReader(LAIBaseReader):
    """
    A reader class for parsing Local Ancestry Inference (LAI) data from `.msp` or `msp.tsv` files
    and constructing a `snputils.ancestry.genobj.LocalAncestryObject` that contains haplotype-level
    ancestry assignments across genomic windows.

    **Expected MSP content:**

    - **Chromosome numbers** (`#chm`): Chromosome numbers corresponding to each genomic window.
    - **Physical positions** (`spos`, `epos`): Start and end physical positions (base pairs) for each window (if available).
    - **Genetic positions** (`sgpos`, `egpos`): Start and end centimorgan positions for each window (if available).
    - **Window sizes** (`n snps`): Number of SNPs in each genomic window (if available).
    - **Haplotype identifiers**: Column headers for LAI data columns, starting from the first LAI column.
      Each haplotype is typically denoted by a sample identifier with a `.0` (paternal) or `.1` (maternal) suffix.
    - **LAI data**: Haplotype-level ancestry information for each genomic window.
    """
    def __init__(self, file: Union[str, pathlib.Path]) -> None:
        """
        Args:
            file (str or pathlib.Path): 
                Path to the file to be read. It should end with `.msp` or `.msp.tsv`.
        """
        self.__file = pathlib.Path(file)

    @property
    def file(self) -> pathlib.Path:
        """
        Retrieve `file`.

        Returns:
            **file (str or pathlib.Path):** 
                Path to the file to be read. It should end with `.msp` or `.msp.tsv`.
        """
        return self.__file

    def _get_samples(self, msp_df: pd.DataFrame, first_lai_col_indx: int) -> List[str]:
        """
        Extract unique sample identifiers from the pandas DataFrame.

        Args:
            msp_df (pd.DataFrame): 
                The DataFrame representing the `.msp` data, including LAI columns.
            first_lai_col_indx (int): 
                Index of the first column containing LAI data.

        Returns:
            **list:** List of unique sample identifiers.
        """
        # Get all columns starting from the first LAI data column
        query_samples_dub = msp_df.columns[first_lai_col_indx:]

        # Select only one of the maternal/paternal samples by taking every second sample
        single_ind_idx = np.arange(0, len(query_samples_dub), 2)
        query_samples_sing = query_samples_dub[single_ind_idx]

        # Remove the suffix from sample names to get clean identifiers
        query_samples = [qs[:-2] for qs in query_samples_sing]

        return query_samples

    def _get_ancestry_map_from_comment(self, comment: str) -> Dict[str, str]:
        """
        Construct an ancestry map from the comment line of the `.msp` file.

        This method parses the comment string to create a mapping of ancestry numerical identifiers 
        to their corresponding ancestry names.

        Args:
            comment (str): 
                The comment line containing ancestry mapping information.

        Returns:
            dict: A dictionary mapping names to ancestry numbers (as strings).
        """
        # Split the comment to extract relevant ancestry mapping data
        comment = comment.split(' ')[-1].replace('\n', '').split('\t')
        ancestry_map = {}

        # Populate the ancestry map with parsed values
        for elem in comment:
            x, y = elem.split('=')
            ancestry_map[y] = x

        return ancestry_map

    def read(self) -> 'LocalAncestryObject':
        """
        Read data from the provided `file` and construct a `snputils.ancestry.genobj.LocalAncestryObject`.

        Returns:
            **LocalAncestryObject:**
                A LocalAncestryObject instance.
        """
        log.info(f"Reading '{self.file}'...")

        # Open the file and read the first two lines to identify the header and comment
        with open(self.file) as f:
            first_line = f.readline()
            second_line = f.readline()

        first_line_ = [h.replace('\n', '') for h in first_line.split("\t")]
        second_line_ = [h.replace('\n', '') for h in second_line.split("\t")]

        # Determine which line contains the #chm reference value for the header
        if "#chm" in first_line_:
            comment = None
            header = first_line_
        elif "#chm" in second_line_:
            comment = first_line
            header = second_line_
        else:
            raise ValueError(
                f"Header not found. Expected '#chm' in the first two lines. "
                f"First line: {first_line.strip()} | Second line: {second_line.strip()}"
            )

        # Ensure no duplicate columns exist in the header
        if len(header) != len(set(header)):
            raise ValueError("Duplicate columns detected in the header.")

        # Read the main data into a DataFrame, skipping comment lines
        msp_df = pd.read_csv(self.file, sep="\t", comment="#", names=header)

        # Extract chromosomes data
        chromosomes = msp_df['#chm'].to_numpy()

        # Extract physical positions (if available)
        column_counter = 1
        try:
            physical_pos = msp_df[['spos', 'epos']].to_numpy()
            column_counter += 2
        except KeyError:
            physical_pos = None
            log.warning("Physical positions ('spos' and 'epos') not found.")
        
        # Extract centimorgan positions (if available)
        try:
            centimorgan_pos = msp_df[['sgpos', 'egpos']].to_numpy()
            column_counter += 2
        except KeyError:
            centimorgan_pos = None
            log.warning("Genetic (centimorgan) positions ('sgpos' and 'egpos') not found.")
        
        # Extract window sizes (if available)
        try:
            window_sizes = msp_df['n snps'].to_numpy()
            column_counter += 1
        except KeyError:
            window_sizes = None
            log.warning("Window sizes ('n snps') not found.")
        
        # Extract LAI data (haplotype-level)
        lai = msp_df[msp_df.columns[column_counter:]].to_numpy()

        # Extract haplotype identifiers
        haplotypes = msp_df.columns[column_counter:].to_list()

        # Extract haplotype identifiers and sample identifiers
        samples = self._get_samples(msp_df, column_counter)

        # Validate the number of samples matches the LAI data dimensions
        n_samples = len(samples)
        if n_samples != int(lai.shape[1] / 2):
            raise ValueError(
                "Mismatch between the number of sample identifiers and the expected number of samples in the LAI array. "
                f"Expected {int(lai.shape[1] / 2)} samples (derived from LAI data); found {n_samples}."
            )
        
        # Count number of unique ancestries in the LAI data
        n_ancestries = len(np.unique(lai))

        # Parse ancestry map from the comment (if available)
        ancestry_map = None
        if comment is not None:
            ancestry_map = self._get_ancestry_map_from_comment(comment)
            if len(ancestry_map) != n_ancestries:
                warnings.warn(
                    "Mismatch between the number of unique ancestries in the LAI data "
                    f"({n_ancestries}) and the number of classes in the ancestry map "
                    f"({len(ancestry_map)})."
                )
        else:
            # Provide default ancestry mapping if no comment is provided
            ancestry_map = None
            warnings.warn(
                "Ancestry map not found. It is recommended to provide an .msp file that contains the ancestry "
                "map as a comment in the first line."
            )

        return LocalAncestryObject(
            haplotypes=haplotypes,
            lai=lai,
            samples=samples,
            ancestry_map=ancestry_map,
            window_sizes=window_sizes,
            centimorgan_pos=centimorgan_pos,
            chromosomes=chromosomes,
            physical_pos=physical_pos
        )

LAIBaseReader.register(MSPReader)
