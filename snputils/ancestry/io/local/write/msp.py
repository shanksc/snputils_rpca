import logging
from pathlib import Path
from typing import Union
import pandas as pd

from .base import LAIBaseWriter
from snputils.ancestry.genobj.local import LocalAncestryObject

log = logging.getLogger(__name__)


class MSPWriter(LAIBaseWriter):
    """
    A writer class for exporting local ancestry data from a `snputils.ancestry.genobj.LocalAncestryObject` 
    into an `.msp` or `.msp.tsv` file.
    """
    def __init__(self, laiobj: LocalAncestryObject, file: Union[str, Path]) -> None:
        """
        Args:
            laiobj (LocalAncestryObject):
                A local ancestry object instance.
            file (str or pathlib.Path): 
                Path to the file where the data will be saved. It should end with `.msp` or `.msp.tsv`. 
                If the provided path does not have one of these extensions, the `.msp` extension will be appended.
        """
        self.__laiobj = laiobj
        self.__file = Path(file)

    @property
    def laiobj(self) -> LocalAncestryObject:
        """
        Retrieve `laiobj`. 

        Returns:
            **LocalAncestryObject:** 
                A local ancestry object instance.
        """
        return self.__laiobj

    @property
    def file(self) -> Path:
        """
        Retrieve `file`.

        Returns:
            **pathlib.Path:** 
                Path to the file where the data will be saved. It should end with `.msp` or `.msp.tsv`. 
                If the provided path does not have one of these extensions, the `.msp` extension will be appended.
        """
        return self.__file
    
    def write(self):
        """
        Write the data contained in the `laiobj` instance to the specified output `file`. 
        If the file already exists, it will be overwritten.

        **Output MSP content:**

        The output `.msp` file will contain local ancestry assignments for each haplotype across genomic windows.
        Each row corresponds to a genomic window and includes the following columns:

        - `#chm`: Chromosome numbers corresponding to each genomic window.
        - `spos`: Start physical position for each window (if available).
        - `epos`: End physical position for each window (if available).
        - `sgpos`: Start centimorgan position for each window (if available).
        - `egpos`: End centimorgan position for each window (if available).
        - `n snps`: Number of SNPs in each genomic window (if available).
        - `SampleID.0`: Local ancestry for the first haplotype of the sample for each window.
        - `SampleID.1`: Local ancestry for the second haplotype of the sample for each window.
        """
        log.info(f"LAI object contains: {self.laiobj.n_samples} samples, {self.laiobj.n_ancestries} ancestries.")

        # Define the required file extensions
        valid_extensions = ('.msp', '.msp.tsv')

        # Append '.msp' extension if not already present
        if not self.file.name.endswith(valid_extensions):
            self.file = self.file.with_name(self.file.name + '.msp')
        
        # Initialize the columns and data dictionary for the DataFrame
        columns = []
        lai_dic = {}

        # Add chromosome numbers
        lai_dic["#chm"] = self.laiobj.chromosomes
        columns.append("#chm")

        # Add physical positions if available
        if self.laiobj.physical_pos is not None:
            lai_dic["spos"] = self.laiobj.physical_pos[:, 0]
            lai_dic["epos"] = self.laiobj.physical_pos[:, 1]
            columns.extend(["spos", "epos"])
        else:
            log.warning("Physical positions ('spos' and 'epos') are not available in the LAI object.")

        # Add genetic positions if available
        if self.laiobj.centimorgan_pos is not None:
            lai_dic["sgpos"] = self.laiobj.centimorgan_pos[:, 0]
            lai_dic["egpos"] = self.laiobj.centimorgan_pos[:, 1]
            columns.extend(["sgpos", "egpos"])
        else:
            log.warning("Genetic positions ('sgpos' and 'egpos') are not available in the LAI object.")

        # Add window sizes if available
        if self.laiobj.window_sizes is not None:
            lai_dic["n snps"] = self.laiobj.window_sizes
            columns.append("n snps")
        else:
            log.warning("Window sizes ('n snps') are not available in the LAI object.")

        # Add haplotype-level LAI data
        ilai = 0
        for ID in self.laiobj.samples:
            # First haplotype
            lai_dic[f"{ID}.0"] = self.laiobj.lai[:, ilai]
            columns.append(f"{ID}.0")
            # Second haplotype
            lai_dic[f"{ID}.1"] = self.laiobj.lai[:, ilai + 1]
            columns.append(f"{ID}.1")
            # Move to the next pair of haplotypes
            ilai += 2

        # Create a DataFrame from the dictionary containing all data
        lai_df = pd.DataFrame(lai_dic, columns=columns)

        log.info(f"Writing '{self.file}'...")

        # Save the DataFrame to the .msp file in tab-separated format without header and index
        lai_df.to_csv(self.file, sep="\t", index=False, header=False)
        
        # Construct the header line for the output file containing the column headers
        header_line = "\t".join(columns)
        
        # Open the file for reading and prepend the header lines
        with open(self.file, "r+") as f:
            content = f.read()
            f.seek(0, 0)

            # If an ancestry map is available, include it as the first line
            if self.laiobj.ancestry_map is not None:
                ancestries_codes = list(self.laiobj.ancestry_map.keys())  # Get ancestry codes
                ancestries = list(self.laiobj.ancestry_map.values())      # Get ancestry names

                # Create the first line with ancestry mapping
                first_line = "#Subpopulation order/codes: " + "\t".join(
                    f"{ancestries[ai]}={ancestries_codes[ai]}" for ai in range(len(ancestries))
                )
                f.write(first_line.rstrip('\r\n') + '\n' + header_line + '\n' + content)
            else:
                # If no ancestry map, only add the header line
                f.write(header_line + '\n' + content)

        log.info(f"Finished writing '{self.file}'.")

        return None

LAIBaseWriter.register(MSPWriter)
