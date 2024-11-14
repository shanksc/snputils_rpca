import logging
from pathlib import Path
from typing import Union
import pandas as pd

from .base import LAIBaseWriter
from snputils.ancestry.genobj.local import LocalAncestryObject

log = logging.getLogger(__name__)


class MSPWriter(LAIBaseWriter):
    """
    A class for writing data stored in a `LocalAncestryObject` instance into an `.msp` file.
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
    
    def write(self):
        """
        Write the data contained in the `laiobj` instance into an `.msp` file.

        This method constructs the `.msp` file where each row corresponds to a genomic 
        window and includes the following columns:

            - `#chm`: Chromosome number corresponding to each window.
            - `spos`: Start physical position for each window.
            - `epos`: End physical position for each window.
            - `sgpos`: Start centimorgan position for each window.
            - `egpos`: End centimorgan position for each window.
            - `n snps`: Number of SNPs within each genomic window.
            - `SampleID.0`: Local ancestry for the first haplotype of the sample for each window.
            - `SampleID.1`: Local ancestry for the second haplotype of the sample for each window.
        
        If the specified file does not end with `.msp`, the extension will be appended.
        """
        log.info(f"LAI object contains: {self.laiobj.n_samples} samples, {self.laiobj.n_ancestries} ancestries.")

        # Define the required file extension for the output
        file_extension = ".msp"

        # Append '.msp' extension to __file if not already present
        if not self.__file.name.endswith(file_extension):
            self.__file = self.__file.with_name(self.__file.name + file_extension)
        
        # Prepare columns for the DataFrame
        columns = ["spos", "epos", "sgpos", "egpos", "n snps"]
        lai_dic = {
            "#chm" : self.laiobj.chromosomes,
            "spos" : self.laiobj.physical_pos[:,0],
            "epos" : self.laiobj.physical_pos[:,1],
            "sgpos" : self.laiobj.centimorgan_pos[:,0],
            "egpos" : self.laiobj.centimorgan_pos[:,1],
            "n snps" : self.laiobj.window_sizes
        }

         # Populate the dictionary with sample data
        ilai = 0
        for ID in self.laiobj.samples:
            # First haplotype
            lai_dic[ID+".0"] = self.laiobj.lai[:,ilai]
            columns.append(ID+".0")
            # Second haplotype
            lai_dic[ID+".1"] = self.laiobj.lai[:,ilai+1]   
            columns.append(ID+".1")
            # Move to the next pair of haplotypes
            ilai += 2
            
        # Create a DataFrame from the dictionary containing all data
        lai_df = pd.DataFrame(lai_dic)

        log.info(f"Writing MSP file to '{self.file}'...")

        # Save the DataFrame to the .msp file in tab-separated format
        lai_df.to_csv(self.__file, sep="\t", index=False, header=False)
        
        # Construct the second line for the output file containing the column headers
        second_line = "#chm" + "\t" + "\t".join(columns)
        
        # If an ancestry map is available, prepend it to the output file
        if self.laiobj.ancestry_map is not None:
            ancestries_codes = list(self.laiobj.ancestry_map.keys()) # Get corresponding codes
            ancestries = list(self.laiobj.ancestry_map.values()) # Get ancestry names
            
            # Create the first line for the ancestry information, detailing subpopulation codes
            first_line = "#Subpopulation order/codes: " + "\t".join(
                f"{a}={ancestries_codes[ai]}" for ai, a in enumerate(ancestries)
            )

            # Open the file for reading and prepend the first line       
            with open(self.__file, "r+") as f:
                content = f.read()
                f.seek(0,0)
                f.write(first_line.rstrip('\r\n') + '\n' + second_line + '\n' + content)

        log.info(f"Finished writing MSP file to '{self.file}'.")

        return None

LAIBaseWriter.register(MSPWriter)
