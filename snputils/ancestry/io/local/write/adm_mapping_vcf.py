import logging
import numpy as np
from typing import Dict
from pathlib import Path

from snputils.ancestry.genobj.local import LocalAncestryObject
from snputils.snp.genobj.snpobj import SNPObject
from snputils.snp.io.write.vcf import VCFWriter

log = logging.getLogger(__name__)


class AdmixtureMappingVCFWriter:
    """
    A class that converts and writes Local Ancestry Inference (LAI) data into ancestry-specific 
    VCF files for ADMIXTURE mapping.
    """
    def __init__(
            self, 
            laiobj: LocalAncestryObject, 
            file: str, 
            ancestry_map: Dict
        ):
        """
        Args:
            laiobj (LocalAncestryObject): 
                A `LocalAncestryObject` instance.
            file (str): 
                Path to write VCF file.
            ancestry_map (dict): 
                A dictionary of str that represents the ancestry mapping
                as it is in the LocalAncestryObject.
        """
        self.__laiobj = laiobj
        self.__file = file
        self.__ancestry_map = ancestry_map

    def write(self):
        """
        Write VCF files for each ancestry type defined in the ancestry map.

        For each ancestry, this class maps LAI data to SNP alleles, using `1` to indicate positions 
        that match the specified ancestry and `0` for other ancestries at that position. Separate 
        VCF files are generated for each ancestry, allowing ancestry-specific genetic analysis in 
        a SNP-compatible format.

        The output consists of separate VCF files for each ancestry type, with filenames
        suffixed by the corresponding ancestry label, as defined in `ancestry_map`.

        If the specified `file` does not end with `.vcf` or `.bcf`, the `.vcf` extension will be appended.
        """
        # Process the list of positions to include both the start and end coordinates for each window
        # Iterate over each ancestry key in the ancestry mapping
        for key in self.__ancestry_map:
            ancestry = int(key)
            anc_string = self.__ancestry_map[key]

            # Define the output file format, ensuring it has the correct ancestry-specific suffix
            file_extension = (".vcf", ".bcf")
            file_path = Path(self.__file)
            
            # Check if file has one of the specified extensions
            if file_path.suffix not in file_extension:
                # If file does not have the correct extension, default to ".vcf"
                output_file = file_path.with_name(f"{file_path.stem}_{anc_string}.vcf")
            else:
                # If file has the correct extension, insert the ancestry string before the extension
                output_file = file_path.with_name(f"{file_path.stem}_{anc_string}{file_path.suffix}")

            # Format start and end positions for the VCF file
            if self.__laiobj.physical_pos is not None:
                pos_list = [f"{val1}_{val2}" for val1, val2 in self.__laiobj.physical_pos]
            else:
                pos_list = None

            # Modify LAI data values to simulate a SNP file
            # The positions in LAI corresponding to the current ancestry key are mapped to 1, and the rest to 0
            match = (self.__laiobj.lai == ancestry).astype(int)
            match = match.reshape(len(self.__laiobj.lai),int(len(self.__laiobj.lai[0])/2), 2 )

            # Set up VCF-related data
            calldata_gt = match
            samples = np.array(self.__laiobj.samples)
            variants_chrom = self.__laiobj.chromosomes
            variants_list = [str(i+1) for i in range(len(self.__laiobj.lai))]
            variants_id = np.array(variants_list)
            variants_ref = np.full(calldata_gt.shape[0], 'A', dtype='U5')
            variants_alt = np.full((calldata_gt.shape[0], 3), ['T', 'G', 'C'], dtype='U1')

            # Create the SNPObject
            variant_data_obj = SNPObject(
                calldata_gt=calldata_gt,
                samples=samples,
                variants_chrom=variants_chrom,
                variants_id=variants_id,
                variants_ref = variants_ref,
                variants_alt = variants_alt,
                variants_pos = pos_list,
            )

            # Log the start of the VCF file writing process
            log.info(f"Writing VCF file for ancestry '{anc_string}' to '{output_file}'...")

            # Write the VCF file
            vcf_writer = VCFWriter(variant_data_obj, output_file)
            vcf_writer.write()

            log.info(f"Finished writing VCF file for ancestry '{anc_string}' to '{output_file}'...")

        return