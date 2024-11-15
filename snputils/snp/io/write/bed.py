import logging
import pandas as pd
import numpy as np
import pgenlib as pg
from pathlib import Path

from snputils.snp.genobj import SNPObject

log = logging.getLogger(__name__)


class BEDWriter:
    """Writes an object in bed/bim/fam formats in the specified output path.

    Args:
        snpobj: The SNPObject to be written.
        file: The output file path.

    """

    def __init__(self, snpobj: SNPObject, filename: str):
        self.__snpobj = snpobj.copy()
        self.__filename = Path(filename)

    def write(self):
        """Writes the SNPObject to bed/bim/fam formats."""

        # Save .bed file
        if self.__filename.suffix != '.bed':
            self.__filename = self.__filename.with_suffix('.bed')

        log.info(f"Writing .bed file: {self.__filename}")

        # If the input matrix has three dimensions, it indicates that the data is divided into two strands.
        if len(self.__snpobj.calldata_gt.shape) == 3:
            # Sum the two strands
            self.__snpobj.calldata_gt = self.__snpobj.calldata_gt.transpose(1, 0, 2).sum(axis=2, dtype=np.int8)

        # Infer the number of samples and variants from the matrix
        samples, variants = self.__snpobj.calldata_gt.shape

        # Define the PgenWriter to save the data
        data_save = pg.PgenWriter(filename=str(self.__filename).encode('utf-8'),
                                  sample_ct=samples,
                                  variant_ct=variants,
                                  nonref_flags=True,
                                  hardcall_phase_present=False,
                                  dosage_present=True,
                                  dosage_phase_present=False)

        # Fill the data_save object with the matrix of individuals x variants
        for snp_i in range(0, variants):
            data_save.append_biallelic(np.ascontiguousarray(self.__snpobj.calldata_gt[:, snp_i]))

        # Save the .bed file
        data_save.close()

        log.info(f"Finished writing .bed file: {self.__filename}")

        # Remove .bed from the file name
        if self.__filename.suffix == '.bed':
            self.__filename = self.__filename.with_suffix('')

        # Save .fam file
        log.info(f"Writing .fam file: {self.__filename}")

        # Fill .fam file
        fam_file = pd.DataFrame(columns=['fid', 'iid', 'father', 'mother', 'gender', 'trait'])
        fam_file['iid'] = self.__snpobj.samples
        fam_file['fid'] = self.__snpobj.samples

        # Save .fam file
        fam_file.to_csv(self.__filename.with_suffix('.fam'), sep='\t', index=False, header=False)
        log.info(f"Finished writing .fam file: {self.__filename}")

        # Save .bim file
        log.info(f"Writing .bim file: {self.__filename}")

        # Fill .bim file
        bim_file = pd.DataFrame(columns=['chrom', 'snp', 'cm', 'pos', 'a0', 'a1'])
        bim_file['chrom'] = self.__snpobj.variants_chrom
        bim_file['snp'] = self.__snpobj.variants_id
        bim_file['cm'] = 0  # TODO: read, save and write too if available?
        log.warning("The .bim file is being saved with 0 cM values.")
        bim_file['pos'] = self.__snpobj.variants_pos
        bim_file['a0'] = self.__snpobj.variants_alt[:, 0]
        bim_file['a1'] = self.__snpobj.variants_ref

        # Save .bim file
        bim_file.to_csv(self.__filename.with_suffix('.bim'), sep='\t', index=False, header=False)
        log.info(f"Finished writing .bim file: {self.__filename}")
