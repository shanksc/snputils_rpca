import logging
import numpy as np
import polars as pl
import pgenlib as pg
from pathlib import Path

from snputils.snp.genobj.snpobj import SNPObject

log = logging.getLogger(__name__)


class PGENWriter:
    """
    Writes a genotype object in PGEN format (.pgen, .psam, and .pvar files) in the specified output path.
    """

    def __init__(self, snpobj: SNPObject, filename: str):
        """
        Initializes the PGENWriter instance.

        Parameters
        ----------
        snpobj : SNPObject
            The SNPObject containing genotype data to be written.
        file : str
            Base path for the output files (excluding extension).
        TODO: add support for parallel writing by chromosome.
        """
        self.__snpobj = snpobj
        self.__filename = Path(self.__filename)

    def write(self):
        """
        Writes the SNPObject data to .pgen, .psam, and .pvar files.
        """
        file_extensions = (".pgen", ".psam", ".pvar")
        if self.__filename.suffix in file_extensions:
            self.__filename = self.__filename.with_suffix('')
        self.__file_extension = ".pgen"

        self.write_pvar()
        self.write_psam()
        self.write_pgen()

    def write_pvar(self):
        """
        Writes variant data to the .pvar file.
        """
        log.info(f"Writing to {self.__filename}.pvar")
        df = pl.DataFrame(
            {
                "#CHROM": self.__snpobj.variants_chrom,
                "POS": self.__snpobj.variants_pos,
                "ID": self.__snpobj.variants_id,
                "REF": self.__snpobj.variants_ref,
                "ALT": self.__snpobj.variants_alt[:, 0],
                "FILTER": self.__snpobj.variants_filter_pass,
                # TODO: add INFO column to SNPObject and write it to the .pvar file? (if not it's lost)
            }
        )
        # TODO: add header to the .pvar file, if not it's lost
        df.write_csv(f"{self.__filename}.pvar", separator="\t")

    def write_psam(self):
        """
        Writes sample metadata to the .psam file.
        """
        log.info(f"Writing {self.__filename}.psam")
        df = pl.DataFrame(
            {
                "IID": self.__snpobj.samples,
                "SEX": "NA",  # Add SEX as nan for now
                # TODO: add SEX as Optional column to SNPObject and write it to the .psam file (if not it's lost)
            }
        )
        df.write_csv(f"{self.__filename}.psam", separator="\t")

    def write_pgen(self):
        """
        Writes the genotype data to a .pgen file.
        """
        log.info(f"Writing to {self.__filename}.pgen")
        phased = True if self.__snpobj.calldata_gt.ndim == 3 else False
        if phased:
            num_variants, num_samples, num_alleles = self.__snpobj.calldata_gt.shape
            # Flatten the genotype matrix for pgenlib
            flat_genotypes = self.__snpobj.calldata_gt.reshape(
                num_variants, num_samples * num_alleles
            )
        else:
            num_variants, num_samples = self.__snpobj.calldata_gt.shape
            flat_genotypes = self.__snpobj.__calldata_gt

        with pg.PgenWriter(
            f"{self.__filename}.pgen".encode("utf-8"),
            sample_ct=num_samples,
            variant_ct=num_variants,
            hardcall_phase_present=phased,
        ) as writer:
            for variant_index in range(num_variants):
                writer.append_alleles(
                    flat_genotypes[variant_index].astype(np.int32), all_phased=phased
                )
