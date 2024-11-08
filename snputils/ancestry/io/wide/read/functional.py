from pathlib import Path
from typing import Union, Optional

from snputils.ancestry.genobj.wide import GlobalAncestryObject


def read_admixture(
    Q_file: Union[str, Path],
    P_file: Union[str, Path],
    sample_file: Optional[Union[str, Path]] = None,
    snp_file: Optional[Union[str, Path]] = None,
    ancestry_file: Optional[Union[str, Path]] = None,
) -> 'GlobalAncestryObject':
    """
    Read admixture files into a GlobalAncestryObject.

    This function processes the provided Q and P matrix files to extract per-sample ancestry proportions
    and per-ancestry SNP frequencies. Optional files for sample identifiers, SNP identifiers, and 
    ancestry labels are also read, if provided.

    Args:
        Q_file (str or pathlib.Path):
            Path to the file containing the Q matrix (per-sample ancestry proportions).
        P_file (str or pathlib.Path):
            Path to the file containing the P/F matrix (per-ancestry SNP frequencies).
        sample_file (str or pathlib.Path, optional):
            Path to the file containing sample identifiers. If None, sample identifiers are not loaded.
        snp_file (str or pathlib.Path, optional):
            Path to the file containing SNP identifiers. If None, SNP identifiers are not loaded.
        ancestry_file (str or pathlib.Path, optional):
            Path to the file containing ancestry labels for each sample. If None, ancestries are not loaded.

    Returns:
            GlobalAncestryObject: 
                A global ancestry object instance.
    """
    from snputils.ancestry.io.wide.read.admixture import AdmixtureReader

    return AdmixtureReader(
        Q_file=Q_file,
        P_file=P_file,
        sample_id_file=sample_file,
        snp_file=snp_file,
        ancestry_file=ancestry_file
    ).read()


read_adm = read_admixture
