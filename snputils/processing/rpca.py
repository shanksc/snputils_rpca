import numpy as np
import copy
from typing import Tuple, Optional, Union, List

from snputils.snp.genobj.snpobj import SNPObject
from ._utils.srpca import *

class RPCA: 
    """
    A class for Robust Principal Component Analysis (RPCA) on SNP data stored in `snputils.snp.genobj.SNPObject`.
    This class uses nonconvex principal component pursuit for recovering a low-rank matrix from missing 
    and corrupted data. 
    """

    def __init__(
        self, 
        snpobj: Optional['SNPObject'] = None, 
        rank: int = 5, 
        average_strands: bool = False, 
        samples_subset: Optional[Union[int, List]] = None, 
        snps_subset: Optional[Union[int, List]] = None,
        method: Optional[str] = 'nonconvex_srpcp',
        max_iters: int = 500,
        verbose: bool = True,
        eps_abs: float = 1e-4, 
        eps_rel: float = 1e-4 
    ):
        """
        Args:
            snpobj (SNPObject, optional): 
                A SNPObject instance.
            rank (int, default=5): 
                The number of principal components. Default is 5.
            average_strands (bool, default=True): 
                True if the haplotypes from the two parents are to be combined (averaged) for each individual, or False otherwise.
            samples_subset (int or list of int, optional): 
                Subset of samples to include, as an integer for the first samples or a list of sample indices.
            snps_subset (int or list of int, optional): 
                Subset of SNPs to include, as an integer for the first SNPs or a list of SNP indices.
        """
        self.__snpobj = snpobj
        self.__rank = rank
        self.__average_strands = average_strands
        self.__samples_subset = samples_subset
        self.__snps_subset = snps_subset
        self.__method = method 
        self.__X_ = None
        self.__L_ = None  # Store low rank matrix 
        self.__S_ = None # Store sparse matrix 
        self.__max_iters = max_iters
        self.__verbose = verbose 
        self.__eps_abs = eps_abs
        self.__eps_rel = eps_rel
        
        # Handle sample and SNP subsets
        if isinstance(samples_subset, int):
            X = X[:samples_subset]
        elif isinstance(samples_subset, list):
            X = X[samples_subset]
        
        if isinstance(snps_subset, int):
            X = X[:, :snps_subset]
        elif isinstance(snps_subset, list):
            X = X[:, snps_subset] 
            
        if self.snpobj is not None:
            self.fit_transform(snpobj)


    @property
    def snpobj(self) -> Optional['SNPObject']:
        """
        Retrieve `snpobj`.

        Returns:
            **SNPObject:** A SNPObject instance.
        """
        return self.__snpobj

    @snpobj.setter
    def snpobj(self, x) -> None:
        """
        Update `snpobj`.
        """
        self.__snpobj = x

    @property
    def rank(self) -> int:
        """
        Retrieve `rank`.
        
        Returns:
            **int:** rank for srpcp.
        """
        return self.__rank
    
    @rank.setter
    def rank(self, x: int) -> None:
        """
        Update `rank`.
        """
        self.__rank = x

    @property
    def average_strands(self) -> bool:
        """
        Retrieve `average_strands`.

        Returns:
            **bool:** 
                True if the haplotypes from the two parents are to be combined (averaged) for each individual, or False otherwise.
        """
        return self.__average_strands

    @average_strands.setter
    def average_strands(self, x: bool) -> None:
        """
        Update `average_strands`.
        """
        self.__average_strands = x

    @property
    def samples_subset(self) -> Optional[Union[int, List[int]]]:
        """
        Retrieve `samples_subset`.

        Returns:
            **int or list of int:** 
                Subset of samples to include, as an integer for the first samples or a list of sample indices.
        """
        return self.__samples_subset

    @samples_subset.setter
    def samples_subset(self, x: Optional[Union[int, List[int]]]) -> None:
        """
        Update `samples_subset`.
        """
        self.__samples_subset = x

    @property
    def snps_subset(self) -> Optional[Union[int, List[int]]]:
        """
        Retrieve `snps_subset`.

        Returns:
            **int or list of int:** Subset of SNPs to include, as an integer for the first SNPs or a list of SNP indices.
        """
        return self.__snps_subset

    @snps_subset.setter
    def snps_subset(self, x: Optional[Union[int, List[int]]]) -> None:
        """
        Update `snps_subset`.
        """
        self.__snps_subset = x

    @property
    def method(self) -> str:
        """
        Retrieve `method`. 

        Returns: 
            **str:** RPCA or matrix recovery method being used. 
        """
        return self.__method

    @method.setter 
    def method(self, x: Optional[str]) -> None:
        """
        Update `method`. 
        """
        self.__method = x 

    @property
    def rank_(self) -> Optional[int]:
        """
        Retrieve `n_components_`.

        Returns:
            **int:** The rank to use for srpcp
        """
        return self.__rank_

    @property
    def X_(self) -> Optional[np.ndarray]:
        """
        Retrieve `X_`.

        Returns:
            **array of shape (n_samples, n_snps):** 
                The SNP data matrix used to fit the model.
        """
        return self.__X_

    @X_.setter
    def X_(self, x: np.ndarray) -> None:
        """
        Update `X_`.
        """
        self.__X_ = x

    @property
    def L_(self) -> Optional[np.ndarray]:
        """
        Retrieve `L_`.

        Returns:
            **array of shape (n_samples, n_snps):** 
                The low rank SNP data matrix.
        """
        return self.__L_

    @L_.setter
    def L_(self, x: np.ndarray) -> None:
        """
        Update `X_new_`.
        """
        self.__L_ = x 

    @property
    def S_(self) -> Optional[np.ndarray]:
        """
        Retrieve `S_`.

        Returns:
            **array of shape (n_samples, n_snps):** 
                The sparse SNP data matrix.
        """
        return self.__S_

    @S_.setter
    def S_(self, x: np.ndarray) -> None:
        """
        Update `X_new_`.
        """
        self.__S_ = x

    @property
    def max_iters(self) -> int:
        """
        Retrieve `max_iters`.

        Returns:
            **int:** Maximum number of iterations for the RPCA algorithm.
        """
        return self.__max_iters

    @max_iters.setter
    def max_iters(self, x: int) -> None:
        """
        Update `max_iters`.

        Args:
            x (int): The new maximum number of iterations for the RPCA algorithm.
        """
        self.__max_iters = x

    @property
    def verbose(self) -> bool:
        """
        Retrieve `verbose`.

        Returns:
            **bool:** Whether to print detailed logs during RPCA computation.
        """
        return self.__verbose

    @verbose.setter
    def verbose(self, x: bool) -> None:
        """
        Update `verbose`.

        Args:
            x (bool): Whether to enable verbose logging.
        """
        self.__verbose = x

    @property
    def eps_abs(self) -> float:
        """
        Retrieve `eps_abs`.

        Returns:
            **float:** Absolute tolerance for convergence in the RPCA algorithm.
        """
        return self.__eps_abs

    @eps_abs.setter
    def eps_abs(self, x: float) -> None:
        """
        Update `eps_abs`.

        Args:
            x (float): New absolute tolerance for convergence.
        """
        self.__eps_abs = x

    @property
    def eps_rel(self) -> float:
        """
        Retrieve `eps_rel`.

        Returns:
            **float:** Relative tolerance for convergence in the RPCA algorithm.
        """
        return self.__eps_rel

    @eps_rel.setter
    def eps_rel(self, x: float) -> None:
        """
        Update `eps_rel`.

        Args:
            x (float): New relative tolerance for convergence.
        """
        self.__eps_rel = x


    def copy(self) -> 'RPCA':
        """
        Create and return a copy of `self`.

        Returns:
            **RPCA:** 
                A new instance of the current object.
        """
        return copy.copy(self)

    def _get_data_from_snpobj(
            self, 
            snpobj: Optional['SNPObject'] = None, 
            average_strands: Optional[bool] = None, 
            samples_subset: Optional[Union[int, List]] = None, 
            snps_subset: Optional[Union[int, List]] = None
        ) -> np.ndarray:
        """
        Retrieve and prepare SNP data for RPCA analysis, with options for selecting subsets and handling strands.

        This method processes SNP data stored in an `SNPObject`, which may include averaging of paternal 
        and maternal strands or selecting subsets of samples and SNPs.

        Args:
            snpobj (SNPObject, optional):
                A SNPObject object instance. If None, defaults to `self.snpobj`.
            average_strands (bool, optional):
                True if the haplotypes from the two parents are to be combined (averaged) for each individual, or False otherwise.
                If None, defaults to `self.average_strands`.
            samples_subset (int or list of int, optional):
                Subset of samples to include, as an integer for the first samples or a list of sample indices.
                If None, defaults to `self.samples_subset`.
            snps_subset (int or list of int, optional):
                Subset of SNPs to include, as an integer for the first SNPs or a list of SNP indices.
                If None, defaults to `self.snps_subset`.

        Returns:
            numpy.ndarray:
                The processed SNP data.
        """
        if snpobj is None:
            snpobj = self.snpobj
        if average_strands is None:
            average_strands = self.average_strands
        if samples_subset is None:
            samples_subset = self.samples_subset
        if snps_subset is None:
            snps_subset = self.snps_subset
            
        
        if snpobj.calldata_gt.ndim == 2:
            X = np.transpose(snpobj.calldata_gt.astype(float), (1,0))
        elif snpobj.calldata_gt.ndim == 3:
            X = np.transpose(snpobj.calldata_gt.astype(float), (1,0,2))
        
            if average_strands:
                X = np.mean(X, axis=2)
            else:
                X = np.reshape(X, (-1, X.shape[1]))
        else:
            raise ValueError(f"Invalid shape for `calldata_gt`: expected a 2D or 3D array, but got {snpobj.calldata_gt.ndim}D array.")
    



        # Handle sample and SNP subsets
        if isinstance(samples_subset, int):
            X = X[:samples_subset]
        elif isinstance(samples_subset, list):
            X = X[samples_subset]
        
        if isinstance(snps_subset, int):
            X = X[:, :snps_subset]
        elif isinstance(snps_subset, list):
            X = X[:, snps_subset]
        
        return X
        
    
    def nonconvex_srpcp(self, X_incomplete, rank, max_iters, verbose, eps_abs, eps_rel):
        omega = ~np.isnan(X_incomplete)
        n, p = X_incomplete.shape 
        D = np.nan_to_num(X_incomplete, nan=-1)
        
        if verbose:
            print(f'dimensions: {n}, {p}')

        rho = 0.1
        res_primal = 0
        res_dual = 0
        thresh_primal = 0
        thresh_dual = 0
        
        L1 = np.zeros((n, p))
        L2 = np.zeros((n, p))
        L3 = np.zeros((n, p))

        S1 = np.zeros((n, p))
        S2 = np.zeros((n, p))

        Z = np.zeros((n, p))
        Y1 = np.zeros((n, p))
        Y2 = np.zeros((n, p))
        Y3 = np.zeros((n, p))
        Y4 = np.zeros((n, p))

        mu = np.sqrt(p / 2)     
        lambd = 1 / np.sqrt(n)

        if verbose:
            print(f'mu: {mu:.5f}')
            print(f'lambda: {lambd:.5f}')

        converged = False

        for i in range(0, max_iters):
            if i % 50 == 0 and verbose:
                print(f"Current iteration {i}")
                print(f"res dual {res_dual}")
                print(f"res primal {res_primal}")
                print(f"thresh dual {thresh_dual}")
                print(f"thresh primal {thresh_primal}")
            
            L2_old = L2.copy()
            S2_old = S2.copy()
            L3_old = L3.copy()
            
            #update first primal variables 
            L1 = proj_rank_r((L2 + L3 - Y1 / rho - Y4 / rho) / 2, k=rank)
            S1 = prox_l1(S2 - Y2 / rho, lambd / rho)

            Z = prox_fro((((L2 + S2 - Y3 / rho) - D)), mu / rho) + D
            L3 = np.maximum(L1 + Y4 / rho, 0)

            #second primal variables 
            L2 = omega * (1/3 * (2 * L1 - S1 + Z + (2 * Y1 - Y2 + Y3) / rho))
            L2_na = (1 - omega) * (L1 + Y1 / rho)
            L2 += L2_na 

            S2 = omega * (1/3 * (2 * S1 - L1 + Z + (2 * Y2 - Y1 + Y3) / rho))
            S2_na = (1 - omega) * (S1 + Y2 / rho) 
            S2 += S2_na 

            #update dual variables 
            Y1 += rho * (L1 - L2)
            Y2 += rho * (S1 - S2)
            Y3 += rho * omega * (Z - (L2 + S2))
            Y4 += rho * (L1 - L3)

            #calc residuals
            res_primal = np.sqrt(np.linalg.norm(L1 - L2, 'fro')**2 + 
                            np.linalg.norm(S1 - S2, 'fro')**2 +
                            np.linalg.norm(omega * (Z - L2 - S2), 'fro')**2 +
                            np.linalg.norm(L1 - L3, 'fro')**2)
            res_dual = rho * np.sqrt(np.linalg.norm(L2 + L3 - L2_old - L3_old, 'fro')**2 + 
                            np.linalg.norm(S2 - S2_old, 'fro')**2 +
                            np.linalg.norm(omega * (L2 - L2_old + S2 - S2_old), 'fro')**2)
            

            if res_primal > 10 * res_dual:
                rho *= 2 
            elif res_dual > 10 * res_primal:
                rho /= 2

            #stopping criterias
            thresh_primal = (
                eps_abs * np.sqrt(4 * n * p) +
                eps_rel * max(
                    np.sqrt(2 * np.linalg.norm(L1, 'fro')**2 +
                            np.linalg.norm(S1, 'fro')**2 +
                            np.linalg.norm(Z, 'fro')**2),
                    np.sqrt(np.linalg.norm(L2, 'fro')**2 +
                            np.linalg.norm(S2, 'fro')**2 +
                            np.linalg.norm(omega * (L2 + S2), 'fro')**2 +
                            np.linalg.norm(L3, 'fro')**2)
                )
            )

            thresh_dual = (
                eps_abs * np.sqrt(3 * n * p) +
                eps_rel * np.sqrt(
                    np.linalg.norm(Y1 + Y4, 'fro')**2 +
                    np.linalg.norm(Y2, 'fro')**2 +
                    np.linalg.norm(Y3, 'fro')**2
                )
            )
            
            if res_primal < thresh_primal and res_dual < thresh_dual:
                converged = True
                print(f"Converged in {i+1} iterations.")
                break 

        L = (L1 + L2 + L3) / 3
        S = (S1 + S2) / 2  

        if not converged and verbose:
            print(f'Did not converge in {i+1} iterations')

        #clip everything
        clip = np.amax(D)
        print(clip)
        L = np.clip(L, 0, clip)
        S = np.clip(S, 0, clip)
        
        if self.average_strands:
            L = np.round(2 * L) / 2
            S = np.round(2 * S) / 2 
        else:
            L = np.round(L) 
            S = np.round(L)

        return (L, S) 



    def fit_transform(
            self,
            snpobj: Optional['SNPObject'] = None,
             average_strands: Optional[bool] = None, 
            samples_subset: Optional[Union[int, List]] = None, 
            snps_subset: Optional[Union[int, List]] = None
        ) -> 'RPCA':
        """
        Fit the model to the input SNP data stored in the provided `snpobj`.
        
        Args:
            snpobj (SNPObject, optional): 
                A SNPObject instance. If None, defaults to `self.snpobj`.
            average_strands (bool, optional): 
                True if the haplotypes from the two parents are to be combined (averaged) for each individual, or False otherwise.
                If None, defaults to `self.average_strands`.
            samples_subset (int or list of int, optional): 
                Subset of samples to include, as an integer for the first samples or a list of sample indices.
                If None, defaults to `self.samples_subset`.
            snps_subset (int or list of int, optional): 
                Subset of SNPs to include, as an integer for the first SNPs or a list of SNP indices.
                If None, defaults to `self.snps_subset`.
 
        Returns:
            **RPCA:**
                The fitted instance of `self`.
        """
        self.X_ = self._get_data_from_snpobj(snpobj, average_strands, samples_subset, snps_subset)
        
        #if self.method == 'nonconvex_srpcp':int
        L, S = self.nonconvex_srpcp(self.X_, self.rank, self.max_iters, self.verbose, self.eps_abs, self.eps_rel)
        
        self.L = L
        self.S = S 

        return L, S 



