{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module LAPACK (
  -- * General functions
  globalFAbstol
  -- * Utility arguments
  , JOBZ
  , jzEigvals
  , jzEigvalsEigvecs

  , RANGE
  , rngAll
  , rngInterval
  , rngEigNums

  , UPLO
  , uploUpper
  , uploLower
  -- * Eigenproblem methods
  , hszheevr
  , hszhbevx
  ) where

import           Data.Complex.Utils

import           Foreign                      (Storable, alloca, with)
import           Foreign.C.String             (castCharToCChar,withCAString)
import           Foreign.C.Types              (CChar, CDouble(..), CInt(..))
import           Foreign.Marshal.Array        (allocaArray)
import           Foreign.Ptr                  (Ptr, nullPtr)

import           System.IO.Unsafe             (unsafePerformIO)

import qualified Data.Vector.Storable.Mutable as SVM

-- | 'CChar' wrapper type for indicating job types
newtype JOBZ = JOBZ CChar deriving (Show, Eq, Storable)

charToJobz :: Char -> JOBZ
charToJobz = JOBZ . castCharToCChar

-- | Pass to 'JOBZ' arguments to only calculate eigenvalues
jzEigvals :: JOBZ
jzEigvals = charToJobz 'N'
-- | Pass to 'JOBZ' arguments to calculate eigenvectors and eigenvalues
jzEigvalsEigvecs :: JOBZ
jzEigvalsEigvecs = charToJobz 'V'

-- | 'CChar' wrapper type for indicating range parameters for calculating
-- eigenvectors/eigenvalues
newtype RANGE = RANGE CChar deriving (Show, Eq, Storable)

charToRange :: Char -> RANGE
charToRange = RANGE . castCharToCChar

-- | Pass to 'RANGE' arguments to calculate for all eigenvalues.
rngAll :: RANGE
rngAll = charToRange 'A'
-- | Pass to 'RANGE' arguments to calculate for all eigenvalues whose value
-- falls in some interval.
rngInterval :: RANGE
rngInterval = charToRange 'V'
-- | Pass to 'RANGE' arguments to calculate for some interval of eigenvalues
-- numbered in ascending order.
rngEigNums :: RANGE
rngEigNums = charToRange 'I'

-- | 'CChar' wrapper type for indicating upper/lower triangular matrix representation.
newtype UPLO = UPLO CChar deriving (Show, Eq, Storable)

charToUplo :: Char -> UPLO
charToUplo = UPLO . castCharToCChar

-- | Pass to 'UPLO' arguments if input matrix is upper triangular.
uploUpper :: UPLO
uploUpper = charToUplo 'U'
-- | Pass to 'UPLO' arguments if input matrix is lower triangular.
uploLower :: UPLO
uploLower = charToUplo 'L'

foreign import ccall unsafe "dlamch_"
        f_dlamch :: Ptr CChar -> CDouble

-- | Value calculated from Fortran function dlamch('S') that is the tolerance
-- that produces the most accurate values for eigenvalue/eigenvector problems
-- by routines in LAPACK.
globalFAbstol :: Double
globalFAbstol = realToFrac . unsafePerformIO $ with (castCharToCChar 'S') $ return . f_dlamch

foreign import ccall unsafe "ilaenv_"
  f_ilaenv :: Ptr CInt
           -> Ptr CChar
           -> Ptr CChar
           -> Ptr CInt
           -> Ptr CInt
           -> Ptr CInt
           -> Ptr CInt
           -> CInt

hsilaenv :: (Num b) =>
            CInt -- ^ ISPEC
         -> String -- ^ routine name
         -> String -- ^ routine options
         -> CInt
         -> CInt
         -> CInt
         -> CInt
         -> b
hsilaenv ispec name opts n1 n2 n3 n4 =
  fromIntegral . unsafePerformIO $
  with ispec $ \pispec ->
  withCAString name $ \pname ->
  withCAString opts $ \popts ->
  with n1 $ \pn1 -> with n2 $ \pn2 -> with n3 $ \pn3 -> with n4 $ \pn4 ->
  return $ f_ilaenv pispec pname popts pn1 pn2 pn3 pn4
  
zheevrBS :: Int
zheevrBS = max (hsilaenv 1 "zhetrd" "" (-1) (-1) (-1) (-1))
               (hsilaenv 1 "zunmtr" "" (-1) (-1) (-1) (-1))

condWith :: Storable a => Bool -> a -> (Ptr a -> IO b) -> IO b
condWith False _ f = f nullPtr
condWith True v f = with v f

condAllocaArray :: Storable a => Bool -> Int -> (Ptr a -> IO b) -> IO b
condAllocaArray False _ f = f nullPtr
condAllocaArray True d f = allocaArray d f

foreign import ccall unsafe "zhbevx_" f_zhbevx
        :: Ptr JOBZ -- JOBZ
        -> Ptr RANGE -- RANGE
        -> Ptr UPLO -- UPLO
        -> Ptr CInt -- N (order of matrix)
        -> Ptr CInt -- KD - Number of super/subdiagonals
        -> Ptr FDComplex -- AB - Upper/lower Hermitian band matrix
        -> Ptr CInt -- LDAB - Leading dimension of array AB
        -> Ptr FDComplex -- Q
        -> Ptr CInt -- LDQ
        -> Ptr CDouble -- VL
        -> Ptr CDouble -- VU
        -> Ptr CInt -- IL
        -> Ptr CInt -- IU
        -> Ptr CDouble -- ABSTOL
        -> Ptr CInt -- M
        -> Ptr CDouble -- W
        -> Ptr FDComplex -- Z
        -> Ptr CInt -- LDZ
        -> Ptr FDComplex -- Work
        -> Ptr CDouble -- RWORK
        -> Ptr CInt -- IWORK
        -> Ptr CInt -- IFAIL
        -> Ptr CInt -- INFO
        -> IO ()

-- | Haskell wrapper for FORTRAN ZHBEVX function to find eigenvalues and
-- eigenvectors of Hermitian band matrices.
hszhbevx :: JOBZ -- ^ the "type" of job, one of:
                 -- * 'jzEigvals' for eigenvalues
                 -- * 'jzEigvalsEigvecs' for eigenvalues and eigenvectors
         -> RANGE -- ^ the "type" of range specification, one of:
                  -- * 'rngAll' for all eigenvalues/eigenvectors
                  -- * 'rngEigNums' for a range of eigenvalues/eigenvectors
                  -- * 'rngInterval' for eigenvalues/eigenvectors whose eigenvalues lie in an interval
         -> UPLO -- ^ the orientation ('uploUpper' or 'uploLower') of the
                 -- input matrix
         -> Int -- ^ the order or dimension of the input matrix
         -> Int -- ^ the number of super/sub-diagonals in the matrix
                -- (generally 1 less than the number of rows)
         -> Ptr FDComplex -- ^ the input matrix (will get destroyed)
         -> Int -- ^ the leading dimension of the input matrix (generally super/sub-diagonals + 1)
         -> Double -- ^ used to determine lower bound for eigenvalue interval
         -> Double -- ^ used to determine uppwer bound for eigenvalue interval
         -> Int -- ^ used to determine lowest eigenvalue index
         -> Int -- ^ used to determine highest eigenvalue index
         -> Double -- ^ the eigenvalue tolerance -- most accurate results are
                   -- achieved by using 'globalFAbstol'
         -> IO (SVM.IOVector CDouble, Maybe (SVM.IOVector FDComplex))
            -- ^ returns a tuple of mutable Data.Vectors, the first
            -- containing the eigenvalues (in ascending order), the second
            -- containing the corresponding eigenvectors in column-major
            -- order (i.e. eigenvectors are stored contiguously)
hszhbevx jz rng ul mo sd pab ldab vl vu il iu tol =
    with jz $ \pjz ->
    with rng $ \prng ->
    with ul $ \pul ->
    with (fromIntegral mo) $ \pmo ->
    with (fromIntegral sd) $ \psd ->
    with (fromIntegral ldab) $ \pldab ->
    condWith (rng == rngInterval) (realToFrac vl) $ \pvl ->
    condWith (rng == rngInterval) (realToFrac vu) $ \pvu ->
    condWith (rng == rngEigNums) (fromIntegral il) $ \pil ->
    condWith (rng == rngEigNums) (fromIntegral iu) $ \piu ->
    with (realToFrac tol) $ \ptol -> do
      w <- SVM.new mo
      z <- if jz == jzEigvalsEigvecs
           then fmap Just (SVM.new (mo * mo))
           else return Nothing
      SVM.unsafeWith w $ \pw ->
          (\f -> maybe (f nullPtr) (`SVM.unsafeWith` f) z) $ \pz ->
          condAllocaArray (jz == jzEigvalsEigvecs) (mo * mo) $ \pq ->
          allocaArray mo $ \pwork ->
          allocaArray (7 * mo) $ \prwork ->
          allocaArray (5 * mo) $ \piwork ->
          allocaArray mo $ \pifail ->
              alloca $ \pm ->
              alloca $ \pinfo -> do
                f_zhbevx pjz
                         prng
                         pul
                         pmo
                         psd
                         pab
                         pldab
                         pq
                         pmo
                         pvl
                         pvu
                         pil
                         piu
                         ptol
                         pm
                         pw
                         pz
                         pmo
                         pwork
                         prwork
                         piwork
                         pifail
                         pinfo
                return (w, z)

foreign import ccall unsafe "zheevr_" f_zheevr
  :: Ptr JOBZ
  -> Ptr RANGE
  -> Ptr UPLO -- ^ whether to use upper/lower triangle of matrix
  -> Ptr CInt -- ^ N (order of matrix)
  -> Ptr FDComplex -- ^ Hermitian matrix A
  -> Ptr CInt -- ^ Leading dimension of array A
  -> Ptr CDouble -- ^ VL
  -> Ptr CDouble -- ^ VU
  -> Ptr CInt -- ^ IL
  -> Ptr CInt -- ^ IU
  -> Ptr CDouble -- ^ ABSTOL
  -> Ptr CInt -- ^ M
  -> Ptr CDouble -- ^ W
  -> Ptr FDComplex -- ^ Z
  -> Ptr CInt -- ^ LDZ
  -> Ptr CInt -- ^ ISUPPZ
  -> Ptr FDComplex -- ^ WORK
  -> Ptr CInt -- ^ LWORK
  -> Ptr CDouble -- ^ RWORK
  -> Ptr CInt -- ^ LRWORK
  -> Ptr CInt -- ^ IWORK
  -> Ptr CInt -- ^ LIWORK
  -> Ptr CInt -- ^ INFO
  -> IO ()

hszheevr :: JOBZ -- ^ the "type" of job
         -> RANGE -- ^ the "type" of range specification
         -> UPLO -- ^ the portion of the input matrix to use
         -> Int -- ^ the order of the input matrix
         -> Ptr FDComplex -- ^ the input matrix (gets destroyed)
         -> Int -- ^ the leading dimension of the input matrix (generally
                -- equal to order)
         -> Double -- ^ lower bound for eigenvalue interval
         -> Double -- ^ upper bound for eigenvalue interval
         -> Int -- ^ lowest eigenvalue index
         -> Int -- ^ highest eigenvalue index
         -> Double -- ^ eigenvalue tolerance -- most accurate results
                   -- achieved by using 'globalFAbstol'
         -> IO (SVM.IOVector CDouble, Maybe (SVM.IOVector FDComplex))
         -- ^ returns a tuple of mutable vectors, the first containing
         -- eigenvalues in ascending order, the second containing the
         -- corresponding eigenvalues in column-major order
hszheevr jz rng ul mo pmat ldm vl vu il iu tol =
  with jz $ \pjz ->
  with rng $ \prng ->
  with ul $ \pul ->
  with (fromIntegral mo) $ \pmo ->
  with (fromIntegral ldm) $ \pldm ->
  condWith (rng == rngInterval) (realToFrac vl) $ \pvl ->
  condWith (rng == rngInterval) (realToFrac vu) $ \pvu ->
  condWith (rng == rngEigNums) (fromIntegral il) $ \pil ->
  condWith (rng == rngEigNums) (fromIntegral iu) $ \piu ->
  with (realToFrac tol) $ \ptol ->
  alloca $ \pm ->
  allocaArray (2 * mo) $ \pisupz ->
  allocaArray lwork $ \pwork ->
  with (fromIntegral lwork) $ \plwork ->
  allocaArray lrwork $ \prwork ->
  with (fromIntegral lrwork) $ \plrwork ->
  allocaArray liwork $ \piwork ->
  with (fromIntegral liwork) $ \pliwork ->
  alloca $ \pinfo -> do
    w <- SVM.new mo
    z <- if jz == jzEigvalsEigvecs
            then fmap Just (SVM.new ldz)
            else return Nothing
    SVM.unsafeWith w $ \pw ->
      (\f -> maybe (f nullPtr) (`SVM.unsafeWith` f) z) $ \pz -> do
        f_zheevr pjz
                 prng
                 pul
                 pmo
                 pmat
                 pldm
                 pvl
                 pvu
                 pil
                 piu
                 ptol
                 pm
                 pw
                 pz
                 pmo
                 pisupz
                 pwork
                 plwork
                 prwork
                 plrwork
                 piwork
                 pliwork
                 pinfo
        return (w, z) 
  where lwork = (zheevrBS + 1) * mo
        lrwork = 24 * mo
        liwork = 10 * mo
        ldz = mo * mo
