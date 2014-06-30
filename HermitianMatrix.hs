{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

module HermitianMatrix
       (DenseHMatrix(..)
       ,BandedHMatrix(..)
       ,TriangularRepr
       ,upperRep
       ,lowerRep)
where

import Matrix
import DenseMatrix
import LAPACK
import ComplexUtils
import Operator

import Data.Maybe (fromMaybe)

import qualified Data.List as LST
import qualified Data.Traversable as T

import Control.Monad (forM_)
import Control.Monad.ST (runST)
import Control.Monad.Primitive (PrimMonad,PrimState)

import System.IO.Unsafe (unsafePerformIO)

import Foreign.C.Types (CDouble)

import qualified Data.Vector.Generic as GV
import qualified Data.Vector.Generic.Mutable as GVM
import qualified Data.Vector.Storable as SV
import qualified Data.Vector.Storable.Mutable as SVM

data DenseHMatrix s a
  = DenseHMatrix
    { dhmOrder :: Int
    , dhmFull :: Bool
    , dhmData :: s a }
    
instance (Conjugable a, GV.Vector v a) => Matrix DenseHMatrix v a where
  rows = dhmOrder
  cols = rows

  row m r = GV.generate (rows m) (ij m r)

  col m@(DenseHMatrix mo f d) c =
    if f
       then GV.unsafeSlice (c * mo) mo d
       else GV.generate mo (flip (ij m) c)

  ij (DenseHMatrix mo True d) r c = d GV.! ((c * mo) + r)
  ij (DenseHMatrix mo False d) r c
    | r <= c = d GV.! ((c * mo) + r)
    | otherwise = cconj $ d GV.! ((r * mo) + c)

  transpose m@(DenseHMatrix _ _ d) = m { dhmData = GV.map cconj d }

  generate mo _ f =
    DenseHMatrix { dhmOrder = mo, dhmData = d, dhmFull = False }
    where d = runST $ do
            mv <- GVM.new (mo * mo)
            forM_ [0..(mo - 1)] $ \c ->
              forM_ [0..c] $ \r ->
                GVM.write mv ((c * mo) + r) (f r c)
            GV.unsafeFreeze mv

  generateM mo _ f = do
    d <- GV.generateM (mo * mo) (f' . flip quotRem mo)
    return DenseHMatrix 
           { dhmOrder = mo, dhmData = d, dhmFull = False }
    where fill = f 0 0
          f' (c,r) | c <= r = f r c
                   | otherwise = fill

  fromList l = DenseHMatrix
    { dhmOrder = length $ head l
    , dhmData = GV.fromList . LST.concat . LST.transpose $ l
    , dhmFull = True }
    
instance (Functor s) => Functor (DenseHMatrix s) where
  fmap f m@(DenseHMatrix _ _ d) = m { dhmData = fmap f d }

freezeLAPACKES
  :: (Fractional b, Real a1, SVM.Storable a1, SVM.Storable a,
      SVM.Storable b, PrimMonad m) =>
     (SVM.MVector (PrimState m) a1, Maybe (SVM.MVector (PrimState m) a))
     -> m (SV.Vector b, SV.Vector a)
freezeLAPACKES (mvals, maybeMVecs) = do
  vals <- SV.unsafeFreeze mvals
  maybeVecs <- T.mapM SV.unsafeFreeze maybeMVecs
  let evecs = fromMaybe GV.empty maybeVecs
  return (SV.map realToFrac vals, evecs)

dhmFullEigensystem :: (FComplexable a CDouble
                      , GV.Vector v a
                      , GV.Vector v FDComplex) =>
                      DenseHMatrix v a -- ^ input matrix
                   -> Bool -- ^ calculate eigenvectors
                   -> RANGE -- ^ eigenvalue calculation range type
                   -> Double -- ^ lower bound of eigenvalue interval
                   -> Double -- ^ upper bound of eigenvalue interval
                   -> Int -- ^ lower eigenvalue index
                   -> Int -- ^ upper eigenvalue index
                   -> (SV.Vector Double, SV.Vector FDComplex)
dhmFullEigensystem (DenseHMatrix mo _ d) vecs rng vl vu il iu =
  unsafePerformIO $ do
    mdat <- SV.unsafeThaw . GV.convert . GV.map toFComplex $ d
    SVM.unsafeWith mdat $ \pmat -> 
      hszheevr jz rng uploUpper mo pmat mo vl vu il iu globalFAbstol >>= freezeLAPACKES
  where jz = if vecs
                then jzEigvalsEigvecs
                else jzEigvals

instance (FComplexable a CDouble, GV.Vector v a, GV.Vector v FDComplex) =>
  EndoOperator (DenseHMatrix v a) where
  type Eigenvalue (DenseHMatrix v a) = Double
  type EigenvalueStorage (DenseHMatrix v a) = SV.Vector
  type EigenvectorStorage (DenseHMatrix v a) = DenseMatrix SV.Vector FDComplex

  eigvals m (Just (lo,hi)) = GV.unsafeSlice 0 (hi - lo + 1) . fst $
    dhmFullEigensystem m False rngEigNums 0 0 lo hi
  eigvals m Nothing = fst $ dhmFullEigensystem m False rngAll 0 0 0 0

  eigvecs m (Just (lo,hi)) =
    DenseMatrix { dmRows = dhmOrder m
                , dmCols = hi - lo + 1
                , dmRep = ColMajor
                , dmData = snd $ dhmFullEigensystem m True rngEigNums 0 0 lo hi }
  eigvecs m Nothing =
    DenseMatrix { dmRows = dhmOrder m
                , dmCols = dhmOrder m
                , dmRep = ColMajor
                , dmData = snd $ dhmFullEigensystem m True rngAll 0 0 0 0 }

  adjoint = id

-- | Type for representations of triangular matrices, either upper or lower.
newtype TriangularRepr = TriRepr Bool deriving (Eq,Show)

-- | Upper triangular representation indicator
upperRep :: TriangularRepr
upperRep = TriRepr True

-- | Lower triangular representation indicator
lowerRep :: TriangularRepr
lowerRep = TriRepr False

-- | Basic FORTRAN-compatible type (i.e. col-major storage) type for banded
-- Hermitian matrices.
data BandedHMatrix v a
  = BHMatrix
    { bhmOrder :: !Int -- ^ order of the (square) matrix
    , bhmSDs :: !Int -- ^ super/sub-diagonals of the triangular matrix
    , bhmRep :: !TriangularRepr -- ^  array storage convention
    , bhmData :: v a -- ^ (column-major) data
    }
    deriving (Eq, Show)

instance (Conjugable a, GV.Vector v a) => Matrix BandedHMatrix v a where
  rows = bhmOrder
  cols = rows

  row m r = GV.generate (bhmOrder m) (ij m r)
  col m c = GV.generate (bhmOrder m) (flip (ij m) c)

  ij (BHMatrix _ s rep d) r c | abs (r - c) > s = 0
                              | inRep = d `GV.unsafeIndex` vi r c
                              | otherwise = cconj $ d `GV.unsafeIndex` vi c r
    where inRep = if rep == upperRep 
                     then c >= r
                     else c <= r
          vi = if rep == upperRep
                  then \i j -> s * (j + 1) + i
                  else \i j -> s * j + i

  transpose m@(BHMatrix _ _ _ d) = m { bhmData = GV.map cconj d }

  fromList bs = BHMatrix { bhmOrder = length . head $ bs
                              , bhmSDs = length bs - 1
                              , bhmRep = upperRep
                              , bhmData = GV.fromList . LST.concat . LST.transpose $ bs }

  generate mo sd f =
    BHMatrix
     { bhmOrder = mo
     , bhmSDs = sd
     , bhmRep = upperRep
     , bhmData = GV.generate ((sd + 1) * mo) (uncurry f . rc) } 
    where rc i = let (c,r) = quotRem i (sd + 1)
                 in (c - sd + r, c)

  generateM mo sd f = do
    d <- GV.generateM ((sd + 1) * mo) (uncurry f . rc)
    return BHMatrix { bhmOrder = mo
                    , bhmSDs = sd
                    , bhmRep = upperRep
                    , bhmData = d }
    where rc i = let (c,r) = quotRem i (sd + 1)
                 in (c - sd + r, c)

-- | Wrapper to "raw" haskell function that calculates a full eigensystem for
-- banded hermitian matrices.
bhmFullEigensystem :: (FComplexable a CDouble,
                       GV.Vector v a, 
                       GV.Vector v FDComplex) =>
                      BandedHMatrix v a -- ^ the input matrix
                   -> Bool -- ^ calculate eigenvectors?
                   -> RANGE -- ^ eigenvalue calculation range type
                   -> Double -- ^ used as lower bound for eigenvalue interval
                   -> Double -- ^ used as upper bound for eigenvalue interval
                   -> Int -- ^ used as lower eigenvalue number
                   -> Int -- ^ used as upper eigenvalue number
                   -> (SV.Vector Double, SV.Vector FDComplex)
bhmFullEigensystem bhm vecs rng vl vu il iu =
  unsafePerformIO $ do
    mdat <- SV.unsafeThaw . GV.convert . GV.map toFComplex . bhmData $ bhm
    SVM.unsafeWith mdat $ \pab ->
        hszhbevx 
          jz rng ul
          mo sd pab (sd + 1)
          vl vu il iu globalFAbstol >>= freezeLAPACKES
  where jz = if vecs then jzEigvalsEigvecs else jzEigvals 
        ul = if bhmRep bhm == upperRep
                then uploUpper
                else uploLower
        mo = bhmOrder bhm
        sd = bhmSDs bhm

instance (FComplexable a CDouble, GV.Vector v a, GV.Vector v FDComplex) =>
  EndoOperator (BandedHMatrix v a) where
  type Eigenvalue (BandedHMatrix v a) = Double
  type EigenvalueStorage (BandedHMatrix v a) = SV.Vector
  type EigenvectorStorage (BandedHMatrix v a) = DenseMatrix SV.Vector FDComplex

  eigvals m (Just (lo,hi)) = GV.unsafeSlice 0 (hi - lo + 1) $ 
                             fst $ bhmFullEigensystem m False rngEigNums 0 0 lo hi
  eigvals m Nothing = fst $ bhmFullEigensystem m False rngAll 0 0 0 0

  eigvecs m (Just (lo,hi)) =
    DenseMatrix { dmRows = bhmOrder m
                , dmCols = hi - lo + 1
                , dmRep = ColMajor
                , dmData = snd $ bhmFullEigensystem m True rngEigNums 0 0 lo hi }
  eigvecs m Nothing =
    DenseMatrix { dmRows = bhmOrder m
                , dmCols = bhmOrder m
                , dmRep = RowMajor
                , dmData = snd $ bhmFullEigensystem m True rngAll 0 0 0 0 }

  adjoint = id

instance (Functor s) => Functor (BandedHMatrix s) where
  fmap f m@(BHMatrix _ _ _ d) = m { bhmData = fmap f d }
