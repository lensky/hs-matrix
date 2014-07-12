{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

module Data.Matrix.Hermitian.Banded
       (BandedHMatrix(..)
       ,vecChangeRep
       ,TriangularRepr
       ,upperRep
       ,lowerRep
       ,fullEigensystem)
where

import Data.Matrix
import Data.Matrix.Utils
import Data.Matrix.Dense (DenseMatrix,fromColMajorVector)
import LAPACK
import LAPACK.Utils
import Data.Complex
import Data.Complex.Utils
import Operator

import qualified Data.List as LST

import Data.Maybe (fromMaybe)

import Control.Monad (forM_)
import System.IO.Unsafe (unsafePerformIO)

import Foreign.C.Types (CDouble)

import qualified Data.Vector.Generic as GV
import qualified Data.Vector.Generic.Mutable as GVM
import qualified Data.Vector.Storable as SV
import qualified Data.Vector.Storable.Mutable as SVM
import qualified Data.Vector.Unboxed as UV

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

-- | Change the in-memory vector representation for a 'BandedHMatrix'.                 
vecChangeRep :: (Conjugable a, GV.Vector v a) => BandedHMatrix v a -> BandedHMatrix v a
vecChangeRep m@(BHMatrix mo sd rep d)
  = m { bhmRep = newRep, bhmData = d' }
  where newRep = if rep == upperRep
                    then lowerRep
                    else upperRep
        swapsU i = let i' = i + sd + (sd - 1) * rem i (sd + 1)
                       l = GV.length d
                   in if i' >= l
                         then (i,i)
                         else (i,i')
        swapsL i = let (j,j') = swapsU i
                       l' = GV.length d - 1
                   in (l' - j, l' - j')
        swaps = if rep == upperRep
                   then swapsU
                   else swapsL 
        d' = GV.modify (\v -> 
               forM_ [0 .. ((sd + 1) * mo) - 1 - sd] $ \n ->
                 uncurry (GVM.swap v) $ swaps n)
              (GV.map cconj d)

instance (Conjugable a, GV.Vector v a, GV.Vector v Int) =>
  MatrixRing BandedHMatrix v a BandedHMatrix v a where
  type MSumMatrix BandedHMatrix v a BandedHMatrix v a = BandedHMatrix
  type MSumElement BandedHMatrix v a BandedHMatrix v a = a
  type MSumStore BandedHMatrix v a BandedHMatrix v a = v

  mp m@(BHMatrix mo sd rep d) m'@(BHMatrix _ sd' rep' d')
    | (rep == rep) && (sd == sd') = m { bhmData = GV.zipWith (+) d d'}
    | rep == rep' = let si = if rep == upperRep
                                then (sd - sd') * mo
                                else (sd' + 1) * mo
                        d'' = switchMajorAxis d' mo (sd' + 1)
                        (ds, de) = GV.splitAt si (switchMajorAxis d mo (sd + 1))
                        df = if rep == upperRep
                                then ds GV.++ GV.zipWith (+) de d''
                                else GV.zipWith (+) ds d' GV.++ de
                    in if sd >= sd'
                          then m { bhmData = switchMajorAxis df (sd + 1) mo }
                          else mp m' m
    | otherwise = mp m (vecChangeRep m')

-- | Wrapper to "raw" Haskell function 'hszhbevx' for the eigensystem of a
-- banded Hermitian matrix.
fullEigensystem :: (FComplexable a CDouble 
                   , GV.Vector v a
                   , GV.Vector v (FComplex CDouble)
                   , GV.Vector v' Double
                   , GV.Vector v' CDouble
                   , GV.Vector v'' (FComplex CDouble)
                   , GV.Vector v'' (Complex Double)) =>
                   BandedHMatrix v a -- ^ the input matrix
                   -> Bool -- ^ calculate eigenvectors?
                   -> RANGE -- ^ eigenvalue calculation range type
                   -> Double -- ^ used as lower bound for eigenvalue interval
                   -> Double -- ^ used as upper bound for eigenvalue interval
                   -> Int -- ^ used as lower eigenvalue number
                   -> Int -- ^ used as upper eigenvalue number
                   -> (v' Double, Maybe (v'' (Complex Double)))
fullEigensystem bhm vecs rng vl vu il iu =
  unsafePerformIO $ do
    mdat <- SV.unsafeThaw . GV.convert . GV.map toFComplex . bhmData $ bhm
    SVM.unsafeWith mdat $ \pab ->
        hszhbevx 
          jz rng ul
          mo sd pab (sd + 1)
          vl vu (il + 1) (iu + 1) globalFAbstol >>= freezeHermitianES
  where jz = if vecs then jzEigvalsEigvecs else jzEigvals 
        ul = if bhmRep bhm == upperRep
                then uploUpper
                else uploLower
        mo = bhmOrder bhm
        sd = bhmSDs bhm

instance (FComplexable a CDouble
         ,GV.Vector v a
         ,GV.Vector v (FComplex CDouble)) =>
  EndoOperator (BandedHMatrix v a) where
  type Eigenvalue (BandedHMatrix v a) = Double
  type EigenvalueStorage (BandedHMatrix v a) = UV.Vector
  type EigenvectorStorage (BandedHMatrix v a) = DenseMatrix UV.Vector (Complex Double)

  eigvals m (Just (lo,hi)) = GV.unsafeTake (hi - lo + 1) . fst $ 
                             (fullEigensystem m False rngEigNums 0 0 lo hi 
                               :: (UV.Vector Double, Maybe (UV.Vector (Complex Double))))
  eigvals m Nothing = fst
    (fullEigensystem m False rngAll 0 0 0 0 
    :: (UV.Vector Double, Maybe (UV.Vector (Complex Double))))

  eigvecs m (Just (lo,hi)) = fromColMajorVector rs cs . 
                             maybe GV.empty (GV.unsafeTake (rs * cs)) $ 
                             snd (fullEigensystem m True rngEigNums 0 0 lo hi 
                               :: (UV.Vector Double, Maybe (UV.Vector (Complex Double))))
    where rs = bhmOrder m
          cs = hi - lo + 1
          
  eigvecs m Nothing = fromColMajorVector (bhmOrder m) (bhmOrder m) . 
                      fromMaybe GV.empty . snd $
                      (fullEigensystem m True rngAll 0 0 0 0
                        :: (UV.Vector Double, Maybe (UV.Vector (Complex Double))))

  adjoint = id

instance (Functor s) => Functor (BandedHMatrix s) where
  fmap f m@(BHMatrix _ _ _ d) = m { bhmData = fmap f d }
