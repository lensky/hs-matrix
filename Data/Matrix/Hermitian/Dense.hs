{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

module Data.Matrix.Hermitian.Dense
       (DenseHMatrix(..)
       ,fullEigensystem) 
where

import Data.Matrix
import Data.Matrix.Dense
import LAPACK
import LAPACK.Utils
import Data.Complex
import Data.Complex.Utils
import Data.Maybe (fromMaybe)
import Operator

import qualified Data.List as LST

import Control.Monad (forM_)
import Control.Monad.ST (runST)

import System.IO.Unsafe (unsafePerformIO)

import Foreign.C.Types (CDouble)

import qualified Data.Vector.Generic as GV
import qualified Data.Vector.Generic.Mutable as GVM
import qualified Data.Vector.Storable as SV
import qualified Data.Vector.Storable.Mutable as SVM
import qualified Data.Vector.Unboxed as UV

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
  
instance (Conjugable a, GV.Vector v a) => MatrixRing DenseHMatrix v a DenseHMatrix v a where 
  type MSumMatrix DenseHMatrix v a DenseHMatrix v a = DenseHMatrix
  type MSumStore DenseHMatrix v a DenseHMatrix v a = v
  type MSumElement DenseHMatrix v a DenseHMatrix v a = a

  mp m@(DenseHMatrix _ f d) (DenseHMatrix _ f' d')
    = m { dhmFull = f && f'
        , dhmData = GV.zipWith (+) d d'
        }

fullEigensystem :: (FComplexable a CDouble 
                   , GV.Vector v a
                   , GV.Vector v (FComplex CDouble)
                   , GV.Vector v' Double
                   , GV.Vector v' CDouble
                   , GV.Vector v'' (FComplex CDouble)
                   , GV.Vector v'' (Complex Double)) =>
                   DenseHMatrix v a -- ^ input matrix
                   -> Bool -- ^ calculate eigenvectors
                   -> RANGE -- ^ eigenvalue calculation range type
                   -> Double -- ^ lower bound of eigenvalue interval
                   -> Double -- ^ upper bound of eigenvalue interval
                   -> Int -- ^ lower eigenvalue index
                   -> Int -- ^ upper eigenvalue index
                   -> (v' Double, Maybe (v'' (Complex Double)))
fullEigensystem (DenseHMatrix mo _ d) vecs rng vl vu il iu =
  unsafePerformIO $ do
    mdat <- SV.unsafeThaw . GV.convert . GV.map toFComplex $ d
    SVM.unsafeWith mdat $ \pmat -> 
      hszheevr jz rng uploUpper mo pmat mo vl vu (il + 1) (iu + 1) globalFAbstol >>= freezeHermitianES
  where jz = if vecs
                then jzEigvalsEigvecs
                else jzEigvals

instance (FComplexable a CDouble
         ,GV.Vector v a
         ,GV.Vector v (FComplex CDouble)) =>
  EndoOperator (DenseHMatrix v a) where
  type Eigenvalue (DenseHMatrix v a) = Double
  type EigenvalueStorage (DenseHMatrix v a) = UV.Vector
  type EigenvectorStorage (DenseHMatrix v a) = DenseMatrix UV.Vector (Complex Double)

  eigvals m (Just (lo,hi)) = GV.unsafeTake (hi - lo + 1) . fst $
                             (fullEigensystem m False rngEigNums 0 0 lo hi
                               :: (UV.Vector Double, Maybe (UV.Vector (Complex Double))))
  eigvals m Nothing = fst
                      (fullEigensystem m False rngAll 0 0 0 0
                        :: (UV.Vector Double, Maybe (UV.Vector (Complex Double))))

  eigvecs m (Just (lo,hi)) =
    DenseMatrix { dmRows = rs
                , dmCols = cs
                , dmRep = ColMajor
                , dmData = maybe GV.empty (GV.unsafeTake (rs*cs)) . snd $
                           (fullEigensystem m True rngEigNums 0 0 lo hi
                             :: (UV.Vector Double, Maybe (UV.Vector (Complex Double)))) }
    where rs = dhmOrder m
          cs = hi - lo + 1
  eigvecs m Nothing =
    DenseMatrix { dmRows = dhmOrder m
                , dmCols = dhmOrder m
                , dmRep = ColMajor
                , dmData = fromMaybe GV.empty . snd $ 
                           (fullEigensystem m True rngAll 0 0 0 0
                             :: (UV.Vector Double, Maybe (UV.Vector (Complex Double)))) }

  adjoint = id
