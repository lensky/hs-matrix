{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module DenseMatrix
       (DenseMatrix(..))
 where

import Matrix

import Data.Tuple (swap)

import qualified Data.List as LST

import qualified Data.Vector.Generic as GV

data DenseMatrix s a
  = DenseMatrix
    { dmRows :: Int
    , dmCols :: Int
    , dmRep :: ArrayRep
    , dmData :: s a 
    }

instance (GV.Vector v a) => Matrix DenseMatrix v a where 
  rows = dmRows
  cols = dmCols

  row (DenseMatrix _ cs RowMajor d) r = GV.unsafeSlice (r * cs) cs d
  row m@(DenseMatrix rs _ ColMajor _) r = GV.generate rs (ij m r)

  col (DenseMatrix rs _ ColMajor d) c = GV.unsafeSlice (c * rs) rs d
  col m@(DenseMatrix _ cs RowMajor _) c = GV.generate cs (flip (ij m) c)

  ij (DenseMatrix _ cs RowMajor d) r c = d `GV.unsafeIndex` (r * cs + c)
  ij (DenseMatrix rs _ ColMajor d) r c = d `GV.unsafeIndex` (c * rs + r)

  transpose m@(DenseMatrix _ _ RowMajor _) = m { dmRep = ColMajor }
  transpose m@(DenseMatrix _ _ ColMajor _) = m { dmRep = RowMajor }  

  generate rs cs f =
    DenseMatrix { dmRows = rs
                , dmCols = cs
                , dmData = GV.generate (rs * cs) (uncurry f . rc)
                , dmRep = ColMajor }
    where rc = swap . flip quotRem rs

  generateM rs cs f = do
    d <- GV.generateM (rs * cs) (uncurry f . rc)
    return DenseMatrix { dmRows = rs, dmCols = cs, dmRep = ColMajor, dmData = d }
    where rc = swap . flip quotRem rs

  fromList l = DenseMatrix
    { dmRows = length l
    , dmCols = length $ head l
    , dmRep = RowMajor
    , dmData = GV.fromList $ LST.concat l}
    
instance (Functor s) => Functor (DenseMatrix s) where
  fmap f m@(DenseMatrix _ _ _ d) = m { dmData = fmap f d }
