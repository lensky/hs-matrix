{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}

module Data.Matrix.Dense
       (DenseMatrix(..)
       ,vecChangeRep)
 where
 
import Data.Matrix
import Data.Matrix.Utils

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

  transpose m@(DenseMatrix rs cs RowMajor _) = m { dmRows = cs
                                                 , dmCols = rs
                                                 , dmRep = ColMajor }
  transpose m@(DenseMatrix rs cs ColMajor _) = m { dmRows = cs
                                                 , dmCols = rs
                                                 , dmRep = RowMajor }  

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
    
vecChangeRep :: (GV.Vector v Int, GV.Vector v a) => DenseMatrix v a -> DenseMatrix v a
vecChangeRep m@(DenseMatrix rs cs rep d)
  = m { dmRep = newRep, dmData = d' }
  where newRep = case rep of
                   ColMajor -> RowMajor
                   RowMajor -> ColMajor
        (majm, subm) = case rep of
                         ColMajor -> (cs, rs)
                         RowMajor -> (rs, cs)
        d' = switchMajorAxis d majm subm

instance (Functor s) => Functor (DenseMatrix s) where
  fmap f m@(DenseMatrix _ _ _ d) = m { dmData = fmap f d }

instance (Num a, GV.Vector v a, GV.Vector v Int) => 
  MatrixRing DenseMatrix v a DenseMatrix v a where 
  type MSumMatrix DenseMatrix v a DenseMatrix v a = DenseMatrix
  type MSumElement DenseMatrix v a DenseMatrix v a = a
  type MSumStore DenseMatrix v a DenseMatrix v a = v
  mp m@(DenseMatrix rs cs rep d) m'@(DenseMatrix _ _ rep' d')
     | rep == rep' = DenseMatrix { dmRows = rs
                                 , dmCols = cs
                                 , dmRep = rep
                                 , dmData = GV.zipWith (+) d d' }
     | otherwise = mp m (vecChangeRep m')
