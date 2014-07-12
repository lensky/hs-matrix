{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}

module Data.Matrix
       (
       -- * Matrix classes
       Matrix(..)
       ,MatrixRing(..)
       -- * Utility methods
       ,showAsMatrix
       ,putAsMatrix)
where

import Data.Tensor hiding (generate,generateM)
import qualified Data.Tensor as T
import Data.Tensor.MultiVector (MultiVector(..), cacheCoeffs, deepReIndex)
import qualified Data.Permutation as P

import qualified Data.List as LST

import qualified Data.Vector.Generic as GV

-- | Class for an object that can be seen as a matrix.
class Matrix t s a where
  rows :: t s a -> Int
  cols :: t s a -> Int

  row :: t s a -> Int -> s a
  col :: t s a -> Int -> s a

  -- | Row-column indexing.
  ij :: t s a -> Int -> Int -> a
  transpose :: t s a -> t s a

  -- | Generate a matrix from a function that takes a row and column as input.
  generate :: Int -> Int -> (Int -> Int -> a) -> t s a
  -- | Generate a matrix from a monadic function that takes a row and column as input.
  generateM :: (Monad m) => Int -> Int -> (Int -> Int -> m a) -> m (t s a)

  -- | Generate the matrix from a list. Generally not the most efficient way
  -- to create a matrix, and mostly useful for testing.
  fromList :: [[a]] -> t s a
  
  toColList :: t s a -> [s a]
  toColList m = map (col m) [0..cols m - 1]
  
  toRowList :: t s a -> [s a]
  toRowList m = map (row m) [0..rows m - 1]

-- | Class to implement matrix addition.
class (Matrix m s a, Matrix m' s' a') => MatrixRing m s a m' s' a' where
  -- | Type for the result of addition. This is important, for example, when
  -- adding banded to sparse matrices, when the result is not necessarily of
  -- either type.
  type MSumMatrix m s a m' s' a' :: (* -> *) -> * -> *
  type MSumElement m s a m' s' a' :: *
  type MSumStore m s a m' s' a' :: * -> *

  -- | Matrix element-wise addition.
  mp :: (Matrix (MSumMatrix m s a m' s' a') (MSumStore m s a m' s' a') (MSumElement m s a m' s' a')) => 
        m s a -> m' s' a' -> (MSumMatrix m s a m' s' a') (MSumStore m s a m' s' a') (MSumElement m s a m' s' a') 

-- | Utility function to produce a string representation of a matrix more
-- suitable for human consumption. Mainly useful for debugging.
showAsMatrix :: (Show a, Matrix t s a) => t s a -> String
showAsMatrix m = LST.intercalate "\n" . map showRow $ [0..(rs - 1)]
  where rs = rows m
        cs = cols m
        showRow r = unwords . map (show . ij m r) $ [0..(cs - 1)]

-- | Utility function to show a representation of a matrix more suitable for
-- human consumption. Mainly useful for debugging.
putAsMatrix :: (Show a, Matrix t s a) => t s a -> IO ()
putAsMatrix = putStrLn . showAsMatrix

instance (GV.Vector v a, GV.Vector v Int) => Matrix (MultiVector MatrixTS) v a where 
  rows = fst . mtsAsPair . shape
  cols = snd . mtsAsPair . shape

  row mv r = mvData (mv ?? pairAsMSlice (Single r, All))
  col mv c = mvData (mv ?? pairAsMSlice (All, Single c))

  ij mv r c = mv ? listAsSel [r,c]

  transpose = reIndex (P.fromList [1,0])

  generate rs cs f = T.generate (pairAsMts (rs,cs)) (uncurry f . mSelAsPair)

  generateM rs cs f = T.generateM (pairAsMts (rs,cs)) (uncurry f . mSelAsPair)

  fromList l = cacheCoeffs MultiVector { mvShape = pairAsMts (rs,cs)
                                       , mvLinearStorageOrder = P.identity 2
                                       , mvData = GV.fromList $ LST.concat l}
    where rs = length l
          cs = length $ head l

instance (Num a, GV.Vector v a, GV.Vector v Int) => 
  MatrixRing (MultiVector MatrixTS) v a (MultiVector MatrixTS) v a where 
  type MSumMatrix (MultiVector MatrixTS) v a (MultiVector MatrixTS) v a = (MultiVector MatrixTS)
  type MSumElement (MultiVector MatrixTS) v a (MultiVector MatrixTS) v a = a
  type MSumStore (MultiVector MatrixTS) v a (MultiVector MatrixTS) v a = v
  mp m m' | mvLinearStorageOrder m == mvLinearStorageOrder m' = 
            m { mvData = GV.zipWith (+) (mvData m) (mvData m')}
          | otherwise = mp m (deepReIndex (P.fromList [1,0]) m')
