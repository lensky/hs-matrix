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


import qualified Data.List as LST
import qualified Data.Permutation as P

import qualified Data.Tensor as T
import Data.Tensor hiding (generate,generateM)
import Data.Tensor.Dense.VTensor (MD_VTensor, VTensor(..))

import qualified Data.Tensor.Dense.VTensor as VT

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

instance (GV.Vector v a, GV.Vector v Int
          ,LinearStorageInfo li (Int,Int) (Int,Int)
          ,P.Permutable (li (Int,Int) (Int,Int)))
         => Matrix (MD_VTensor li) v a where 
  rows = fst . shape
  cols = snd . shape

  row mv r = _vData (mv *? [Inx r, All])
  col mv c = _vData (mv *? [All, Inx c])

  ij mv r c = mv *! (r,c)

  transpose = reIndex (P.fromList [1,0])

  generate rs cs f = T.generate (rs,cs) (uncurry f)

  generateM rs cs f = T.generateM (rs,cs) (uncurry f)

  fromList l = VTensor { _vStorageScheme = lss
                       , _vData = GV.fromList $ LST.concat l }
    where rs = length l
          cs = length $ head l
          lss = fromShapeOrder (rs,cs) (P.identity 2)

instance (Num a
         ,LinearStorageInfo li (Int,Int) (Int,Int)
         ,P.Permutable (li (Int, Int) (Int,Int))
         ,GV.Vector v a, GV.Vector v Int) 
         => MatrixRing (MD_VTensor li) v a (MD_VTensor li) v a where 
  type MSumMatrix (MD_VTensor li) v a (MD_VTensor li) v a = (MD_VTensor li)
  type MSumElement (MD_VTensor li) v a (MD_VTensor li) v a = a
  type MSumStore (MD_VTensor li) v a (MD_VTensor li) v a = v
  mp = VT.zipWith (+)
