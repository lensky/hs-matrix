{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}

module Data.Matrix
       (Matrix(..)
       ,showAsMatrix
       ,putAsMatrix
       ,ArrayRep(..)
       ,MatrixRing(..))
 where

import qualified Data.List as LST

class Matrix t s a where
  rows :: t s a -> Int
  cols :: t s a -> Int

  row :: t s a -> Int -> s a
  col :: t s a -> Int -> s a

  ij :: t s a -> Int -> Int -> a
  transpose :: t s a -> t s a

  generate :: Int -> Int -> (Int -> Int -> a) -> t s a
  generateM :: (Monad m) => Int -> Int -> (Int -> Int -> m a) -> m (t s a)

  fromList :: [[a]] -> t s a
  
  toColList :: t s a -> [s a]
  toColList m = map (col m) [0..cols m - 1]
  
  toRowList :: t s a -> [s a]
  toRowList m = map (row m) [0..rows m - 1]

class (Matrix m s a, Matrix m' s' a') => MatrixRing m s a m' s' a' where
  type MSumMatrix m s a m' s' a' :: (* -> *) -> * -> *
  type MSumElement m s a m' s' a' :: *
  type MSumStore m s a m' s' a' :: * -> *

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

data ArrayRep = RowMajor | ColMajor deriving (Eq, Show)
