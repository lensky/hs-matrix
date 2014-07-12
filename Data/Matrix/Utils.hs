{-# LANGUAGE FlexibleContexts #-}

module Data.Matrix.Utils 
       (switchMajorAxis)
where

import qualified Data.Vector.Generic as GV

-- | Exchanges the major axes of a 2D vector.
switchMajorAxis :: (GV.Vector v a, GV.Vector v Int) => v a -> Int -> Int -> v a
switchMajorAxis v majm subm =
  GV.unsafeBackpermute v $ GV.generate (majm * subm) $ \n ->
    let (mn, sn) = quotRem n majm
    in sn * subm + mn
