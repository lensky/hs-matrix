{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}

module Data.Matrix.Dense
       (DenseMatrix
       ,fromColMajorVector)
 where
 
import Data.Tensor hiding (generate,generateM)
import Data.Tensor.MultiVector (MatrixLike, MultiVector(..), cacheCoeffs)
import qualified Data.Permutation as P

import qualified Data.Vector.Generic as GV

type DenseMatrix v a = MatrixLike v a

-- | Create a dense matrix from column-major raw vector data.
fromColMajorVector :: (GV.Vector v a) => Int -> Int -> v a -> DenseMatrix v a
fromColMajorVector rs cs v = cacheCoeffs MultiVector { mvShape = pairAsMts (rs,cs)
                                                     , mvLinearStorageOrder = P.fromList [1,0]
                                                     , mvData = v }
