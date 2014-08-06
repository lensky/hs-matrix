{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}

module Data.Matrix.Dense
       (DenseMatrix
       ,fromColMajorVector)
 where
 
import Data.Tensor hiding (generate,generateM)
import Data.Tensor.Dense.VTensor (VTensor(..),MD_VTensor)
import qualified Data.Permutation as P

import qualified Data.Vector.Generic as GV

-- | Synonym for a matrix-like multivector.
type DenseMatrix v a = MD_VTensor LinearStorageC v a

-- | Create a dense matrix from column-major raw vector data.
fromColMajorVector :: (GV.Vector v a) => Int -> Int -> v a -> DenseMatrix v a
fromColMajorVector rs cs v = 
  VTensor { _vStorageScheme = fromShapeOrder (rs,cs) (P.fromList [1,0]) 
          , _vData = v }
