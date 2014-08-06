{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE MultiParamTypeClasses #-}

{-# OPTIONS_GHC -fno-warn-orphans #-}

module LAPACK.Utils
       (freezeConvert
       ,freezeHermitianES)
where

import Control.Applicative

import Control.Monad (liftM)
import Control.Monad.Primitive (PrimMonad,PrimState)

import Data.Complex
import Data.Complex.Utils

import Foreign.Storable
import Foreign.C.Types (CDouble(..))

import qualified Data.Traversable as T
import qualified Data.Vector.Generic as GV
import qualified Data.Vector.Generic.Mutable as MV
import qualified Data.Vector.Unboxed as UV
import qualified Data.Vector.Storable as SV

-- | Freeze and convert a storable, mutable vector.
freezeConvert :: (PrimMonad m, Functor m, Storable a, GV.Vector v b, GV.Vector v a)
              => (a -> b) 
              -> GV.Mutable SV.Vector (PrimState m) a 
              -> m (v b)
freezeConvert f mvals = (GV.map f . GV.convert) <$> SV.unsafeFreeze mvals

-- | Freeze a Hermitian eigensystem as returned by typical LAPACK calls.
freezeHermitianES :: (Real a, Storable a, RealFloat a', Storable a'
                     ,Fractional b, Fractional b'
                     ,PrimMonad m, Functor m
                     ,GV.Vector v a
                     ,GV.Vector v b
                     ,GV.Vector v' (FComplex a')
                     ,GV.Vector v' (Complex b'))
                  => (GV.Mutable SV.Vector (PrimState m) a
                     ,Maybe (GV.Mutable SV.Vector (PrimState m) (FComplex a')))
                  -> m (v b, Maybe (v' (Complex b')))
freezeHermitianES (mvals, maybeMVecs) = do
  vals <- freezeConvert realToFrac mvals
  maybeVecs <- T.mapM (freezeConvert (complexToFrac . fromFComplex)) maybeMVecs
  return (vals, maybeVecs)

newtype instance UV.MVector s CDouble = MV_CDouble (UV.MVector s Double)
newtype instance UV.Vector CDouble = V_CDouble (UV.Vector Double)

instance MV.MVector UV.MVector CDouble where
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicOverlaps #-}
  {-# INLINE basicUnsafeNew #-}
  {-# INLINE basicUnsafeReplicate #-}
  {-# INLINE basicUnsafeRead #-}
  {-# INLINE basicUnsafeWrite #-}
  {-# INLINE basicClear #-}
  {-# INLINE basicSet #-}
  {-# INLINE basicUnsafeCopy #-}
  {-# INLINE basicUnsafeGrow #-}
  basicLength (MV_CDouble v) = MV.basicLength v
  basicUnsafeSlice i n (MV_CDouble v) = MV_CDouble $ MV.basicUnsafeSlice i n v
  basicOverlaps (MV_CDouble v1) (MV_CDouble v2) = MV.basicOverlaps v1 v2
  basicUnsafeNew n = MV_CDouble `liftM` MV.basicUnsafeNew n
  basicUnsafeReplicate n (CDouble z) = MV_CDouble `liftM` MV.basicUnsafeReplicate n z
  basicUnsafeRead (MV_CDouble v) i = CDouble `liftM` MV.basicUnsafeRead v i
  basicUnsafeWrite (MV_CDouble v) i (CDouble z) = MV.basicUnsafeWrite v i z
  basicClear (MV_CDouble v) = MV.basicClear v
  basicSet (MV_CDouble v) (CDouble z) = MV.basicSet v z
  basicUnsafeCopy (MV_CDouble v1) (MV_CDouble v2) = MV.basicUnsafeCopy v1 v2
  basicUnsafeMove (MV_CDouble v1) (MV_CDouble v2) = MV.basicUnsafeMove v1 v2
  basicUnsafeGrow (MV_CDouble v) n = MV_CDouble `liftM` MV.basicUnsafeGrow v n

instance GV.Vector UV.Vector CDouble where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_CDouble v) = V_CDouble `liftM` GV.basicUnsafeFreeze v
  basicUnsafeThaw (V_CDouble v) = MV_CDouble `liftM` GV.basicUnsafeThaw v
  basicLength (V_CDouble v) = GV.basicLength v
  basicUnsafeSlice i n (V_CDouble v) = V_CDouble $ GV.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_CDouble v) i
                = CDouble `liftM` GV.basicUnsafeIndexM v i
  basicUnsafeCopy (MV_CDouble mv) (V_CDouble v)
                = GV.basicUnsafeCopy mv v
  elemseq (V_CDouble v) (CDouble x) = GV.elemseq v x

instance UV.Unbox CDouble
