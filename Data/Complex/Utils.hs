{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE OverlappingInstances #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE TypeFamilies #-}

module Data.Complex.Utils 
  (FComplex(..)
  ,FDComplex
  ,FComplexable(..)
  ,Conjugable(..)
  ,fromFComplex
  ,complexToFrac)
where

import Foreign.Ptr (castPtr)
import Foreign.C.Types (CDouble)
import Foreign.Storable

import Data.Complex

import Data.Typeable (Typeable)
import Data.Data (Data)

import Control.Monad (liftM)

import qualified Data.Vector.Generic as GV
import qualified Data.Vector.Generic.Mutable as MV
import qualified Data.Vector.Unboxed as UV

-- | Wrapping of 'Data.Complex.Complex' type for explicit memory structure
-- for Fortran FFI compatibility.
--
-- The explicit structure is the simplest; two memory-adjacent storable
-- numbers. Note that for most uses of the FFI the value will have to be
-- converted to doubles for writing to memory.
newtype FComplex a = FComplex (Complex a) deriving (Show,Eq,Num,Fractional,Floating,Typeable,Data)

instance (RealFloat a, Storable a) => Storable (FComplex a) where
    sizeOf (FComplex z) = 2 * sizeOf (realPart z)
    alignment (FComplex z) = alignment (realPart z)
    peek p = do
      let q = castPtr p
      r <- peek q
      c <- peekElemOff q 1
      return $ FComplex (r :+ c)
    poke p (FComplex (r :+ c)) = do
      let q = castPtr p
      poke q r
      pokeElemOff q 1 c

newtype instance UV.MVector s (FComplex a) = MV_FComplex (UV.MVector s (Complex a))
newtype instance UV.Vector (FComplex a) = V_FComplex (UV.Vector (Complex a))

instance (UV.Unbox (Complex a)) => MV.MVector UV.MVector (FComplex a) where
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
  basicLength (MV_FComplex v) = MV.basicLength v
  basicUnsafeSlice i n (MV_FComplex v) = MV_FComplex $ MV.basicUnsafeSlice i n v
  basicOverlaps (MV_FComplex v1) (MV_FComplex v2) = MV.basicOverlaps v1 v2
  basicUnsafeNew n = MV_FComplex `liftM` MV.basicUnsafeNew n
  basicUnsafeReplicate n (FComplex z) = MV_FComplex `liftM` MV.basicUnsafeReplicate n z
  basicUnsafeRead (MV_FComplex v) i = FComplex `liftM` MV.basicUnsafeRead v i
  basicUnsafeWrite (MV_FComplex v) i (FComplex z) = MV.basicUnsafeWrite v i z
  basicClear (MV_FComplex v) = MV.basicClear v
  basicSet (MV_FComplex v) (FComplex z) = MV.basicSet v z
  basicUnsafeCopy (MV_FComplex v1) (MV_FComplex v2) = MV.basicUnsafeCopy v1 v2
  basicUnsafeMove (MV_FComplex v1) (MV_FComplex v2) = MV.basicUnsafeMove v1 v2
  basicUnsafeGrow (MV_FComplex v) n = MV_FComplex `liftM` MV.basicUnsafeGrow v n
  

instance (RealFloat a, UV.Unbox a) => GV.Vector UV.Vector (FComplex a) where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_FComplex v) = V_FComplex `liftM` GV.basicUnsafeFreeze v
  basicUnsafeThaw (V_FComplex v) = MV_FComplex `liftM` GV.basicUnsafeThaw v
  basicLength (V_FComplex v) = GV.basicLength v
  basicUnsafeSlice i n (V_FComplex v) = V_FComplex $ GV.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_FComplex v) i
                = FComplex `liftM` GV.basicUnsafeIndexM v i
  basicUnsafeCopy (MV_FComplex mv) (V_FComplex v)
                = GV.basicUnsafeCopy mv v
  elemseq (V_FComplex v) (FComplex x) = GV.elemseq v x

instance (RealFloat a, UV.Unbox a) => UV.Unbox (FComplex a)

-- | Default double-precision complex type used in Fortran (complex*16 or double complex).
type FDComplex = FComplex CDouble

-- | Convert from a 'FComplex' value to a Haskell complex number.
fromFComplex :: FComplex b -> Complex b
fromFComplex (FComplex c) = c

-- | Utility class for conversion to storable complex values.
class FComplexable a b where
    toFComplex :: (RealFloat b, Storable b) => a -> FComplex b

instance (Storable a) => FComplexable (FComplex a) a where
    toFComplex = id

instance (Real a) => FComplexable (Complex a) b where
    toFComplex (r :+ c) = FComplex (realToFrac r :+ realToFrac c)

instance (Real a) => FComplexable a b where
    toFComplex r = FComplex (realToFrac r :+ 0)

-- | Utility class for generic conjugation of number types.
class (Num a) => Conjugable a where
  cconj :: a -> a

instance (Real a) => Conjugable a where
  cconj = id

instance (RealFloat a) => Conjugable (Complex a) where
  cconj = conjugate

instance (RealFloat a) => Conjugable (FComplex a) where
  cconj (FComplex c) = FComplex $ conjugate c

-- | Utility function to take complex numbers to arbitrary numeric types in
-- the components.
complexToFrac :: (Real a, Fractional b) => Complex a -> Complex b
complexToFrac (r :+ c) = realToFrac r :+ realToFrac c
