{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE OverlappingInstances #-}
{-# LANGUAGE DeriveDataTypeable #-}

module ComplexUtils (FComplex(..), FDComplex, FComplexable(..), Conjugable(..)) where

import Foreign.Ptr (castPtr)
import Foreign.C.Types (CDouble)
import Foreign.Storable

import Data.Complex

import Data.Typeable (Typeable)
import Data.Data (Data)

-- | Wrapping of 'Data.Complex' 'Complex' type for explicit memory structure
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

-- | Default double-precision complex type used in Fortran (complex*16 or double complex).
type FDComplex = FComplex CDouble

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
