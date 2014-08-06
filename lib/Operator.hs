{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}

module Operator (EndoOperator(..)) where

-- | Class for an operator from a vector space into itself.
class EndoOperator o where
  -- | Type of the eigenvalue. In general (Complex Double), but some
  -- operators may have real, positive, etc. eigenvectors.
  type Eigenvalue o :: *

  -- | Array type for eigenvalue storage. Usually 'Data.Vector.Unboxed' is useful for numeric types.
  type EigenvalueStorage o :: * -> *
  -- | Storage for eigenvectors.
  type EigenvectorStorage o :: *

  -- | Find the eigenvalues, optionally passing eigenvalue indices of interest (both inclusive).
  eigvals :: o -> Maybe (Int, Int) -> (EigenvalueStorage o) (Eigenvalue o)
  -- | Find the eigenvectors, optionally passing eigenvalue indices of interest (both inclusive).
  eigvecs :: o -> Maybe (Int, Int) -> EigenvectorStorage o

  -- | Equivalent to @('eigvals' o n, 'eigvecs' o n)@, but possibly more performant.
  eigsys :: o -> Maybe (Int, Int) -> ((EigenvalueStorage o) (Eigenvalue o), EigenvectorStorage o)
  eigsys o n = (eigvals o n, eigvecs o n)

  -- | In terms of the inner product space.
  adjoint :: o -> o
