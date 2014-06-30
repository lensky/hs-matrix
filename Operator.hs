{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}

module Operator (EndoOperator(..)) where

-- | Class for an operator from a vector space into itself.
class EndoOperator o where
  type Eigenvalue o :: *      

  type EigenvalueStorage o :: * -> *
  type EigenvectorStorage o :: *

  eigvals :: o -> Maybe (Int, Int) -> (EigenvalueStorage o) (Eigenvalue o)
  eigvecs :: o -> Maybe (Int, Int) -> EigenvectorStorage o

  eigsys :: o -> Maybe (Int, Int) -> ((EigenvalueStorage o) (Eigenvalue o), EigenvectorStorage o)
  eigsys o n = (eigvals o n, eigvecs o n)

  adjoint :: o -> o
