#include "preproc.h"

module types
    implicit none

    type tInfo
        integer :: dimen
        real(typeReal), allocatable :: cost (:,:)
        
        integer :: REINSERTION  = 1
        integer :: OR_OPT_2     = 2
        integer :: OR_OPT_3     = 3 
        integer :: SWAP         = 4
        integer :: TWO_OPT      = 5

        real(typeReal) :: fmax = 3.4028235E+38
        
        integer :: reinsert_call = 0

        integer, allocatable :: rnd(:)
        integer :: rnd_index = 1
    end type

    type tSeqInfo
        real(typeReal) :: T, C, W
    end type

    type tSolution
        real(typeReal) :: cost
        integer :: s_size
        integer, allocatable :: s (:)
        type(tSeqInfo), allocatable :: seq (:,:)
    end type

end module

