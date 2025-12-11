!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: AddSourceTerms.F90                             !
!    CONTAINS: subroutine AddSourceTerms                  !
!                                                         !
!    PURPOSE: Adds source terms to advection-diffusion    !
!             equations.                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine AddSourceTerms

    implicit none

    ! Dummy routine - no need to add zeros to explicit terms
    ! Compiler will most probably get rid of this routine during optimization
    ! Useful to offer a "hook" for custom source terms (ex. channel flow forcing)

end subroutine AddSourceTerms  