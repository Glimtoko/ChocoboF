
!     Last change:  JMM  2nd July 2004 
! 13th July saleadvv, saleaadvm added 
! 7th July sale added 
! 2nd July equipotential solver added
! 2nd July timestep changed to do perpendicular distances
! 26th June anti hourglass subroutine added
! 23rd June artificial viscosity edge christensen
! 23rd june volume intdiv put into subroutine logical zintdivvol
! 22nd more calls to femel and intdiv 
! 21st june axisymmetric div changed for artificial viscosity
!     2D measure required only no extra radial
! 7th June debug files removed
! 3rd june total energy test added
! 17th May axisymmetric - artif visc and CFL
!                       - work on sqrt(area) not vol now  
! 12th May axisymmetric added- logical zaxis
!                            - changed volume
!                            - changed div for a.v.
!                            - changed energy 
! 7th May disjoint nodes- momentum and accn changed
! 6th May disjoint nodes
!        - altered mergenodes to deal with 
!           different meshx ratio 2,3.
! 5th May regions finally working for nreg=6, changed artificial
!          visc back to sqrt(area)divu
!  
!21st April added region interface connections
!          - included altering nodelist, node-node array, momentum
!          - mergenodes subroutine added to deal with listing the
!          - merged interface nodes and adding them to node-node
!         - and element-element array
!16th march added namelist for geom variables,
!                              time parameters and eos 
!15th March realsied CFL not working played with tstep
!15th march changed volume calc to jac
!2nd march growth 1.02 dtinit=1e-4,dtmax=1e-3, cl=0.1 cq=1.0
! 1st march div working change cl to 0.06 dsc 0.1
! 1.00pm 26th feb viscosity-   added sqrt(area)*divv    
! 10.00am 26th feb viscosity-  added sqrt(area)deltav/deltax 
! 25th feb viscosity- vel diff
! 24th feb changed artificial viscosity to correct div error 
! 23rd feb corrected sound speed      
! 16th Feb added artificial viscosity bulk     
! 13th feb initial conditions sods shock tube                             
! 12th Feb energy subroutine     
! 11th,10th Feb femel, momentum, eos, accn, density subroutines, subroutines put
! in general dummy variables, predictor corrector drive written.   
! 7th Feb allocation of physical variables
!         initialise physical variables
!         calc volume
!         calc mass
!         calc element sound speed
!         calc stable timestep
! 22nd Jan, node-node array
!          changed variables in line with Andys 
!          element-element: removed one if loop that checked face back from node 
!                          also removed do loop through general cells 
!                          as cartesian ordered.     
! 20th Jan module connectivity added calculates element-element, node-node
! 19th jan geom finished  without bc arrays, with x, y, element-node
! ===================================================================
            
MODULE GLOBALCONSTANTS
! defines double precision
IMPLICIT NONE

INTEGER, PARAMETER        :: DP=selected_real_kind(15,300)

REAL (KIND=DP), PARAMETER :: zero = 0.00000000000000000_DP,             &
                              one = 1.00000000000000000_DP,             &
                              two = 2.00000000000000000_DP,             &
                            three = 3.00000000000000000_DP,             &
                             four = 4.00000000000000000_DP,             &
                             five = 5.00000000000000000_DP,             &
                              six = 6.00000000000000000_DP,             &
                            seven = 7.00000000000000000_DP,             &
                            eight = 8.00000000000000000_DP,             &
                             nine = 9.00000000000000000_DP,             &
                              ten = 10.0000000000000000_DP,             &
                           twelve = 12.0000000000000000_DP,             &  
                           thirty = 30.0000000000000000_DP,             &
                       thirtyfive = 35.0000000000000000_DP,             &
                             half = 5.00000000000000000e-1_DP,          &                       
                          quarter = 2.50000000000000000e-1_DP,          &
                               pi = 3.14159265358979324_DP,             &
                            twopi = two*pi                           
COMPLEX (KIND=DP), PARAMETER :: im =(0.00000000000000000_DP,             &
                                    1.00000000000000000_DP)                               
END MODULE GLOBALCONSTANTS
