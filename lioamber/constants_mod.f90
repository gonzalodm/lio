#include "datatypes/datatypes.fh"
module constants_mod
  implicit NONE
    LIODBLE :: pi32,pi,rpi,pi5,pi52,pis32,piss,rpis,pis5,pis52
    parameter( pi32 = 5.56832799683170698D0, &
            &  pi   = 3.14159265358979312D0, &
            &  rpi  = 1.77245385090551588D0, &
            &  pi5  = 34.9868366552497108D0, &
            &  pi52 = 34.9868366552497108D0 )

    parameter( pis32= 5.56832799683170698E0, &
            &  piss = 3.14159265358979312E0, &
            &  rpis = 1.77245385090551588E0, &
            &  pis5 = 34.9868366552497108E0, &
            &  pis52= 34.9868366552497108E0 )
    LIODBLE, parameter :: bohr = 0.529177d0
    LIODBLE, parameter :: massprot_elec = 1822.88857149d0
    LIODBLE, parameter :: H_to_eV = 27.211396132d0
end module constants_mod
