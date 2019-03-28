module excitedsubs
   implicit none

contains
!
!
#include "ExcProp.f90"
#include "fcaApp.f90"
#include "GenerateDensities.f90"

#include "addInt.f90"
#include "Alpha_calc.f90"
#include "Ap_calculate.f90"
#include "basis_deinit.f90"
#include "basis_init.f90"
#include "calc2eFITT.f90"

#include "ChangeBasisF.f90"
#include "diagonH.f90"
#include "error.f90"
#include "initvec.f90"
#include "matMOtomatAO.f90" 

#include "MtoIANV.f90"
#include "multlr.f90"
#include "new_vectors.f90"
#include "norma.f90"
#include "ObtainOsc.f90"
#include "OscStr.f90"
#include "Pbeta_calc.f90"
#include "PCG_solve.f90"

#include "Prec_calculate.f90"
#include "PrintResults.f90"
#include "QRfactorization.f90"
#include "RCalculate.f90"
#include "RelaxedDensity.f90"
#include "residual.f90"
#include "RitzObtain.f90"

#include "subspaceMat.f90"
#include "tda.f90"
#include "total_fock.f90"
#include "TransDipole.f90"
#include "UnDiffDens.f90"
#include "vecMOtomatMO.f90"
#include "VecToMat.f90"
#include "XmatForm.f90"
#include "Zvector.f90"

#include "forcesexc.f90"
#include "Wcalculate.f90"
#include "WSgradcalc.f90"
#include "HVgradcalc.f90"
#include "CoulombForce.f90"
#include "varcoef_calc.f90"
#include "int2GExc.f90"
#include "int3GExc.f90"


! open shell linear response
#include "open_ExcProp.f90"
#include "open_tda.f90"
#include "open_basis_init.f90"
#include "open_basis_deinit.f90"
#include "open_matMOtomatAO.f90"


end module excitedsubs