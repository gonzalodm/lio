included :=
included += ExcProp.f90
included += fcaApp.f90
included += GenerateDensities.f90

included += addInt.f90
included += Alpha_calc.f90
included += Ap_calculate.f90
included += basis_deinit.f90
included += basis_init.f90
included += calc2eFITT.f90
included += ChangeBasisF.f90
included += diagonH.f90
included += error.f90
included += initvec.f90
included += matMOtomatAO.f90
included += MtoIANV.f90
included += multlr.f90
included += new_vectors.f90
included += norma.f90
included += ObtainOsc.f90
included += OscStr.f90
included += Pbeta_calc.f90
included += PCG_solve.f90
included += Prec_calculate.f90
included += PrintResults.f90
included += QRfactorization.f90
included += RCalculate.f90
included += RelaxedDensity.f90
included += residual.f90
included += RitzObtain.f90
included += subspaceMat.f90
included += tda.f90
included += total_fock.f90
included += TransDipole.f90
included += UnDiffDens.f90
included += vecMOtomatMO.f90
included += VecToMat.f90
included += XmatForm.f90
included += Zvector.f90
included += forcesexc.f90
included += Wcalculate.f90
included += WSgradcalc.f90
included += HVgradcalc.f90
included += CoulombForce.f90
included += varcoef_calc.f90
included += int2GExc.f90
included += int3GExc.f90

# open shell Linear response
included += open_ExcProp.f90
included += open_tda.f90
included += open_basis_init.f90
included += open_basis_deinit.f90
included += open_matMOtomatAO.f90
included += open_initvec.f90
included += open_MtoIANV.f90
included += open_subspaceMat.f90
included += open_RitzObtain.f90
included += open_residual.f90
included += open_new_vectors.f90
included += open_PrintResults.f90
included += open_OscStr.f90
included += open_calc2eFITT.f90



$(OBJPATH)/excitedsubs.o : $(included) excitedsubs.mk
