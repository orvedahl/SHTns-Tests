f90sources += cheb_grids.f90
f90sources += Legendre_Polynomials.F90
ifdef shtns
  f90sources += Legendre_Transforms_SHTns.F90
else
  f90sources += Legendre_Transforms.F90
endif
f90sources += main.f90
f90sources += probin.f90
f90sources += Structures.F90
f90sources += test_suite.f90
