      program main
      
!     test the integral solver
!      call test_dqk15      
!     test the molar concentration gradient     
!      call test_moduleVaporConc
!     test the conductivity
!      call test_conductivity
!     test the function to evaluate the molar latent heat
!      call test_MolarLatentHeat
!     test the function to evaluate derivation of pressure
!     call test_PressureDerivation_dP(T)
!     test the power function in dP(T) calculation
!      call test_powerfunction_L(T)
!     test the vapor pressure of saturated water (something wrong with the results)
!      call test_VaporPressure
!      call test_DiffVaporPressure
!     test the tortuosity(porosity)
!      call test_tortuosity(porosity)
!     test the effective thermal conductivity(error)
!      call test_EffThermalConductivity
!     test the Effective Molar Enthalpy
!      call test_EffectiveMolarEnthalpy
!     test the effective diffusion coefficient
!      call test_ EffectiveDiffusivity
!      call test_DAVG
!      call test_SetStreamIndex
!      call test_AvgMolWeight
!      call test_SpecVolV
!      call test_CalcProfile
      call test_TEHP

      end program main