c     Calculate the QH and QC in the thermoelectric heat pump (TEHP)
      subroutine test_TEHP
        use CommonDef
        real*8 :: Q1, Q2, T1, T2
        type(TE) :: TEC2
        type(StreamProp), dimension(2) :: S1, S2
        call get_TEHP(TEC2)
        call get_Stream(2, S1)
        call get_Stream(2, S2)
        T1 = S1(1)%T
        T2 = S2(1)%T
        call TEHP(TEC2, 1.0d0, T1, T2, Q1, Q2)
        call HeatBalance(Q1, S1(1), S1(2))
        call HeatBalance(Q2, S2(1), S2(2))
        write(*, "(3(A, F8.2))") " Q1 = ", Q1, " S1TIN = ", S1(1)%T, 
     &  " S1TOUT = ", S1(2)%T
        write(*, "(3(A, F8.2))") " Q2 = ", Q2, " S2TIN = ", S2(1)%T, 
     &  " S2TOUT = ", S2(2)%T
        pause
      end subroutine test_TEHP
      
      subroutine HeatBalance(Q, SIN, SOUT)
c     calculate the effluent
        use CommonDef
        type(StreamProp), intent(inout) :: SIN, SOUT
        real*8, intent(in) :: Q
        real*8 :: cp, WIN, WOUT, TIN, TOUT
        WIN = SIN%W
        TIN = SIN%T
        cp = SIN%PhysProp%cp
        TOUT = TIN + Q/(WIN*cp)
        SOUT%T = TOUT
      end subroutine HeatBalance
      
      subroutine get_Stream(numstm, streams)
        use CommonDef
        type(StreamProp), dimension(numstm), intent(inout) :: streams
        integer :: i
        do i = 1, numstm
          call initStream(streams(i))
          streams(i)%W = 1.d0 ! mass flowrate [kg/s]
          streams(i)%T = 3.2315d2 ! temperature [K]
          streams(i)%PhysProp%cp = 4.18d3 ! specific heat [J/kg-K]
        end do
      end subroutine get_Stream
      
      subroutine get_TEHP(TEC)
        use CommonDef
        type(TE), intent(inout) :: TEC
        TEC%ID = 'xxxx'
        TEC%Resistance = 1.2d0
        TEC%alpha = 1.0d0
        TEC%kappa = 1.0d0
      end subroutine get_TEHP
      
      subroutine TEHP(TEC, current, TH, TC, QH, QC)
        use CommonDef
        type(TE), intent(in) :: TEC ! thermoelectric cooler 
        real*8, intent(in) :: current,    ! current (A)
     &                          TH,   ! hot-side temperature (K)
     &                          TC    ! cold-side temperature (K)
        real*8, intent(out) :: QH,  ! heat release (W)
     &                           QC   ! absorbed heat (W)
        real*8 :: JouleHeat, PeltierEffect, ConductiveHeat
        real*8 :: R, a, k
        
        R = TEC%Resistance
        a = TEC%alpha
        k = TEC%kappa
        
        JouleHeat = current*current*R
        ConductiveHeat = k*(TH-TC)
        
        QH =  half*JouleHeat + a*current*TH - ConductiveHeat
        QC = -half*JouleHeat + a*current*TC - ConductiveHeat
        QC = -QC
        
      end subroutine