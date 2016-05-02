  module mSolverPardiso

        integer, allocatable :: listaDosElemsPorNoCRS(:,:)

        integer :: posPonteiro, contPonteiro, posColunas, posCoef
        integer :: lda, nonzeros, nonzerosEst

      contains
!
!=======================================================================
!    
      subroutine solverPardisoEsparso(alhs, b, Ap, Ai, pt_, iparm_, dparm_, neq, nonzeros, &
     &                                            simetria, label,parte)

      
      implicit none 
!
      real*8  :: alhs(*), b(neq)
      integer :: Ap(*), Ai(*)
      integer, intent(in) :: neq, nonzeros
      logical, intent(in) :: simetria
      character(LEN=3) :: label
      character(LEN=4) :: parte
      INTEGER pt_(64), iparm_(64)
      REAL*8  dparm_(64)
! 
      real*8, allocatable  :: xGeo(:), xVel(:)

      if(label=='vel') then
         if(.not.allocated(xVel)) allocate(xVel(neq))
         xVel=0.d0
         call solverPardisoPPD(Ap, Ai, pt_, iparm_, dparm_, alhs, xVel, b, neq, nonzeros, simetria, label,parte)
      endif

      
      if(label=='geo') then
         if(.not.allocated(xGeo)) allocate(xGeo(neq))
         xGeo=0.d0
         call solverPardisoPPD(Ap, Ai, pt_, iparm_, dparm_, alhs, xGeo, b, neq, nonzeros, simetria, label,parte)
      endif   
           

      end subroutine solverPardisoEsparso
!
!=======================================================================
!  
      subroutine solverPardisoPPD(ia, ja, pt_, iparm_, dparm_, a, x, b, neq, nonzeros, simetria, label, parte)

      implicit none
!
        INTEGER, intent(in)  ::  ia(neq+1), ja(nonzeros)
        INTEGER pt_(64), iparm_(64)
        REAL*8  dparm_(64)
        
        REAL*8, INTENT(IN)   :: a(nonzeros)
        REAL*8, INTENT(INOUT):: b(neq)
        REAL*8, INTENT(INOUT)  :: x(neq)
        INTEGER, INTENT(IN)  :: NEQ, NONZEROS
        LOGICAL, INTENT(IN)  :: simetria
        character(LEN=3), INTENT(IN) :: label
        character(LEN=4), INTENT(IN) :: parte
!
        REAL*8,  save :: ddum
        INTEGER, save :: maxfct, mnum, mtype, phase, nrhs, error,  msglvl
        INTEGER, save :: idum, solver
        INTEGER :: i
        REAL*8  :: t1, t2 , tt1, tt2
        integer :: omp_get_num_threads
!
#ifdef withPardiso
!
!  .. Setup Pardiso control parameters und initialize the solvers     
!     internal adress pointers. This is only necessary for the FIRST   
!     call of the PARDISO solver.                                     
!     
! !$OMP PARALLEL

      if(parte=='reor'.or.parte=='full') then

         iparm_=0
         mtype     = -2   ! real and symmetric matrix, inefinite
         iparm_(11) = 0    !Do not use (symmetric matrices).
         if(simetria.eqv..false.)  mtype     = 11   ! unsymmetric matrix, indefinite

         solver     = 0    ! use 0 for sparse direct method or 1 multi-recursive iterative solver
         msglvl     = 0    ! with statistical information
         iparm_(33) = 0    ! compute determinant 
         iparm_(52) = 1    !For OpenMP-threaded solver
         iparm_(2)  = 0 !ou 0?      !Fill-In reduction reordering.
         iparm_(27) = 1
      
         nrhs      = 1
         mnum      = 1
         pt_       = 0
         idum      = 0
         ddum      = 0.0
         dparm_    = 0.0
         maxfct    = 1
         error     = 0
!
!  .. Numbers of Processors ( value of OMP_NUM_THREADS )
!
         iparm_(3) = 1

#ifdef withOMP
!$OMP PARALLEL
         iparm_(3) =   omp_get_num_threads()
!$OMP END PARALLEL
#endif
         print*, "em pardiso com", iparm_(3), "numThreads"
!        
         call timing(tt1)
!
!  .. PARDISO license check and initialize solver
!       call pardisoinit(pt, mtype, solver, iparm, dparm, error)
!
         IF (error .NE. 0) THEN
            IF (error.EQ.-10 ) WRITE(*,*) 'No license file found'
            IF (error.EQ.-11 ) WRITE(*,*) 'License is expired'
            IF (error.EQ.-12 ) WRITE(*,*) 'Wrong username or hostname'
            STOP ' (error .NE. 0) '
         ELSE
            WRITE(*,*) '[PARDISO]: License check was successful ... '
         END IF
! 
!..   Reordering and Symbolic Factorization, This step also allocates
!     all memory that is necessary for the factorization
!
         WRITE(*,*) 'Begining reordering ... ', label
         call timing(t1)
         phase     = 11     ! only reordering and symbolic factorization

         CALL pardiso (pt_, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                    idum, nrhs, iparm_, msglvl, ddum, ddum, error, dparm_)
         call timing(t2)
! #ifdef mostrarTempos
!        write(*,*) "reordering: ", label, ", tempo de parede = ", t2 - t1
! #endif

         IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP
         END IF

!        WRITE(*,*) 'Number of nonzeros in factors   = ',iparm(18)
!        WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

      end if  !if(parte=='reor'.or.parte=='full') then


      if(parte=='fact'.or. parte=='full') then
!
!.. Factorization.
!
         WRITE(*,*) 'Begining factorization  ... '
         call timing(t1)
         phase     = 22  ! only factorization
         CALL pardiso (pt_, maxfct, mnum, mtype, phase, neq, a, ia, ja,  &
                    idum, nrhs, iparm_, msglvl, ddum, ddum, error, dparm_) 
         call timing(t2)
! #ifdef mostrarTempos
!       write(*,*) "factorization: ", label, ", tempo = ", t2 - t1
! #endif

         IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP
         ENDIF 

     endif  !     if(parte=='fact'.or. parte=='full') then


      if(parte=='back'.or. parte=='full') then

!.. Back substitution and iterative refinement
         WRITE(*,*) 'Begining backsubstitution  ... '
         call timing(t1)
         iparm_(8) = 1   ! max numbers of iterative refinement steps
         phase     = 33  ! only solve
         iparm_(6) = 1 
         CALL pardiso (pt_, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                   idum, nrhs, iparm_, msglvl, b, x, error, dparm_) 
         call timing(t2)
         
         
! #ifdef mostrarTempos
!       write(*,*) "backsubstitution: ",label, ", tempo = ", t2 - t1
! #endif

         call timing(tt2)
! #ifdef mostrarTempos
!       WRITE(*,*) 'Solve completed ...  ',label, ", tempo =", tt2-tt1
! #endif


!!!!!!!!!!  ITERATIVO  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! 
!        WRITE(*,*) 'Precondicionador  ... '

!        msglvl = 0 !0 para nao imprimir as estatisticas
!        solver    = 1
!        iparm(4)=0
!        phase = 12 ! set up of the preconditioner
!        CALL pardiso (pt , maxfct , mnum , mtype , phase , neq, a, ia , ja ,& 
!         idum , nrhs , iparm , msglvl , ddum , ddum , error , dparm )
! 
! 
! ! .. Iterative Solve
!        WRITE(*,*) 'solver ... '
!        iparm(4) = 61
!        phase = 23 ! only solve
!        CALL pardiso (pt , maxfct , mnum , mtype , phase , neq, a, ia , ja , &
!           idum , nrhs , iparm , msglvl , b, x, error , dparm )
! 
!          tt2 = omp_get_wtime( );
!       WRITE(*,*) 'Solve completed ...  ', tt2-tt1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       desalocando memoria
!       WRITE(*,*) 'Begining memory desalocation  ... '
!       call timing(t1)
!       phase = -1
!       CALL pardiso (pt_, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
!                    idum, nrhs, iparm_, msglvl, b, x, error, dparm_)
!       call timing(t2)
! #ifdef mostrarTempos
!       write(*,*) "memory desalocation: ", label, ", tempo = ", t2 - t1
! #endif

!       b(1:neq) = x(1:neq)

      endif !if(parte=='back'.or. parte=='full') then

#endif

      end subroutine solverPardisoPPD
!
!=======================================================================
!    
      subroutine criarPonteirosMatEsparsa_CRS(nsd, ndof, neq, numCoefPorLinha,  &
              conectsElem, listaDosElems, id, LMstencilEq_, Ap_, Ai_, numConexoes,  &
              nen, numConexoesPorElem, nonzeros, simetria)

      implicit none 
!
      integer, intent(in) :: nsd, neq, numConexoes, nen, numConexoesPorElem,nonzeros,ndof
      logical, intent(in) :: simetria
      integer, intent(out) :: numCoefPorLinha
      integer :: conectsElem(nen,*), listaDosElems(numConexoesPorElem,*)
      integer :: id(ndof,*)
      integer, allocatable :: LMstencilEq_(:,:), Ap_(:), Ai_(:)
      
!
       if(numCoefPorLinha==26.or.numCoefPorLinha==18) then
!
         if(.not.allocated(LMstencilEq_))allocate(LMstencilEq_(neq,numCoefPorLinha)); 
         LMstencilEq_=0

         call montarLmStencilGeo_CRS    (LMstencilEq_,listaDosElems, id, &
                   conectsElem, numCoefPorLinha, ndof, numConexoes, nen, numConexoesPorElem, neq, simetria)

         allocate(Ap_(neq+1));    Ap_=0  
         call montarPonteiroAp_CRS(Ap_, LMstencilEq_, numCoefPorLinha, neq, nonzeros)
!

         allocate(Ai_(nonzeros)); Ai_=0
         call montarPonteiroAi_CRS(Ai_, LMstencilEq_, numCoefPorLinha, neq, nonzeros)
!
      else
!
         if(.not.allocated(LMstencilEq_))allocate(LMstencilEq_(neq,numCoefPorLinha))
         LMstencilEq_=0

         call montarLmStencilVel_CRS    (LMstencilEq_,listaDosElems, id, &
                   conectsElem, numCoefPorLinha, ndof, numConexoes, numConexoesPorElem, neq, simetria)
!
         allocate(Ap_(neq+1));    Ap_=0  
         call montarPonteiroAp_CRS(Ap_, LMstencilEq_, numCoefPorLinha, neq, nonzeros)
!
         allocate(Ai_(nonzeros)); Ai_=0
         call montarPonteiroAi_CRS(Ai_, LMstencilEq_, numCoefPorLinha, neq, nonzeros)
!
      endif
! 
      contains
      !
      !**** new *************************************************************
      !
      subroutine montarLmStencilGeo_CRS(LMstencilEq, listaDosElems, id, &
        conectsElem, numCoefPorLinha, ndof, numConexoes, nen, numConexoesPorElem, neq, simetria) 
!  
      use mGlobaisEscalares, only: novaMalha
!
      implicit none
!
      integer, intent(in)    :: numCoefPorLinha, numConexoes,numConexoesPorElem,ndof,neq,nen
      integer, intent(inout) :: LMstencilEq(neq,numCoefPorLinha)
      integer, intent(in)    :: listaDosElems(numConexoesPorElem,*), conectsElem(nen,*)
      integer, intent(in)    :: id(ndof,*)
      logical, intent(in)    :: simetria
!
      integer :: nc, cont, nel, i, k, l, m,  numEq, dir
      logical :: propriaEq
      integer :: numCoefAux, numViz
      integer, allocatable :: LMstencilEqAux(:)
!
      if(novaMalha.eqv..true.)  numCoefAux=numCoefPorLinha+14
      if(novaMalha.eqv..false.) numCoefAux=numCoefPorLinha+8
      allocate(LMstencilEqAux(numCoefAux))

      if(novaMalha.eqv..true. ) numViz=6
      if(novaMalha.eqv..false.) numViz=nen

      numEq=0
      LMstencilEq=0
      LMstencilEqAux=0

      do nc=1, numConexoes

         do dir=1, ndof 
   
         if(id(dir,nc).ne.0) then
            numEq=numEq+1
            LMstencilEqAux=0
            cont=1
            propriaEq=.false.

            do k=1, numViz
               if(listaDosElems(k,nc).ne.0) then
                  nel=listaDosElems(k,nc)

                  do i=1, nen!4!numConexoesPorElem

                     if (id(1,conectsElem(i,nel)).ne.0.or.id(2,conectsElem(i,nel)).ne.0) then

                        if(id(1,conectsElem(i,nel)).ne.numEq .and. id(2,conectsElem(i,nel)).ne.numEq) then
                           
                           do m=1, ndof

                              if(simetria.eqv..true.) then
                                 if(id(m,conectsElem(i,nel))>=numEq) then
                                    LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                    cont=cont+1
                                 endif
                              else
                                 if(id(m,conectsElem(i,nel)).ne.0) then
                                    LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                    cont=cont+1
                                 endif
                              endif
                           enddo

                        else
                            if(propriaEq.eqv..false.) then
                               do m=1, ndof
                                  propriaEq=.true.
                                  if(simetria.eqv..true.) then
                                     if(id(m,conectsElem(i,nel))>=numEq) then
                                        LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                        cont=cont+1
                                     endif
                                  else
                                     if(id(m,conectsElem(i,nel)).ne.0) then
                                        LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                        cont=cont+1
                                     endif
                                  endif
                                  
                               enddo
                            endif

                        end if
                     end if
                  end do
!
               end if
            end do !k

         call ordenarLMstencil(LMstencilEqAux(:),numCoefAux)

         do i=1, numCoefAux-1
             if(LMstencilEqAux(i)==LMstencilEqAux(i+1)) LMstencilEqAux(i)=0
         end do

         call ordenarLMstencil(LMstencilEqAux(:),numCoefAux)

         if(novaMalha.eqv..true.)LMstencilEq(numEq,1:numCoefPorLinha)=LMstencilEqAux(15:numCoefAux)
         if(novaMalha.eqv..false.)LMstencilEq(numEq,1:numCoefPorLinha)=LMstencilEqAux(9:numCoefAux)

!            write(*,'(a,i8,a,26i8)'), "LmStencil", numEq, "->",LMstencilEq(numEq,1:numCoefPorLinha)
         end if

       end do !dir
      end do !gl

      deallocate(LMstencilEqAux)
      
      end subroutine montarLmStencilGeo_CRS

     !
      !**** new *************************************************************
      !
      subroutine montarLmStencilVel_CRS(LMstencilEq, listaDosElems, id, &
          conectsElem, numCoefPorLinha, ndof, numConexoes, numConexoesPorElem, neqV, simetria) 
!  
      implicit none
!
      integer, intent(in)    :: numCoefPorLinha, numConexoes,numConexoesPorElem,ndof,neqV
      integer, intent(inout) :: LMstencilEq(neqV,numCoefPorLinha)
      integer, intent(in)    :: listaDosElems(numConexoesPorElem,*), conectsElem(numConexoesPorElem,*)
      integer, intent(in)    :: id(ndof,*)
      logical, intent(in)    :: simetria
!
      integer :: nc, cont, nel, i, k, l, numEq
      integer :: LMstencilEqTemp(numCoefPorLinha)
      integer :: contLados, dir
      logical :: propriaEq
!
      numEq=0
      LMstencilEq=0

      do nc=1, numConexoes

         do dir=1, ndof 
   
         if(id(dir,nc).ne.0) then
            numEq=numEq+1
            cont=1
            contLados=0
            propriaEq=.false.


            do k=1, numConexoesPorElem
               if(listaDosElems(k,nc).ne.0) then
                  nel=listaDosElems(k,nc)

                  do i=1, numConexoesPorElem

                     if (id(dir,conectsElem(i,nel))==0) then
                        LMstencilEq(numEq,cont)=0
                     else
                        if(id(dir,conectsElem(i,nel)).ne.numEq) then
                           LMstencilEq(numEq,cont)=id(dir,conectsElem(i,nel))
                           cont=cont+1
                        end if
                     end if
                  end do
                  LMstencilEq(numEq,numCoefPorLinha)=numEq
               end if
            end do
         end if

       end do !dir
      end do


      do l=1, neqV

         LMstencilEqTemp=0
         cont=1
         LMstencilEqTemp(cont)=l
         do i = 1, numCoefPorLinha
            if(simetria.eqv..true.) then
               if( LMstencilEq(l,i)>=l )then !considera a simetria
                   LMstencilEqTemp(cont)=LMstencilEq(l,i)
                   cont=cont+1
               end if
            else        
               if( LMstencilEq(l,i).ne.0 )then !nao considera simetria!           
                  LMstencilEqTemp(cont)=LMstencilEq(l,i)
                  cont=cont+1
               end if
            endif
         end do
         LMstencilEq(l,:)=LMstencilEqTemp
         call ordenarLMstencil(LMstencilEq(l,:),numCoefPorLinha)

         do i=1, numCoefPorLinha-1
             if(LMstencilEq(l,i)==LMstencilEq(l,i+1)) LMstencilEq(l,i)=0
         end do

         call ordenarLMstencil(LMstencilEq(l,:),numCoefPorLinha)

!           write(*,'(a,i4,a,18i4)'), "LmStencil", l, "->",LMstencilEq(l,1:numCoefPorLinha)

      end do
!      
      end subroutine montarLmStencilVel_CRS
      !
      !----------------------------------------------------------------------
      !
      subroutine ordenarLMstencil(LMstencilEq,numCoefPorLinha)
      
      implicit none
!
      integer, intent(in)    :: numCoefPorLinha      
      integer, intent(inout) :: LMstencilEq(numCoefPorLinha)
!
      integer :: menorEq, n, nn , tmp

      do n = 1, numCoefPorLinha
        menorEq=n
        do nn = n+1, numCoefPorLinha
           if(LMstencilEq(nn)<LMstencilEq(menorEq)) menorEq = nn
        end do

        if(n == menorEq) cycle
        tmp                  = LMstencilEq(n)
        LMstencilEq(n)       = LMstencilEq(menorEq)
        LMstencilEq(menorEq) = tmp
                
      enddo

      end subroutine ordenarLMstencil
      !
      !**** new *************************************************************
      !
      subroutine montarPonteiroAp_CRS (Ap, LMstencilEq, numCoefPorLinha, neq, nonzeros)

      implicit none

      integer :: Ap(*)
      integer :: LMstencilEq(neq,numCoefPorLinha)
      integer :: numCoefPorLinha, nonzeros, neq
!
      integer :: l,j

      ! Montando Ap

      call montarListaPonteiros(Ap, LMstencilEq,neq,numCoefPorLinha)

      ! Contando os valores nao nulos
      nonzeros=0
      do l=2, neq+1
            nonzeros=nonzeros+(Ap(l)-Ap(l-1))
      end do

       end subroutine montarPonteiroAp_CRS
      !
      !**** new *************************************************************
      !
      subroutine montarPonteiroAi_CRS (Ai, LMstencilEq, numCoefPorLinha, neq, nonzeros)

      implicit none

      integer, intent(out) :: Ai(*)
      integer :: LMstencilEq(neq,numCoefPorLinha)
      integer :: numCoefPorLinha, nonzeros, neq
!
      integer :: p

      call montarListaIndices(Ai, LMstencilEq, neq, nonzeros, numCoefPorLinha)

      end subroutine montarPonteiroAi_CRS
      !
      !**** new *************************************************************
      !
      subroutine montarListaPonteiros(Ap, LMstencilEq, neq, numCoefPorLinha)

      implicit none
      
      integer, intent(out) :: Ap(*)
      integer, intent(in) :: neq, numCoefPorLinha,LMstencilEq(neq,numCoefPorLinha)
      integer :: LMstencilEqTemp(0:numCoefPorLinha)
!
      integer :: n, l

      posPonteiro=0
      contPonteiro=0
      do l=1, neq
         posPonteiro=posPonteiro+1
         Ap(posPonteiro)=contPonteiro+1
         LMstencilEqTemp=0
         LMstencilEqTemp(1:numCoefPorLinha)=LMstencilEq(l,:)

         do n = 1, numCoefPorLinha
            if(LMstencilEqTemp(n) == LMstencilEqTemp(n-1)) cycle
            contPonteiro= contPonteiro+1
         enddo    

         if(l==neq) then
            posPonteiro=posPonteiro+1
            Ap(posPonteiro)=contPonteiro+1
         end if
      end do

      end subroutine montarListaPonteiros
      !
      !**** new *************************************************************
      !
      subroutine montarListaIndices(Ai, LMstencilEq, neq, nonzeros, numCoefPorLinha)
!
      implicit none
      integer, intent(inout) :: Ai(nonzeros)
      integer, intent(in)    :: nonzeros, neq, numCoefPorLinha
      integer, intent(in)    :: LMstencilEq(neq,numCoefPorLinha)
!
      integer :: LMstencilEqTemp(0:numCoefPorLinha)
      integer :: i, n, posColunas
!
      posColunas=0
      do i=1, neq
         LMstencilEqTemp=0
         LMstencilEqTemp(1:numCoefPorLinha)=LMstencilEq(i,:)
         do n = 1, numCoefPorLinha
             if(LMstencilEqTemp(n).ne.LMstencilEqTemp(n-1).and.LMstencilEqTemp(n).ne.0 ) then
                 posColunas=posColunas+1
                 Ai(posColunas)= LMstencilEqTemp(n)
             end if
         enddo
      end do
!
      end subroutine montarListaIndices
!
      end subroutine criarPonteirosMatEsparsa_CRS


!
! **** new *********************************************************************
!
       subroutine addlhsCRS(alhs,eleffm,lm,Ap,Ai,nee)
! 
! .... program to add element left-hand-side matrix to
!         global left-hand-side matrix
! 
! 
       implicit none
! 
! .... deactivate above card(s) for single-precision operation
! 
       real*8  :: alhs(*),eleffm(nee,*)
       integer :: lm(*), Ap(*), Ai(*)
       integer :: nee
!
       integer :: i, j, linha, coluna
       integer :: inicio, fim, jj
! 
          do 400 j=1,nee
          coluna = lm(j)
          if (coluna.gt.0) then
! 
             do 200 i=1,nee
             linha = lm(i)

             if (linha.gt.0) then

                inicio=Ap(coluna)
                fim=Ap(coluna+1)-1
                 do jj=inicio, fim
                      if(Ai(jj)==linha) then
                          alhs(jj)= alhs(jj)+eleffm(i,j)
                      endif
                 enddo
             endif
  200       continue
! 
          endif
  400    continue
! 
       return
       end subroutine
!
!=======================================================================
!    
      subroutine criarPonteirosMatEsparsa_CRS_MR(nsd, ndof, neq, numCoefPorLinha,  &
              conectsElem, listaDosElems, id, LMstencilEq_, Ai_, Ap_, numConexoes, numConexoesPorElem,  &
              nonzeros, simetria)

      implicit none 
!
      integer, intent(in) :: nsd, neq, numConexoes, numConexoesPorElem,nonzeros,ndof
      logical, intent(in) :: simetria
      integer, intent(out) :: numCoefPorLinha
      integer :: conectsElem(numConexoesPorElem,*), listaDosElems(numConexoesPorElem,*)
      integer :: id(ndof,*)
      integer, allocatable :: LMstencilEq_(:,:), Ap_(:), Ai_(:)
!
       if(numCoefPorLinha==18) then
!
         if(.not.allocated(LMstencilEq_))allocate(LMstencilEq_(neq,numCoefPorLinha)); 
         LMstencilEq_=0

         call montarLmStencilGeo_CRS    (LMstencilEq_,listaDosElems, id, &
                   conectsElem, numCoefPorLinha, ndof, numConexoes, numConexoesPorElem, neq, simetria)

         allocate(Ap_(neq+1));    Ap_=0  
         call montarPonteiroAp_CRS(Ap_, LMstencilEq_, numCoefPorLinha, neq, nonzeros)
!
         allocate(Ai_(nonzeros)); Ai_=0
         call montarPonteiroAi_CRS(Ai_, LMstencilEq_, numCoefPorLinha, neq, nonzeros)
!
      else
!
         if(.not.allocated(LMstencilEq_))allocate(LMstencilEq_(neq,numCoefPorLinha))
         LMstencilEq_=0

         call montarLmStencilVel_CRS    (LMstencilEq_,listaDosElems, id, &
                   conectsElem, numCoefPorLinha, ndof, numConexoes, numConexoesPorElem, neq, simetria)

         allocate(Ap_(neq+1));    Ap_=0  
         call montarPonteiroAp_CRS(Ap_, LMstencilEq_, numCoefPorLinha, neq, nonzeros)
!
         allocate(Ai_(nonzeros)); Ai_=0
         call montarPonteiroAi_CRS(Ai_, LMstencilEq_, numCoefPorLinha, neq, nonzeros)
!
      endif
! 
      contains
      !
      !**** new *************************************************************
      !
      subroutine montarLmStencilGeo_CRS(LMstencilEq, listaDosElems, id, &
        conectsElem, numCoefPorLinha, ndof, numConexoes, numConexoesPorElem, neq, simetria) 
!  
      implicit none
!
      integer, intent(in)    :: numCoefPorLinha, numConexoes,numConexoesPorElem,ndof,neq
      integer, intent(inout) :: LMstencilEq(neq,numCoefPorLinha)
      integer, intent(in)    :: listaDosElems(numConexoesPorElem,*), conectsElem(numConexoesPorElem,*)
      integer, intent(in)    :: id(ndof,*)
      logical, intent(in)    :: simetria
!
      integer :: nc, cont, nel, i, k, l, m,  numEq, dir
      logical :: propriaEq
      integer :: numCoefAux
      integer, allocatable :: LMstencilEqAux(:)
!
      numCoefAux=numCoefPorLinha+8
      allocate(LMstencilEqAux(numCoefAux))

      numEq=0
      LMstencilEq=0
      LMstencilEqAux=0

      do nc=1, numConexoes

         do dir=1, ndof 
   
         if(id(dir,nc).ne.0) then
            numEq=numEq+1
            LMstencilEqAux=0
            cont=1
            propriaEq=.false.

            do k=1, numConexoesPorElem
               if(listaDosElems(k,nc).ne.0) then
                  nel=listaDosElems(k,nc)

                  do i=1, numConexoesPorElem

                     if (id(1,conectsElem(i,nel)).ne.0.or.id(2,conectsElem(i,nel)).ne.0) then

                        if(id(1,conectsElem(i,nel)).ne.numEq .and. id(2,conectsElem(i,nel)).ne.numEq) then
                           
                           do m=1, ndof

                              if(simetria.eqv..true.) then
                                 if(id(m,conectsElem(i,nel))>=numEq) then
                                    LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                    cont=cont+1
                                 endif
                              else
                                 if(id(m,conectsElem(i,nel)).ne.0) then
                                    LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                    cont=cont+1
                                 endif
                              endif
                           enddo

                        else
                            if(propriaEq.eqv..false.) then
                               do m=1, ndof
                                  propriaEq=.true.
                                  if(simetria.eqv..true.) then
                                     if(id(m,conectsElem(i,nel))>=numEq) then
                                        LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                        cont=cont+1
                                     endif
                                  else
                                     if(id(m,conectsElem(i,nel)).ne.0) then
                                        LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                        cont=cont+1
                                     endif
                                  endif
                                  
                               enddo
                            endif

                        end if
                     end if
                  end do
!
               end if
            end do !k

         call ordenarLMstencil(LMstencilEqAux(:),numCoefAux)

         do i=1, numCoefAux-1
             if(LMstencilEqAux(i)==LMstencilEqAux(i+1)) LMstencilEqAux(i)=0
         end do

         call ordenarLMstencil(LMstencilEqAux(:),numCoefAux)

         LMstencilEq(numEq,1:numCoefPorLinha)=LMstencilEqAux(9:numCoefAux)

         end if

       end do !dir
      end do !gl

      deallocate(LMstencilEqAux)
      
      end subroutine montarLmStencilGeo_CRS

     !
      !**** new *************************************************************
      !
      subroutine montarLmStencilVel_CRS(LMstencilEq, listaDosElems, id, &
          conectsElem, numCoefPorLinha, ndof, numConexoes, numConexoesPorElem, neqV, simetria) 
!  
      implicit none
!
      integer, intent(in)    :: numCoefPorLinha, numConexoes,numConexoesPorElem,ndof,neqV
      integer, intent(inout) :: LMstencilEq(neqV,numCoefPorLinha)
      integer, intent(in)    :: listaDosElems(numConexoesPorElem,*), conectsElem(numConexoesPorElem,*)
      integer, intent(in)    :: id(ndof,*)
      logical, intent(in)    :: simetria
!
      integer :: nc, cont, nel, i, k, l, numEq
      integer :: LMstencilEqTemp(numCoefPorLinha)
      integer :: contLados, dir
      logical :: propriaEq
!
      numEq=0
      LMstencilEq=0

      do nc=1, numConexoes

         do dir=1, ndof 
   
         if(id(dir,nc).ne.0) then
            numEq=numEq+1
            cont=1
            contLados=0
            propriaEq=.false.


            do k=1, numConexoesPorElem
               if(listaDosElems(k,nc).ne.0) then
                  nel=listaDosElems(k,nc)

                  do i=1, numConexoesPorElem

                     if (id(dir,conectsElem(i,nel))==0) then
                        LMstencilEq(numEq,cont)=0
                     else
                        if(id(dir,conectsElem(i,nel)).ne.numEq) then
                           LMstencilEq(numEq,cont)=id(dir,conectsElem(i,nel))
                           cont=cont+1
                        end if
                     end if
                  end do
                  LMstencilEq(numEq,numCoefPorLinha)=numEq
               end if
            end do
         end if

       end do !dir
      end do


      do l=1, neqV

         LMstencilEqTemp=0
         cont=1
         LMstencilEqTemp(cont)=l
         do i = 1, numCoefPorLinha
            if(simetria.eqv..true.) then
               if( LMstencilEq(l,i)>=l )then !considera a simetria
                   LMstencilEqTemp(cont)=LMstencilEq(l,i)
                   cont=cont+1
               end if
            else        
               if( LMstencilEq(l,i).ne.0 )then !nao considera simetria!           
                  LMstencilEqTemp(cont)=LMstencilEq(l,i)
                  cont=cont+1
               end if
            endif
         end do
         LMstencilEq(l,:)=LMstencilEqTemp
         call ordenarLMstencil(LMstencilEq(l,:),numCoefPorLinha)

         do i=1, numCoefPorLinha-1
             if(LMstencilEq(l,i)==LMstencilEq(l,i+1)) LMstencilEq(l,i)=0
         end do

         call ordenarLMstencil(LMstencilEq(l,:),numCoefPorLinha)

!           write(*,'(a,i4,a,18i4)'), "LmStencil", l, "->",LMstencilEq(l,1:numCoefPorLinha)

      end do
!      
      end subroutine montarLmStencilVel_CRS
      !
      !----------------------------------------------------------------------
      !
      subroutine ordenarLMstencil(LMstencilEq,numCoefPorLinha)
      
      implicit none
!
      integer, intent(in)    :: numCoefPorLinha      
      integer, intent(inout) :: LMstencilEq(numCoefPorLinha)
!
      integer :: menorEq, n, nn , tmp

      do n = 1, numCoefPorLinha
        menorEq=n
        do nn = n+1, numCoefPorLinha
           if(LMstencilEq(nn)<LMstencilEq(menorEq)) menorEq = nn
        end do

        if(n == menorEq) cycle
        tmp                  = LMstencilEq(n)
        LMstencilEq(n)       = LMstencilEq(menorEq)
        LMstencilEq(menorEq) = tmp
                
      enddo

      end subroutine ordenarLMstencil
      !
      !**** new *************************************************************
      !
      subroutine montarPonteiroAp_CRS (Ap, LMstencilEq, numCoefPorLinha, neq, nonzeros)

      implicit none

      integer :: Ap(*)
      integer :: LMstencilEq(neq,numCoefPorLinha)
      integer :: numCoefPorLinha, nonzeros, neq
!
      integer :: l,j

      ! Montando Ap

      call montarListaPonteiros(Ap, LMstencilEq,neq,numCoefPorLinha)

      ! Contando os valores nao nulos
      nonzeros=0
      do l=2, neq+1
            nonzeros=nonzeros+(Ap(l)-Ap(l-1))
      end do

       end subroutine montarPonteiroAp_CRS
      !
      !**** new *************************************************************
      !
      subroutine montarPonteiroAi_CRS (Ai, LMstencilEq, numCoefPorLinha, neq, nonzeros)

      implicit none

      integer, intent(out) :: Ai(*)
      integer :: LMstencilEq(neq,numCoefPorLinha)
      integer :: numCoefPorLinha, nonzeros, neq
!
      integer :: p

      call montarListaIndices(Ai, LMstencilEq, neq, nonzeros, numCoefPorLinha)

      end subroutine montarPonteiroAi_CRS
      !
      !**** new *************************************************************
      !
      subroutine montarListaPonteiros(Ap, LMstencilEq, neq, numCoefPorLinha)

      implicit none
      
      integer, intent(out) :: Ap(*)
      integer, intent(in) :: neq, numCoefPorLinha,LMstencilEq(neq,numCoefPorLinha)
      integer :: LMstencilEqTemp(0:numCoefPorLinha)
!
      integer :: n, l

      posPonteiro=0
      contPonteiro=0
      do l=1, neq
         posPonteiro=posPonteiro+1
         Ap(posPonteiro)=contPonteiro+1
         LMstencilEqTemp=0
         LMstencilEqTemp(1:numCoefPorLinha)=LMstencilEq(l,:)

         do n = 1, numCoefPorLinha
            if(LMstencilEqTemp(n) == LMstencilEqTemp(n-1)) cycle
            contPonteiro= contPonteiro+1
         enddo    

         if(l==neq) then
            posPonteiro=posPonteiro+1
            Ap(posPonteiro)=contPonteiro+1
         end if
      end do

      end subroutine montarListaPonteiros
      !
      !**** new *************************************************************
      !
      subroutine montarListaIndices(Ai, LMstencilEq, neq, nonzeros, numCoefPorLinha)
!
      implicit none
      integer, intent(inout) :: Ai(nonzeros)
      integer, intent(in)    :: nonzeros, neq, numCoefPorLinha
      integer, intent(in)    :: LMstencilEq(neq,numCoefPorLinha)
!
      integer :: LMstencilEqTemp(0:numCoefPorLinha)
      integer :: i, n, posColunas
!
      posColunas=0
      do i=1, neq
         LMstencilEqTemp=0
         LMstencilEqTemp(1:numCoefPorLinha)=LMstencilEq(i,:)
         do n = 1, numCoefPorLinha
             if(LMstencilEqTemp(n).ne.LMstencilEqTemp(n-1).and.LMstencilEqTemp(n).ne.0 ) then
                 posColunas=posColunas+1
                 Ai(posColunas)= LMstencilEqTemp(n)
             end if
         enddo
      end do
!
      end subroutine montarListaIndices
!
      end subroutine criarPonteirosMatEsparsa_CRS_MR
      
      !
!=======================================================================
!     
      subroutine criarListaVizinhosCRS(nen,numnp,numel,nVizinMax,conecElem,listaDosElemsPorNo)
      implicit none
!     
!     Objetivo: cria a matriz listaDosElemsPorNo: .... indices dos elementos que possuem o no
!
      integer nen,numnp,numel,nVizinMax
      integer, dimension(nen,numel)  :: conecElem
      integer, dimension(nVizinMax,numnp)  :: listaDosElemsPorNo
      integer :: no,nel,l,numVizAtual
!
      listaDosElemsPorNo(nVizinMax,:) = 0 
      do nel=1, numel
         do l=1, nen
            no=conecElem(l,nel)
                  listaDosElemsPorNo(nVizinMax,no) = listaDosElemsPorNo(nVizinMax,no) + 1 
                  numVizAtual               = listaDosElemsPorNo(nVizinMax,no)  
                  listaDosElemsPorNo(numVizAtual,no) = nel
         end do
      end do
!
!        do no=1, numnp
!           numVizAtual = listaDosElemsPorNo(nVizinMax,no)  
!           write(*,'(8(1X,I3))') no, numVizAtual , listaDosElemsPorNo(1:numVizAtual,no)
!        end do
!           write(*,*) '---------------'
!       stop
!    
      end subroutine


!
     end module mSolverPardiso

