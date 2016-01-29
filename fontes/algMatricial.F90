!==============================================================================
!         programa de elementos finitos  
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
! 
!         + implementacoes de Abimael Loula
!
!         + novo projeto: modular e fortran 90, por
!         Eduardo Garcia,        bidu@lncc.br 
!         Tuane Lopes,           tuane@lncc.br
!
!         LNCC/MCT
!         Petropolis, 07.2013
!==============================================================================
!
      module mAlgMatricial

        integer              :: neqP, nalhsP
        integer              :: neqV, nalhsV
        integer              :: neqD, nalhsD
        integer              :: nedV
        integer              :: meanbwV, meanbwP
        real*8,  allocatable :: alhsP(:,:), brhsP(:)
        real*8,  allocatable :: alhsV(:), brhsV(:)
        real*8,  allocatable :: alhsD(:), clhsD(:), brhsD(:)
        integer, allocatable :: idVeloc(:,:), idDesloc(:,:)
        integer, allocatable :: idiagV(:), idiagD(:)
        integer, allocatable :: lmV(:,:,:), lmD(:,:,:)
        integer :: numCoefPorLinhaVel, numCoefPorLinhaGeo
        integer, allocatable ::listaDosElemsPorNoCRS(:,:)

!funcoes e subrotinas
        public :: backns, factns, back, factor
        public :: diag, load, addnsl, addlhs, addrhs
        public :: btod, kdbc, pivots, ftod, ftodTIME, btdb, colht

        public :: coldot, rowdot
        public :: matadd

      contains
!
!**** new **********************************************************************
!
      subroutine solverDiretoSkyLine(alhs, brhs, idiag, nalhs, neq,label)
!
         implicit none
!
         integer, intent(in)   :: nalhs, neq, idiag(*)
         real*8, intent(inout) :: alhs(*), brhs(*)
         character(len=3) :: label
         real*8 :: t1,t2

         call timing(t1)
         call factor(alhs,idiag,nalhs,neq)
         call timing(t2)
#ifdef mostrarTempos
         print*, "solver skyline: ", label, ", tempo de factorization=", t2-t1
#endif

         call timing(t1)
         call back  (alhs,brhs, idiag,neq)
         call timing(t2)
#ifdef mostrarTempos
         print*, "solver skyline: ", label, ", tempo de backsubstitution=", t2-t1
#endif

      end subroutine solverDiretoSkyLine

!
!----------------------------------------------------------------------
!
      subroutine solverGaussBanda(s,x,r,ns,lb)                                        
      implicit none                       
      
      integer :: ns,lb
      real*8 :: s(ns,lb), r(ns), x(ns)
!    
      integer :: l1,l2,ll,n,i,j,k,l,ls
      real*8 :: c
      real*8 :: part, ti, tf, tfAnt, tempoProc, tempoProcTotEst, somaTempoProcTotEst
      integer :: cont

      write(*,*) "..... Resolvendo o Sistema Pressao 3D, solver interno"
 
      part=0.02
      if(ns > 50000)  part=0.001
!       print* , ' nt(part*ns) =', int(part*ns) 

      call timing(ti) 
      tfAnt=ti
      cont = 0
      somaTempoProcTotEst = 0.0

      ls=(lb-1)/2                                    
      l1=ls+1                                        
      l2=ls+2                                        
      do n=1,ns                                  
!         call estimativasDeDesempenho()
        i=n                                            
        do l=2,l1                                  
          i=i+1                                          
          if(i.gt.ns) cycle ! inicia nova iteracao, proximo l 
          ll=n-i+l1                                      
          if(s(i,ll).ne.0.) then                    
            c=s(i,ll)/s(n,l1)                              
            j=ll-1                                         
            do k=l1,lb                                 
              j=j+1                         
              if(s(n,k).eq.0.) cycle ! inicia nova iteracao, proximo k 
              s(i,j)=s(i,j)-c*s(n,k)                         
            end do                                       
            r(i)=r(i)-c*r(n)           
          end if
        end do                                       
      end do                                       
      r(ns)=r(ns)/s(ns,l1)                                               
      n=ns - 1                                                               
       do n = ns - 1, 1, -1
!        if(mod(n,int(0.05*ns)) == 0) write(*,*)'2, n = ', n
        l=n                                                                
        do k=l2,lb                                                     
          l=l+1                                                              
          if(l.gt.ns) exit ! finaliza o laco de repeticao, k > lb  
          r(n)=r(n)-s(n,k)*r(l)                                              
        end do                                                           
        r(n)=r(n)/s(n,l1)                                                  
       end do
	   x=r
	   
    contains
    subroutine estimativasDeDesempenho()

       if(mod(n,int(part*ns)) == 0) then
             cont = cont + 1
             call timing(tf) 
             tempoProc           = tf - tfAnt
             tempoProcTotEst     = ns * tempoProc/(part*ns)
             somaTempoProcTotEst = somaTempoProcTotEst + tempoProcTotEst 
             write(*,*)' n = ', n
             write(*,*) int(part*ns), 'iteracoes, tempo =', tempoProc
             write(*,*)' tempo total     (inst.) estimado =', tempoProc *      ns  / (part*ns)  
             write(*,*)' tempo decorrido (inst.) estimado =', tempoProc *       n  / (part*ns)  
             write(*,*)' tempo restante  (inst.) estimado =', tempoProc * (ns - n) / (part*ns)
    !         write(*,*)' tempo total (medio)       estimado =', somaTempoProcTotEst/cont  
    !         write(*,*)' tempo restante            estimado =', somaTempoProcTotEst/cont -  (tf - ti)
             !write(*,*)' tempo de estimado total =' , ns * tempoProc/(part*ns)
             !write(*,*)' tempo restante          =' , ns * tempoProc/(part*ns) - (tf - ti)
             tfAnt=tf
       end if
      end subroutine estimativasDeDesempenho
      end  subroutine
!
!**** new **********************************************************************
!
      subroutine diag(idiag,neq,n)
!
      implicit none
!
!.... program to compute diagonal addresses of left-hand-side matrix
!
      integer :: neq, n
      integer :: idiag(neq)
!
      integer :: i

      n = 1
      idiag(1) = 1 
      if (neq.eq.1) return
!
      do 100 i=2,neq
      idiag(i) = idiag(i) + idiag(i-1) + 1
  100 continue
      n = idiag(neq)
!
      return
      end subroutine
!
!**** new **********************************************************************
!
      subroutine load(id,f,brhs,ndof,numnp,nlvect)
!
!.... program to accumulate nodal forces and transfer into
!        right-hand-side vector
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer :: id(ndof,*)
      real*8  :: f(ndof,numnp,*),brhs(*)
      integer :: ndof, numnp, nlvect
!
      integer :: nlv
      integer :: i, j, k
!
      do 300 i=1,ndof
!
      do 200 j=1,numnp
      k = id(i,j)
      if (k.gt.0) then
!
         do 100 nlv=1,nlvect
         brhs(k) = brhs(k) + f(i,j,nlv)
  100    continue
!
      endif
!
  200 continue
!
  300 continue
!
      return
      end subroutine
!
!**** new **********************************************************************
!
      SUBROUTINE LOADTIME(id,f,brhs,ndof,numnp,nlvect,XTIME)
!
!.... program to accumulate nodal forces and transfer into
!        right-hand-side vector
!
      IMPLICIT NONE
!
!.... remove above card for single-precision operation
!
      integer :: id(ndof,*)
      real*8  :: f(ndof,numnp,*),brhs(*)
      integer :: ndof, numnp, nlvect
!
      integer :: nlv
      integer :: i, j, k
      REAL*8  :: XTIME
!
      DO 300 I=1,NDOF
         DO 200 J=1,NUMNP
            K = ID(I,J)
            IF (K.GT.0) THEN
            DO 100 NLV=1,NLVECT
               BRHS(K) = BRHS(K) + F(I,J,NLV)*XTIME
  100       CONTINUE
            ENDIF
  200    CONTINUE
  300 CONTINUE
!
      RETURN
      END SUBROUTINE
!
!**** new **********************************************************************
!
      subroutine load2(id,f,brhs,ndof,numnp,nlvect,nnp,N)
!
!.... program to accumulate nodal forces and transfer into
!        right-hand-side vector
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer :: id(ndof,*)
      real*8  :: f(ndof,numnp,*),brhs(*)
      integer :: ndof, numnp, nlvect,nnp, N
!
      integer :: nlv
      integer :: i, j, k
      real*8 :: teste

       

!       if(mod(mod(nnp,N/5),2)==0) then
      if(mod(nnp,2)==0) then
         print*, "ligado", nnp
         teste=1.0
      else
         print*, "desligado", nnp
         teste=0.0
      endif



!
      do 300 i=1,ndof
!
      do 200 j=1,numnp
      k = id(i,j)
      if (k.gt.0) then
!
         do 100 nlv=1,nlvect
         brhs(k) = brhs(k) + f(i,j,nlv)*teste
  100    continue
!
      endif
!
  200 continue
!
  300 continue
!
      return
      end subroutine
!**** new*********************************************************
      subroutine addnsl(alhs,clhs,eleffm,idiag,lm,nee,ldiag)
!
!         program to add element left-hand-side matrix to          
!                global left-hand-side matrix                      
!                                                                  
!        ldiag = .true.,  add diagonal element matrix              
!                                                                  
!        ldiag = .false, then                                     
!        add full nonsymmetric element matrix                   
!                                                                  
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8  :: alhs(*),clhs(*),eleffm(nee,*)
      integer :: idiag(*),lm(*), nee
      logical ldiag
!
      integer :: i,j,k,l,m
!
      if (ldiag) then
!
         do 100 j=1,nee
             k = iabs(lm(j))
            if (k.gt.0) then
               l = idiag(k)
               alhs(l) = alhs(l) + eleffm(j,j)
            endif
  100    continue
!
      else
!
         do 400 j=1,nee
             k = iabs(lm(j))
            if (k.gt.0) then
               do 200 i=1,nee
                   m = iabs(lm(i))
                  if (m.gt.0) then
                     if (k.gt.m) then
                        l = idiag(k) - k + m
                        alhs(l) = alhs(l) + eleffm(i,j)
                     else
                        l = idiag(m) - m + k
                        clhs(l) = clhs(l) + eleffm(i,j)
                     endif
                     if (k.eq.m) then
                        l = idiag(k)
                        alhs(l) = alhs(l) + eleffm(i,j)
                        clhs(l) = alhs(l)
                     endif
                  endif
  200          continue
            endif
  400    continue
!
      endif
!
      return
      end subroutine

!**** new*********************************************************
      subroutine backns(a,c,b,idiag,neq)
!
!.... program to perform forward reduction and back substitution
!
      implicit none
!
!.... deactivate above card(s) for single-precision operation
!
      real*8  :: a(*),c(*),b(*)
      integer :: idiag(*)
      integer :: neq
!
      integer :: i, j, jj, jcolht, jjlast, jjnext, istart, jtemp
      real*8  :: ajj, bj
!     real*8, external :: coldot
!
!
!.... forward reduction
!
      jj = 0
!
      do 100 j=1,neq
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
      if (jcolht.gt.1) then
           b(j) = b(j) - coldot(c(jjlast+1),b(j-jcolht+1),jcolht-1)
      endif
  100 continue
!
!.... diagonal scaling
!
      do 200 j=1,neq
      ajj = a(idiag(j))
!
!.... warning: diagonal scaling is not performed if ajj equals zero
!
      if (ajj.ne.0.0d0) b(j) = b(j)/ajj
  200 continue
!
!.... back substitution
!
      if (neq.eq.1) return
      jjnext = idiag(neq)
!
      do 400 j=neq,2,-1
      jj     = jjnext
      jjnext = idiag(j-1)
      jcolht = jj - jjnext
      if (jcolht.gt.1) then
         bj = b(j)
         istart = j - jcolht + 1
         jtemp  = jjnext - istart + 1
!
         do 300 i=istart,j-1
         b(i) = b(i) - a(jtemp+i)*bj
  300    continue
!
      endif
!
  400 continue
!
      return
      end subroutine
!**** new *********************************************************
      subroutine factns(a,c,idiag,neq)
!
!.... program to perform crout factorization: a = l * d * u
!
!        a(i):  coefficient matrix stored in compacted column form;
!               after factorization contains d and u
!
!        c(i):  non-symmetric lower triangular coefficient matrix stored in
!                compacted row form; after factorization contains l
! 
!
      implicit none
!
!.... deactivate above card(s) for single-precision operation
!
      real*8  :: a(*),c(*)
      integer :: idiag(*)
      integer :: neq
!
      integer :: i, j, ii, jj, ij, iilast, jjlast
      integer :: istart,  icolht, jcolht, jtemp, jm1
      integer :: length
!     real*8, external :: coldot
!
      jj = 0
!
      do 300 j=1,neq
!
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
!
      if (jcolht.gt.2) then
!
!....... for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j)
!
         istart = j - jcolht + 2
         jm1    = j - 1
         ij     = jjlast + 2
         ii     = idiag(istart-1)
!
         do 100 i=istart,jm1
!
         iilast = ii
         ii     = idiag(i)
         icolht = ii - iilast
         length = min0(icolht-1,i - istart + 1)
         if (length.gt.0)  then
            a(ij) = a(ij) - coldot(a(ij-length),c(ii-length),length)
            c(ij) = c(ij) - coldot(c(ij-length),a(ii-length),length)
         endif
         ij = ij + 1
  100    continue
!
      endif
!
      if (jcolht.ge.2) then
!
!....... for column j and i.le.j-1, replace a(i,j) with u(i,j);
!           replace a(j,j) with d(j,j).
!
         jtemp = j - jj
!
         do 200 ij=jjlast+1,jj-1
!
         ii = idiag(jtemp + ij)
!
!....... warning: the following calculations are skipped 
!                 if a(ii) equals zero
!
         if (a(ii).ne.0.0d0) then
             c(ij) = c(ij)/a(ii)
             a(jj) = a(jj) - c(ij)*a(ij)
             a(ij) = a(ij)/a(ii)
         endif
  200    continue
!
      endif
!
  300 continue
!
      return
      end subroutine

!**** NEW ****************************************************************** 
      subroutine back(a,b,idiag,neq)

! 
!.... program to perform forward reduction and back substitution 
! 
      implicit none
! 
!.... remove above card for single-precision operation 
! 
      real*8  :: a(*)
      integer :: neq
      integer :: idiag(neq)
      real*8  :: b(neq)
!
      integer :: i,j,jcolht,istart,jtemp,jj,jjnext,jjlast
      real*8  :: ajj, bj
      integer :: iniA, iniB, fimA, fimB

!
!.... forward reduction 
! 
      jj = 0
! 

      do j=1,neq
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
      if (jcolht.gt.1) then
           iniA=jjlast+1
           iniB=j-jcolht+1
           fimA=iniA+jcolht-1-1
           fimB=iniB+jcolht-1-1
!            b(j) = b(j) - dot_product(a(iniA:fimA),b(iniB:fimB))
            b(j) = b(j) - coldot(a(jjlast+1),b(j-jcolht+1),jcolht-1)
      end if
      enddo
!
!.... diagonal scaling 
! 
      do j=1,neq
      ajj = a(idiag(j))
      if (ajj.ne.0.0d0) then
             b(j) = b(j)/ajj
      end if
      end do

! 
!.... back substitution 
! 
      if (neq.eq.1) return
      jjnext = idiag(neq)
! 
      do j=neq,2,-1
      jj     = jjnext
      jjnext = idiag(j-1)
      jcolht = jj - jjnext
      if (jcolht.gt.1) then
        bj = b(j)
        istart = j - jcolht + 1
        jtemp  = jjnext - istart + 1
        do i=istart,j-1
           b(i) = b(i) - a(jtemp+i)*bj
        enddo
      endif
! 
      end do
! 
      return
      end subroutine

   
!**** new *************************************************************** 
      subroutine factor(a,idiag,nalhs,neq)
! 
!.... program to perform crout factorization: a = u(transpose) * d * u 
! 
!        a(i):  coefficient matrix stored in compacted column form; 
!               after factorization contains d and u 
! 
      implicit none
! 
!.... remove above card for single-precision operation 
! 
      integer, intent(in) :: neq, nalhs
      real*8, intent(inout)  :: a(nalhs)

      integer, intent(in) :: idiag(neq)
!
      integer :: i, j, jlast, icolht, jcolht, istart, jm1, jtemp
      integer :: ii, ij, jj, jlngth, length, iilast, jjnext
      integer :: jjlast
      real*8  :: ajj, bj, temp 
      integer :: iniA, iniB, fimA, fimB
! 
      jj = 0
      i=0; j=0; jlast=0; icolht=0; jcolht=0; istart=0; jm1=0; jtemp=0;
      ii=0;ij=0;jj=0;jlngth=0;length=0;iilast=0;jjnext=0;jjlast=0;
      ajj=0.0d0;bj=0.0d0;temp=0.0d0;
!
      do 300 j=1,neq
!
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
! 
      if (jcolht.gt.2) then
! 
!....... for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j) 
! 
         istart = j - jcolht + 2
         jm1    = j - 1
         ij     = jjlast + 2
         ii     = idiag(istart-1)
! 
         do 100 i=istart,jm1
! 
         iilast = ii
         ii     = idiag(i)
         icolht = ii - iilast
         jlngth = i - istart + 1
         length = min0(icolht-1,jlngth)
         if (length.gt.0) then
           iniA=ii-length
           iniB=ij-length
           fimA=iniA+length-1
           fimB=iniB+length-1
!            a(ij) = a(ij) - dot_product(a(iniA:fimA),a(iniB:fimB))
           a(ij) = a(ij) - coldot(a(ii-length),a(ij-length),length)

         end if
         ij = ij + 1
  100    continue
!
      endif
! 
      if (jcolht.ge.2) then
! 
!....... for column j and i.le.j-1, replace a(i,j) with u(i,j); 
!           replace a(j,j) with d(j,j). 
! 
         jtemp = j - jj
! 
         do 200 ij=jjlast+1,jj-1
! 
         ii = idiag(jtemp + ij)
         if (a(ii).ne.0.0d0) then
            temp  = a(ij)
            a(ij) = temp/a(ii)
            a(jj) = a(jj) - temp*a(ij)
         endif
  200    continue
! 
      endif
! 
  300 continue

      return
      end subroutine


!**** new **********************************************************************
      subroutine btod(id,d,brhs,ndof,numnp)
!
!.... program to perform transfer from r.h.s. to displacement array
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer :: i, j, k
      integer :: ndof, numnp
      integer, dimension(ndof,numnp) :: id
      real*8,  dimension(ndof,numnp) :: d
      real*8,  dimension(*)          :: brhs
!
      do 200 i=1,ndof
!
         do 100 j=1,numnp
         k = id(i,j)
         if (k.gt.0) then 
             d(i,j) = brhs(k)
         end if
  100    continue
!
  200    continue
!
      return
      end subroutine


!**** new **********************************************************************
      subroutine kdbc(eleffm,elresf,dl,nee)
!
!.... program to adjust load vector for prescribed displacement
!     boundary condition
      use mGlobaisEscalares, only: zero
!
       implicit real*8 (a-h,o-z) 
!
!.... remove above card for single-precision operation
!
      integer :: nee
      real*8  :: eleffm(nee,*),elresf(*),dl(*)
!
      integer :: i,j
      real*8  :: val
!
      do 200 j=1,nee
!
      val=dl(j)
      if(val.eq.zero) go to 200
!
      do 100 i=1,nee
      elresf(i)=elresf(i)-eleffm(i,j)*val
100   continue
!
200   continue
!
      return
      end subroutine

!**** new **********************************************************************
      subroutine pivots(a,idiag,neq,nsq,iecho,*)
!
!.... program to determine the number of zero and negative terms in
!        array d of factorization a = u(transpose) * d * u
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!                
      real*8  :: a(*)                                                       
      integer :: idiag(*)
      integer :: neq, nsq, iecho
!
      integer :: iz, in, n, i
!
      iz = 0
      in = 0
!
      do 100 n=1,neq
      i = idiag(n)
      if (a(i).eq.0.) iz = iz + 1
      if (a(i).lt.0.) in = in + 1
  100 continue
!
      write(iecho,1000) nsq,iz,in
!
      return 1
!
 1000 format(' ',&
     ' zero and/or negative pivots encountered                ', ///5x,&
     ' time sequence number   . . . . . . . . . . . (nsq  ) = ',i10//5x,&
     ' number of zeroes . . . . . . . . . . . . . . . . . . = ',i10//5x,&
     ' number of negatives  . . . . . . . . . . . . . . . . = ',i10//5x)
!
      end subroutine

!**** new **********************************************************************
      subroutine ftod(id,d,f,ndof,numnp,nlvect)
!
!.... program to compute displacement boundary conditions
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer :: id(ndof,*)
      real*8  :: d(ndof,*),f(ndof,numnp,*)
      integer :: ndof, numnp, nlvect
!
      integer :: i, j, k, lv
      real*8  :: val
!
!
      do 300 i=1,ndof
!
            do 200 j=1,numnp
!
            k = id(i,j)
            if (k.gt.0) go to 200
            val = 0.0d0
                  do 100 lv=1,nlvect
                  val = val + f(i,j,lv)
100               continue
!
            d(i,j) = val
!
  200       continue
!
 300  continue
      return
      end subroutine
!
!**** new **********************************************************************
!
      subroutine ftodTIME(id,d,f,ndof,numnp,nlvect,XTIME)
!
!.... program to compute displacement boundary conditions
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer :: id(ndof,*)
      real*8  :: d(ndof,*),f(ndof,numnp,*)
      integer :: ndof, numnp, nlvect
!
      integer :: i, j, k, lv
      real*8  :: val, XTIME
!
      do 300 i=1,ndof
         do 200 j=1,numnp
            k = id(i,j)
            if (k.gt.0) go to 200
            val = 0.0d0
            do 100 lv=1,nlvect
               val = val + f(i,j,lv)*XTIME
100         continue
         d(i,j) = val
200      continue
300   continue
!
      return
!
      end subroutine

!**** new **********************************************************************
      subroutine btdb(elstif,b,db,nee,nrowb,nstr)
!
!.... program to multiply b(transpose) * db taking account of symmetry
!        and accumulate into element stiffness matrix
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8  :: elstif(nee,*),b(nrowb,*),db(nrowb,*)
      integer :: nee,nrowb,nstr
!
      integer :: i,j
!
      do 200 j=1,nee
!
      do 100 i=1,j
      elstif(i,j) = elstif(i,j) + coldot(b(1,i),db(1,j),nstr)   
  100 continue
!
  200 continue
!
      return
      end subroutine
!
!**** new **********************************************************************
!
      subroutine colht(idiag,lm,ned,nen,numel,neq)
!
!.... program to compute column heights in global left-hand-side matrix
!
      implicit none
      integer :: idiag(*),lm(ned,nen,*)
      integer :: ned, nen, numel, neq
!
      integer :: i, j, k
      integer :: m, min, num
!
      do 500 k=1,numel
      min = neq
!
      do 200 j=1,nen
!
      do 100 i=1,ned
      num = lm(i,j,k)
      if (num.gt.0) min = min0(min,num)
  100 continue
!
  200 continue
!
      do 400 j=1,nen
      do 300 i=1,ned
      num = lm(i,j,k)
      if (num.gt.0) then
         m = num - min
         if (m.gt.idiag(num)) then
              idiag(num) = m
         endif
      endif
!
  300 continue

  400 continue
!
  500 continue
!
      return
      end subroutine
! **** new *********************************************************************
       subroutine addlhs(alhs,eleffm,idiag,lm,nee,ldiag,lsym)
! 
! .... program to add element left-hand-side matrix to
!         global left-hand-side matrix
! 
!         ldiag = .true.,  add diagonal element matrix
! 
!         ldiag = .false., then
!           lsym = .true., add upper triangle of full element matrix
!           lsym = .false., add full element matrix
! 
       implicit none
! 
! .... deactivate above card(s) for single-precision operation
! 
       real*8  :: alhs(*),eleffm(nee,*)
       integer :: idiag(*),lm(*)
       integer :: nee
       logical :: ldiag,lsym
!
       integer :: i, j, l, k, m

       if (ldiag) then
! 
          do 100 j=1,nee
          k = lm(j)
          if (k.gt.0) then
             l = idiag(k)
             alhs(l) = alhs(l) + eleffm(j,j)
          endif
  100     continue
! 
       else
! 
          do 400 j=1,nee
          k = lm(j)
          if (k.gt.0) then
! 
             do 200 i=1,j
             m = lm(i)
             if (m.gt.0) then
                if (k.ge.m) then
                   l = idiag(k) - k + m
                else
                   l = idiag(m) - m + k
                endif
                alhs(l) = alhs(l) + eleffm(i,j)

             endif
  200       continue
! 
             if (.not. lsym) then
                do 300 i = j,nee
                m = lm(i)
               if (m .gt. 0) then
                 if (k .ge. m) then
                   l = idiag(k) - k + m
                else
                   l = idiag(m) - m + k
                 endif
               endif
  300         continue
             endif
          endif
  400    continue
! 
       endif
! 
       return
       end subroutine

!      
!:::: new ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
       subroutine addlhsGEO(alhs,clhs,eleffm,idiag,lm,nee,ldiag,lsym) 
! 
! .... program to add element left-hand-side matrix to 
!         global left-hand-side matrix 
! 
!         ldiag = .true.,  add diagonal element matrix 
! 
!         ldiag = .false., then 
!           lsym = .true., add upper triangle of full element matrix 
!           lsym = .false., add full element matrix 
! 
       implicit real*8 (a-h,o-z) 
! 
! .... deactivate above card(s) for single-precision operation 
! 
       integer :: nee
       real*8  :: alhs(*),clhs(*),eleffm(nee,*)
       integer :: idiag(*),lm(*) 
       logical :: ldiag,lsym 
!
       integer :: i,j,k,l,m
! 
       if (ldiag) then 
! 
          do 100 j=1,nee 
          k = lm(j) 
          if (k.gt.0) then 
             l = idiag(k) 
             alhs(l) = alhs(l) + eleffm(j,j) 
          endif 
  100     continue 
! 
       else 
! 
          do 400 j=1,nee 
          k = lm(j) 
          if (k.gt.0) then 
! 
             do 200 i=1,j 
             m = lm(i) 
             if (m.gt.0) then 
                if (k.ge.m) then 
                   l = idiag(k) - k + m 
                else 
                   l = idiag(m) - m + k 
                endif 
                alhs(l) = alhs(l) + eleffm(i,j) 
             endif 
  200       continue 
! 
             if (.not. lsym) then 
                do 300 i = j,nee 
                m = lm(i) 
               if (m .gt. 0) then 
                 if (k .ge. m) then 
                   l = idiag(k) - k + m 
                else 
                   l = idiag(m) - m + k 
                 endif 
                 clhs(l) = clhs(l) + eleffm(i,j) 
               endif 
  300         continue 
             endif 
          endif 
  400    continue 
! 
       endif 
! 
       return 
       end subroutine
!:
!**** new **********************************************************************
      subroutine addrhs (brhs,elresf,lm,nee)
!
!.... program to add element residual-force vector to
!        global right-hand-side vector
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8  :: brhs(*),elresf(*)
      integer :: lm(*)
      integer :: nee
!
      integer :: k, j
!
      do 100 j=1,nee
      k = lm(j)
      if (k.gt.0) brhs(k) = brhs(k) + elresf(j)
  100 continue
!
      return
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine addlhsInterno(l, LMstencilEq, linhamatriz)
        implicit none
!
        integer, intent(in) :: l, LMstencilEq(0:7)
        real*8,  intent(in) :: linhaMatriz(7) 
!
        integer :: n
!
        do n = 1, 7
           alhsP(l,LMstencilEq(n))=alhsP(l, LMstencilEq(n))+linhaMatriz(n)
        end do
!
      end subroutine

!**** new **********************************************************************
      subroutine calcularResiduo(alhs,nalhs,vetX,vetBOriginal,neq,idiag,res)

       implicit none
!
       integer, intent(in)    :: nalhs, neq
       real*8,  intent(in)    :: alhs(nalhs), vetX(neq), vetBOriginal(neq)
       real*8,  intent(inout) :: res(neq)
       integer, intent(in)    :: idiag(neq)
!
       real*8  :: vetBcalculado(neq), norma
       integer :: i

       call matVecMulSkyLine(alhs,nalhs,vetX,neq,idiag,vetBcalculado)
       res=vetBOriginal-vetBcalculado
       norma = 0.0d0
       do i=1, neq
          norma= norma + res(i) * res(i)
       enddo
       norma = sqrt(norma) 
       write(*,'(a,1pe12.4,a)'), "r=Au-b, norma do residuo da solucao, r=" , norma
 
       end subroutine
!**** new **********************************************************************
      subroutine matVecMulSkyLine(alhsV,nalhsV,vetX,neqV,idiagV,vetBcalc)

       implicit none
!
       integer :: nalhsV, neqV
       real*8  :: alhsV(nalhsV), vetX(neqV), vetBcalc(neqV)
       integer :: idiagV(neqV)
!
       integer :: i,j, altura, jcolht, cont
       real*8  :: aij, xj
       integer :: jj, jjnext, istart,jtemp,jjlast
!
       vetBcalc=0.0
       altura=0
       cont=neqV
       jcolht=0
       jj = 0
! 
        do j=1,neqV
        jjlast = jj
        jj     = idiagV(j)
        jcolht = jj - jjlast
        if (jcolht.gt.1) then
              vetBcalc(j) =  vetBcalc(j) + coldot(alhsV(jjlast+1),vetX(j-jcolht+1),jcolht-1)
        end if
        enddo

      do i=1,neqV
          aij=alhsV(idiagV(i))
          xj =vetX(i)
          vetBcalc(i)=vetBcalc(i)+aij*xj
      end do


        jjnext = idiagV(neqV)
        do j=neqV,2,-1
         jj     = jjnext
         jjnext = idiagV(j-1)
         jcolht = jj - jjnext
         if (jcolht.gt.1) then
           xj     = vetX(j)
           istart = j - jcolht + 1
           jtemp  = jjnext - istart + 1
!
!            forma natural, por linhas: a cada novo i tem ui pronto
!                  ui =  ui + aij*xj, j =1, neq, i = 1, neq
! 
!            forma alternativa, por colunas: a cada novo i tem uma
!                                  contribuicao de xi para o valor total ui 
!                  uj =  uj + aji*xi, j =1, neq, i = 1, neq,

             do i=istart,j-1
                vetBcalc(i) = vetBcalc(i) + alhsV(jtemp+i)*xj
             enddo
           endif
! 
        end do
       end subroutine

!**** new **********************************************************************
      function coldot(a,b,n)
!
!.... program to compute the dot product of vectors stored column-wise
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8  :: a(*),b(*)
      integer :: n
!
      real*8  :: coldot
      integer :: i
!

      coldot = 0.0d0
!
!       do 100 i=1,n
!       coldot = coldot + a(i)*b(i)
       coldot=dot_product(a(1:n),b(1:n))

!   100 continue

!
      return
      end function
!**** new **********************************************************************
      subroutine matadd(a,b,c,ma,mb,mc,m,n,iopt)
!
!.... program to add rectangular matrices
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8  :: a(ma,*),b(mb,*),c(mc,*)
      integer :: ma,mb,mc,m,n,iopt
!
      integer :: i,j
!
      go to (1000,2000,3000),iopt
!
!.... iopt = 1, add entire matrices
!
 1000 do 1200 j=1,n
!
      do 1100 i=1,m 
      c(i,j) = a(i,j) + b(i,j)
 1100 continue
!
 1200 continue
      return
!
!.... iopt = 2, add lower triangular and diagonal elements
!
 2000 do 2200 j=1,n
!
      do 2100 i=j,m 
      c(i,j) = a(i,j) + b(i,j)
 2100 continue
!
 2200 continue
      return
!
!.... iopt = 3, add upper triangular and diagonal elements
!
 3000 do 3200 j=1,n
!
      do 3100 i=1,j 
      c(i,j) = a(i,j) + b(i,j)
 3100 continue
!
 3200 continue
      return
!
      end subroutine

!**** new **********************************************************************
      function rowdot(a,b,ma,mb,n)
!
!.... program to compute the dot product of vectors stored row-wise
!
      implicit none
!
!.... remove above card for single precision operation
!
      real*8  :: a(ma,*),b(mb,*)
      integer :: ma, mb, n
!
      real*8  :: rowdot
      integer :: i
!
      rowdot = 0.0d00
!
      do i=1,n
      rowdot = rowdot + a(1,i)*b(1,i)
      enddo
!
      return
      end function

  end module


!
!====NEW MODULE=========================================================
!  

  module mSolucoesExternas
#ifdef withHYPRE
      include 'mpif.h'
#endif
        !umfpack
        integer, allocatable :: AiVel(:), ApVel(:), AiGeo(:), ApGeo(:)
        integer, allocatable :: LMstencilEqVel(:,:), LMstencilEqGeo(:,:)

        integer, allocatable :: listaDosElemsPorNoCRS(:,:)

        INTEGER ptV(64), iparmV(64), ptG(64), iparmG(64)
        REAL*8  dparmV(64), dparmG(64)

        integer :: posPonteiro, contPonteiro, posColunas, posCoef
        integer :: lda, nonzeros, nonzerosEst
        logical :: primeiravezVel, primeiravezGeo

      contains

!
!=======================================================================
!  
#ifdef withlapack
      subroutine solverLapackPPD(Ax, x, b, lda,nbt,numel)
!
      implicit none
!
      real*8        :: Ax(lda,numel), x(1,numel), b(numel)
      integer       :: lda, nbt, numel 
      integer       :: ml, mu, m, n, i
      integer       :: info, ipiv(numel)
!  
      write(*,*) "..... Resolvendo o Sistema de Equacoes 3D, LAPACK"
! 
      ml = (nbt - 1)/2
      mu = ml
      m = ml + mu + 1
      n = numel
      lda = 2 * ml + mu + 1
!  
      write ( *, '(a,4i6)' ) '  Bandwidth is ', m, ml, mu, n 

      !  Factor the matrix.
      call dgbtrf ( n, n, ml, mu, Ax, lda, ipiv, info )

      if ( info .ne. 0 ) then
         write ( *, '(a,i8)' ) '  Factorization failed, INFO = ', info
         return
      end if

      !  Solve the linear system.
      call dgbtrs ( 'n', n, ml, mu, 1, Ax, lda, ipiv, b, n, info )
!      
      x(1,:)=b(:)

      end subroutine solverLapackPPD
#endif


!
!=======================================================================
!  
#ifdef withpardiso

      subroutine solverPardisoPPD_Vel(ia, ja, a, x, b, neq, nonzeros, simetria, label, parte)

      implicit none
!
        INTEGER, intent(in)  ::  ia(neq+1), ja(nonzeros)
        REAL*8, INTENT(IN)   :: a(nonzeros)
        REAL*8, INTENT(INOUT):: b(neq)
        REAL*8, INTENT(OUT)  :: x(neq)
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
!         INTEGER pt(64), iparm(64)
!         REAL*8  dparm(64)
        

!
!  .. Setup Pardiso control parameters und initialize the solvers     
!     internal adress pointers. This is only necessary for the FIRST   
!     call of the PARDISO solver.                                     
!     
! !$OMP PARALLEL


      if(parte=='fact'.or. parte=='full') then

      iparmV=0
      iparmV(0)=0
      mtype     = -2   ! real and symmetric matrix, inefinite
      iparmV(11) = 0    !Do not use (symmetric matrices).
      if(simetria.eqv..false.)  mtype     = 11   ! unsymmetric matrix, indefinite

      solver    = 0    ! use 0 for sparse direct method or 1 multi-recursive iterative solver
      msglvl    = 0    ! with statistical information
      iparmV(33) = 0    ! compute determinant 
      iparmV(52) = 1    !For OpenMP-threaded solver
      iparmV(2)  = 0 !ou 0?      !Fill-In reduction reordering.
      iparmV(27) = 1
      
      nrhs      = 1
      mnum      = 1
      ptV        = 0
      idum      = 0
      ddum      = 0.0
      dparmV     = 0.0
      maxfct    = 1
      error     = 0
!
!  .. Numbers of Processors ( value of OMP_NUM_THREADS )
!
       iparmV(3) = 1
! 
! #ifdef withOMP
!        iparm(3) = omp_get_num_threads()
! #endif

! !$OMP END PARALLEL
!        print*, "em pardiso com", iparm(3), "threads"
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
!        WRITE(*,*) '[PARDISO]: License check was successful ... '
      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!..   Reordering and Symbolic Factorization, This step also allocates
!     all memory that is necessary for the factorization
!
!         WRITE(*,*) 'Begining reordering ... ', label
         call timing(t1)
         phase     = 11     ! only reordering and symbolic factorization

         CALL pardiso (ptV, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                    idum, nrhs, iparmV, msglvl, ddum, ddum, error, dparmV)
         call timing(t2)
! #ifdef mostrarTempos
!         write(*,*) "reordering: ", label, ", tempo de parede = ", t2 - t1
! #endif

      IF (error .NE. 0) THEN
        WRITE(*,*) 'The following ERROR was detected: ', error
        STOP
      END IF

!       WRITE(*,*) 'Number of nonzeros in factors   = ',iparm(18)
!       WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

!
!.. Factorization.
!
!      WRITE(*,*) 'Begining factorization  ... '
       call timing(t1)
      phase     = 22  ! only factorization
      CALL pardiso (ptV, maxfct, mnum, mtype, phase, neq, a, ia, ja,  &
                    idum, nrhs, iparmV, msglvl, ddum, ddum, error, dparmV) 
      call timing(t2)
! #ifdef mostrarTempos
!      write(*,*) "factorization: ", label, ", tempo = ", t2 - t1
! #endif

    IF (error .NE. 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', error
      STOP
    ENDIF 


!       if(label=='geo') then
!         ptG=pt
!         iparmG=iparm
!         dparmG=dparm
!       else
!         ptV=pt
!         iparmV=iparm
!         dparmV=dparm
!       endif


     endif  !     if(parte=='fact'.or. parte=='full') then


      if(parte=='back'.or. parte=='full') then

!       if(label=='geo') then
!         pt=ptG
!         iparm=iparmG
!         dparm=dparmG
!       else
!         pt=ptV
!         iparm=iparmV
!         dparm=dparmV
!       endif


!.. Back substitution and iterative refinement
!      WRITE(*,*) 'Begining backsubstitution  ... '
       call timing(t1)
      iparmV(8)  = 1   ! max numbers of iterative refinement steps
      phase     = 33  ! only solve
      CALL pardiso (ptV, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                   idum, nrhs, iparmV, msglvl, b, x, error, dparmV) 
       call timing(t2)
! #ifdef mostrarTempos
!      write(*,*) "backsubstitution: ",label, ", tempo = ", t2 - t1
! #endif

       call timing(tt2)
! #ifdef mostrarTempos
!      WRITE(*,*) 'Solve completed ...  ',label, ", tempo =", tt2-tt1
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
!      WRITE(*,*) 'Begining memory desalocation  ... '
      call timing(t1)
      phase = -1
      CALL pardiso (ptV, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                   idum, nrhs, iparmV, msglvl, b, x, error, dparmV)
      call timing(t2)
#ifdef mostrarTempos
!      write(*,*) "memory desalocation: ", label, ", tempo = ", t2 - t1
#endif

      b(1:neq) = x(1:neq)

      endif !if(parte=='back'.or. parte=='full') then

      end subroutine solverPardisoPPD_Vel

#endif


!
!=======================================================================
!  
#ifdef withpardiso

      subroutine solverPardisoPPD_Geo(ia, ja, a, x, b, neq, nonzeros, simetria, label, parte)

!        use mPropGeoFisica, only: sinj

      implicit none
!
        INTEGER, intent(in)  :: ia(neq+1), ja(nonzeros)
        REAL*8, INTENT(IN)   :: a(nonzeros)
        REAL*8, INTENT(INOUT):: b(neq)
        REAL*8, INTENT(OUT)  :: x(neq)
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

        integer :: ptAux(64), iparmAux(64)
! 
!
!  .. Setup Pardiso control parameters und initialize the solvers     
!     internal adress pointers. This is only necessary for the FIRST   
!     call of the PARDISO solver.                                     
!     
! !$OMP PARALLEL


      if(parte=='fact'.or. parte=='full') then

      iparmG=0
      iparmG(0)=0
      mtype     = -2   ! real and symmetric matrix, inefinite
      iparmG(11) = 0    !Do not use (symmetric matrices).
      if(simetria.eqv..false.)  mtype     = 11   ! unsymmetric matrix, indefinite

      solver    = 0    ! use 0 for sparse direct method or 1 multi-recursive iterative solver
      msglvl    = 0    ! with statistical information
      iparmG(33) = 0    ! compute determinant 
      iparmG(52) = 1    !For OpenMP-threaded solver
      iparmG(2)  = 0 !ou 0?      !Fill-In reduction reordering.
      iparmG(27) = 1
      
      nrhs      = 1
      mnum      = 1
      ptG       = 0
      idum      = 0
      ddum      = 0.0
      dparmG    = 0.0
      maxfct    = 1
      error     = 0

!
!  .. Numbers of Processors ( value of OMP_NUM_THREADS )
!
       iparmG(3) = 1
! 
! #ifdef withOMP
!        iparm(3) = omp_get_num_threads()
! #endif

! !$OMP END PARALLEL
!        print*, "em pardiso com", iparm(3), "threads"
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
!        WRITE(*,*) '[PARDISO]: License check was successful ... '
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!..   Reordering and Symbolic Factorization, This step also allocates
!     all memory that is necessary for the factorization
!
!         WRITE(*,*) 'Begining reordering ... ', label
         call timing(t1)
         phase     = 11     ! only reordering and symbolic factorization

!       print*, "aqui e", sinj
!        write(4125,*) idum, nrhs, iparmG, msglvl, ddum, ddum, error, dparmG
         CALL pardiso (ptG, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                    idum, nrhs, iparmG, msglvl, ddum, ddum, error, dparmG)
         call timing(t2)
!       print*, "aqui f", sinj

!        stop

! #ifdef mostrarTempos
!         write(*,*) "reordering: ", label, ", tempo de parede = ", t2 - t1
! #endif

      IF (error .NE. 0) THEN
        WRITE(*,*) 'The following ERROR was detected: ', error
        STOP
      END IF

!       WRITE(*,*) 'Number of nonzeros in factors   = ',iparm(18)
!       WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

!
!.. Factorization.
!
!      WRITE(*,*) 'Begining factorization  ... '
       call timing(t1)
      phase     = 22  ! only factorization
      CALL pardiso (ptG, maxfct, mnum, mtype, phase, neq, a, ia, ja,  &
                    idum, nrhs, iparmG, msglvl, ddum, ddum, error, dparmG) 
      call timing(t2)
! #ifdef mostrarTempos
!      write(*,*) "factorization: ", label, ", tempo = ", t2 - t1
! #endif

    IF (error .NE. 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', error
      STOP
    ENDIF 

      endif  !     if(parte=='fact'.or. parte=='full') then
! 
      if(parte=='back'.or. parte=='full') then

!.. Back substitution and iterative refinement
!      WRITE(*,*) 'Begining backsubstitution  ... '
       call timing(t1)
      iparmG(8)  = 1   ! max numbers of iterative refinement steps
      phase     = 33  ! only solve
      CALL pardiso (ptG, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                   idum, nrhs, iparmG, msglvl, b, x, error, dparmG) 
       call timing(t2)
! #ifdef mostrarTempos
!      write(*,*) "backsubstitution: ",label, ", tempo = ", t2 - t1
! #endif

       call timing(tt2)
! #ifdef mostrarTempos
!      WRITE(*,*) 'Solve completed ...  ',label, ", tempo =", tt2-tt1
! #endif

      b(1:neq) = x(1:neq)

      endif !if(parte=='back'.or. parte=='full') then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      if(parte=='full') then

!       desalocando memoria
!      WRITE(*,*) 'Begining memory desalocation  ... '
      call timing(t1)
      phase = -1

      CALL pardiso (ptG, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                     idum, nrhs, iparmG, msglvl, b, x, error, dparmG)
         
      call timing(t2)
#ifdef mostrarTempos
!      write(*,*) "memory desalocation: ", label, ", tempo = ", t2 - t1
#endif

      endif !if(parte=='full') then


      end subroutine solverPardisoPPD_Geo

#endif
!
!=======================================================================
!  
#ifdef withumfpack

      subroutine solverUMFPackPPD(Ap, Ai, Ax, x, b,  neq, nonzeros)

      implicit none

      real*8,   intent(inout) :: Ax(*)
      integer,  intent(inout) :: Ap(*), Ai(*)
      real*8,   intent(inout) :: x(*)
      real*8,   intent(in)    :: b(*)
      integer, intent(in) :: neq, nonzeros
!
      real*8  :: control (20), info (90)
      integer :: j, p, symbolic, numeric, sys

!     ----------------------------------------------------------------
!     convert from 1-based to 0-based
!     ----------------------------------------------------------------

      do 60 j = 1, neq+1
          Ap (j) = Ap (j) - 1
60    continue
      do 70 p = 1, nonzeros
          Ai (p) = Ai (p) - 1
70    continue

!     ----------------------------------------------------------------
!     factor the matrix and save to a file
!     ----------------------------------------------------------------

!     set default parameters
      call umf4def (control)

!     print control parameters.  set control (1) to 1 to print
!     error messages only
      control (1) = 1
      call umf4pcon (control)

!     pre-order and symbolic analysis
      call umf4sym (neq, neq, Ap, Ai, Ax, symbolic, control, info)

  !   check umf4sym error condition
  !   if (info (1) .lt. 0) then
  !       print *, 'Error occurred in umf4sym: ', info (1)
  !       stop
  !   endif

!     numeric factorization
      call umf4num (Ap, Ai, Ax, symbolic, numeric, control, info)

!     print statistics for the numeric factorization
!     call umf4pinf (control, info) could also be done.
!     print 90, info (1), info (66), &
!         (info (41) * info (4)) / 2**20, &
!         (info (42) * info (4)) / 2**20, &
!         info (43), info (44), info (45)
90    format ('numeric factorization:',/, &
            '   status:  ', f5.0, /, &
            '   time:    ', e10.2, /, &
            '   actual numeric LU statistics:', /, &
            '   size of LU:    ', f10.2, ' (MB)', /, &
            '   memory needed: ', f10.2, ' (MB)', /, &
            '   flop count:    ', e10.2, / &
            '   nnz (L):       ', f10.0, / &
            '   nnz (U):       ', f10.0)

!     check umf4num error condition
      if (info (1) .lt. 0) then
          print *, 'Error occurred in umf4num: ', info (1)
          stop
      endif
!     solve Ax=b, without iterative refinement
      sys = 0
      call umf4sol (sys, x, b, numeric, control, info)
      if (info (1) .lt. 0) then
          print *, 'Error occurred in umf4sol: ', info (1)
          stop
      endif

!     print the residual.  x (i) should be 1 + i/n
!     call resid (neq, nonzeros, Ap, Ai, Ax, pressaoElem, brhsP, r)
  
!     ----------------------------------------------------------------
!     solve Ax=b, with iterative refinement
!     ----------------------------------------------------------------

      sys = 0
      call umf4solr (sys, Ap, Ai, Ax, x, b, numeric, control, info)
      if (info (1) .lt. 0) then
          print *, 'Error occurred in umf4solr: ', info (1)
          stop
      endif

!     print the residual.  x (i) should be 7 - i/n
!     call resid (neq, nonzeros, Ap, Ai, Ax, x, b, r)

!     free the numeric factorization
      call umf4fnum (numeric)
!     free the symbolic analysis
      call umf4fsym (symbolic)

      do j = 1, neq+1
          Ap (j) = Ap (j) + 1
      enddo
      do p = 1, nonzeros
          Ai (p) = Ai (p) + 1
      enddo
!
      end subroutine solverUMFPackPPD

#endif



!############# RAVIART-THOMAS ###############################


!
!=======================================================================
!    
      subroutine criarPonteirosMatEsparsa_CRS(nsd, ndof, neq, numCoefPorLinha,  &
              conectsElem, listaDosElems, id, numConexoes, nen, numConexoesPorElem,  &
              nonzeros, simetria)

      use mGlobaisEscalares, only: optSolver

      implicit none 
!
      integer, intent(in) :: nsd, neq, numConexoes, nen, numConexoesPorElem,nonzeros,ndof
      logical, intent(in) :: simetria
      integer, intent(out) :: numCoefPorLinha
      integer :: conectsElem(nen,*), listaDosElems(numConexoesPorElem,*)
      integer :: id(ndof,*)
!
       if(numCoefPorLinha==26.or.numCoefPorLinha==18) then
!
         if(.not.allocated(LMstencilEqGeo))allocate(LMstencilEqGeo(neq,numCoefPorLinha)); 
         LMstencilEqGeo=0

         call montarLmStencilGeo_CRS    (LMstencilEqGeo,listaDosElems, id, &
                   conectsElem, numCoefPorLinha, ndof, numConexoes, nen, numConexoesPorElem, neq, simetria)

         allocate(ApGeo(neq+1));    ApGeo=0  
         call montarPonteiroAp_CRS(ApGeo, LMstencilEqGeo, numCoefPorLinha, neq, nonzeros)
!

         allocate(AiGeo(nonzeros)); AiGeo=0
         call montarPonteiroAi_CRS(AiGeo, LMstencilEqGeo, numCoefPorLinha, neq, nonzeros)
!
      else
!
         if(.not.allocated(LMstencilEqVel))allocate(LMstencilEqVel(neq,numCoefPorLinha))
         LMstencilEqVel=0

         call montarLmStencilVel_CRS    (LMstencilEqVel,listaDosElems, id, &
                   conectsElem, numCoefPorLinha, ndof, numConexoes, numConexoesPorElem, neq, simetria)
!
         allocate(ApVel(neq+1));    ApVel=0  
         call montarPonteiroAp_CRS(ApVel, LMstencilEqVel, numCoefPorLinha, neq, nonzeros)
!
         allocate(AiVel(nonzeros)); AiVel=0
         call montarPonteiroAi_CRS(AiVel, LMstencilEqVel, numCoefPorLinha, neq, nonzeros)
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
!=======================================================================
!    
      subroutine solverUMFPack(alhs, b, Ap, Ai, neq, nonzeros)

      implicit none 
!
      real*8  :: alhs(*), b(neq)
      integer :: Ap(*), Ai(*)
      integer, intent(in) :: neq, nonzeros
! 
      real*8, allocatable  :: x(:)
!
      allocate(x(neq)); x=0.0d0

#ifdef withumfpack
       call solverUMFPackPPD(Ap, Ai, alhs, x, b, neq, nonzeros)
#endif

      b = x
      deallocate(x)

      end subroutine solverUMFPack

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
      subroutine solverPardiso(alhs, b, Ap, Ai, neq, nonzeros, simetria, label,parte)

      implicit none 
!
      real*8  :: alhs(*), b(neq)
      integer :: Ap(*), Ai(*)
      integer, intent(in) :: neq, nonzeros
      logical, intent(in) :: simetria
      character(LEN=3) :: label
      character(LEN=4) :: parte
! 
      real*8, allocatable  :: xGeo(:), xVel(:)

#ifdef withpardiso
        if(label=='geo') then
        if(.not.allocated(xGeo)) then
           allocate(xGeo(neq))
           xGeo=0.d0
        endif
        call solverPardisoPPD_Geo(Ap, Ai, alhs, xGeo, b, neq, nonzeros, simetria, label,parte)
        else
        if(.not.allocated(xVel)) then
           allocate(xVel(neq))
           xVel=0.d0
        endif
        call solverPardisoPPD_Vel(Ap, Ai, alhs, xVel, b, neq, nonzeros, simetria, label,parte)
        endif
#endif

!       deallocate(x)

      end subroutine solverPardiso
!
!=======================================================================
!    
        subroutine solverLapack_RT(alhs, b, neq, idiag, numCoefPorLinha, &
              conecLadaisElem, listaDosElemsPorFace, numLados, numLadosElem, &
               nx, ny, nz)

        use mAlgMatricial,     only : lmV
!
        implicit none 
!
        real*8,  intent(in)    :: alhs(*)
        real*8,  intent(inout) :: b(neq)
        integer, intent(in)    :: idiag(*), neq, numCoefPorLinha, numLados, numLadosElem
        integer :: conecLadaisElem(numLadosElem,*), listaDosElemsPorFace(numlados,*)
! 
        real*8, allocatable  :: x(:)
        integer :: nx,ny,nz,nonzeros
!
        integer :: LMstencilEq(neq,numCoefPorLinha)

        real*8, allocatable :: novoAlhs(:,:)
        integer :: nbt, lda
!
        allocate(x(neq))

        call montarLmStencilLapack    (LMstencilEq,listaDosElemsPorFace, &
                         conecLadaisElem, numCoefPorLinha, numLados, numLadosElem)

        call montarSistemaDeEquacoesLapack(LMstencilEq, numCoefPorLinha, neq, nonzeros)

#ifdef withlapack
        call solverLapackPPD(novoAlhs, x, b, lda,nbt,neq)
#endif

        deallocate(novoAlhs)
        deallocate(x)


        contains

      !
      !**** new *************************************************************
      !
      subroutine montarLmStencilLapack(LMstencilEq, listaDosElemsPorFace, conecLadaisElem, numCoefPorLinha, numlados, numladosElem) 
!  
      use mAlgMatricial, only: idVeloc,neqV   
!
      implicit none
!
      integer, intent(in)    :: numCoefPorLinha, numlados,numladosElem
      integer, intent(inout) :: LMstencilEq(neqV,numCoefPorLinha)
      integer, intent(in)    :: listaDosElemsPorFace(numLadosElem,numLados), conecLadaisElem(numLadosElem,*)
!
      integer :: lado, cont, nel, i, j, k, l, numEq
      integer :: LMstencilEqTemp(numCoefPorLinha)
      integer :: contLados, contL
      logical :: propriaEq
!
      numEq=0
      LMstencilEq=0

      do lado=1, numlados
         cont=1
         contLados=0
         propriaEq=.false.
   
         if(idVeloc(1,lado).ne.0) then
            numEq=numEq+1

            do k=1, numladosElem
               if(listaDosElemsPorFace(k,lado).ne.0) then
                  nel=listaDosElemsPorFace(k,lado)
                  do i=1, numLadosElem
                     if (idVeloc(1,conecLadaisElem(i,nel))==0) then
                        LMstencilEq(numEq,cont)=0
                     else
                        if(idVeloc(1,conecLadaisElem(i,nel)).ne.numEq) then
                           LMstencilEq(numEq,cont)=idVeloc(1,conecLadaisElem(i,nel))
                           cont=cont+1
                        end if
                     end if
                  end do
                  LMstencilEq(numEq,numCoefPorLinha)=numEq
               end if
            end do
         end if
      end do


      do l=1, neqV

         LMstencilEqTemp=0
         cont=1

         LMstencilEqTemp(cont)=l
         do i = 1, numCoefPorLinha
          
           if( LMstencilEq(l,i).ne.0 )then                                !(if( LMstencilEq(l,i).gt.l )then !desta forma para considerar a simetria)
                LMstencilEqTemp(cont)=LMstencilEq(l,i)
                cont=cont+1
            end if
         end do

         LMstencilEq(l,:)=LMstencilEqTemp
         call ordenarLMstencilVelocidade(LMstencilEq(l,:),numCoefPorLinha)


      end do
       
     end subroutine montarLmStencilLapack

      !
      !----------------------------------------------------------------------
      !
      subroutine ordenarLMstencilVelocidade (LMstencilEq,numCoefPorLinha)
      
      implicit none
!
      integer, intent(in)    :: numCoefPorLinha      
      integer, intent(inout) :: LMstencilEq(numCoefPorLinha)
!
      integer :: menorEq, n, nn , tmp, cont

      cont=1
      do n = 1, numCoefPorLinha
        menorEq=n
        do nn = n+1, numCoefPorLinha
           if(LMstencilEq(nn)<LMstencilEq(menorEq)) menorEq = nn
        end do

        if(n == menorEq) cycle
        tmp                  = LMstencilEq(n)
        LMstencilEq(n)       = LMstencilEq(menorEq)
        LMstencilEq(menorEq) = tmp
        cont=cont+1
      enddo

      end subroutine ordenarLMstencilVelocidade


      !
      !**** new *************************************************************
      !
      subroutine montarSistemaDeEquacoesLapack (LMstencilEq, numCoefPorLinha, neq, nonzeros)


      implicit none
      integer :: LMstencilEq(neq,numCoefPorLinha)
      integer :: numCoefPorLinha, nonzeros, neq
!
      integer :: l,i

      ! Montando Ap
      if(.not.allocated(ApVel))    allocate(ApVel(neq+1))
      call montarListaPonteirosVelocidade(ApVel, LMstencilEq,neq,numCoefPorLinha)
 
      ! Contando os valores nao nulos
      nonzeros=0
      do l=2, neq+1
            nonzeros=nonzeros+(ApVel(l)-ApVel(l-1))
      end do

      ! Montando Ai
      if(.not.allocated(AiVel)) allocate(AiVel(nonzeros))
      if(allocated(AiVel))      AiVel=0.d0
      call montarListaIndicesVelocidade(AiVel, LMstencilEq, neq, nonzeros, numCoefPorLinha)

      call reordenarLMStencil(LMstencilEq, neq, numCoefPorLinha)

      call calcularBanda(LMstencilEq, neq, numCoefPorLinha, nbt, lda)
      print*, "nbt=", nbt, "lda=", lda

       
      allocate(novoAlhs(lda,neq)); novoAlhs=0.0d0
! 
      call montarMatrizLapack(novoAlhs, alhs, ApVel, AiVel, LMstencilEq, neq, numCoefPorLinha, nbt, lda)


      end subroutine montarSistemaDeEquacoesLapack

      !
      !**** new *************************************************************
      !
      subroutine montarListaPonteirosVelocidade (Ap, LMstencilEq, neqV, numCoefPorLinha)

      implicit none
      
      integer, intent(out) :: Ap(*)
      integer, intent(in) :: neqV, numCoefPorLinha,LMstencilEq(neqV,numCoefPorLinha)
      integer :: LMstencilEqTemp(0:numCoefPorLinha)
!
      integer :: n, l

      do l=1, neqV
         posPonteiro=posPonteiro+1
         Ap(posPonteiro)=contPonteiro+1
         LMstencilEqTemp=0
         LMstencilEqTemp(1:numCoefPorLinha)=LMstencilEq(l,:)
      
         do n = 1, numCoefPorLinha
            if(LMstencilEqTemp(n) == LMstencilEqTemp(n-1)) cycle
            contPonteiro= contPonteiro+1
         enddo    

         if(l==neqV) then
            posPonteiro=posPonteiro+1
            Ap(posPonteiro)=contPonteiro+1
         end if
      end do

      end subroutine montarListaPonteirosVelocidade
      !
      !**** new *************************************************************
      !
      subroutine montarListaIndicesVelocidade(Ai, LMstencilEq, neqV, nonzeros, numCoefPorLinha)
!
      implicit none
      integer, intent(inout) :: Ai(nonzeros)
      integer, intent(in)    :: nonzeros, neqV, numCoefPorLinha
      integer, intent(in)    :: LMstencilEq(neqV,numCoefPorLinha)
!
      integer :: LMstencilEqTemp(0:numCoefPorLinha)
      integer :: i, n, posColunas
!
      posColunas=0
      do i=1, neqV
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
      end subroutine montarListaIndicesVelocidade

      !
      !----------------------------------------------------------------------
      !
      subroutine reordenarLMStencil(LMstencilEq, neq, numCoefPorLinha)
!
      implicit none
!
      integer, intent(in) :: numCoefPorLinha, neq
      integer, intent(inout) :: LMstencilEq(neq,numCoefPorLinha)
!
      integer :: LMstencilEqTemp(numCoefPorLinha)
      integer :: numEq, contL, j
!        
        do numEq=1, neq
        contL=1
        LMstencilEqTemp=0
        do j=1, numCoefPorLinha
              if(LMstencilEq(numEq,j).ne.0) then
                LMstencilEqTemp(contL)=LMstencilEq(numEq,j)
                contL=contL+1
              endif
        end do

        LMstencilEq(numEq,:)=LMstencilEqTemp
        end do

      end subroutine reordenarLMStencil

      !
      !----------------------------------------------------------------------
      !
       subroutine calcularBanda(LMstencilEq, neq, numCoefPorLinha, nbt, lda)
       
       implicit none

       integer, intent(in)    :: LMstencilEq(neq,numCoefPorLinha), neq, numCoefPorLinha
       integer, intent(inout) :: nbt, lda
!
       integer :: elementMeio, numEqMeio, maior, menor
       integer :: ml, mu, i, meiaBanda
!

       elementMeio=(nx/2+1)+(ny/2+1)+(nz/2+1)
       numEqMeio=maxval(lmv(1,1:numLadosElem,elementMeio))

       meiaBanda = 0
       do i= 1, neq
           maior = maxval(LMstencilEq(i, :))
           if(maior-i > meiaBanda) meiaBanda=maior-i
       end do
       write(*,*) " meiaBanda = ",  meiaBanda
       !meiaBanda = 1600
       !write(*,*) " meiaBanda = ",  meiaBanda

       nbt= meiaBanda*2+1
       !nbt= meiaBanda+(meiaBanda/4)+1 !teste para melhora de desempenho.
       
       ml = (nbt - 1)/2
       mu = ml
       lda = 2 * ml + mu + 1

       end subroutine calcularBanda

      !
      !----------------------------------------------------------------------
      !
       subroutine montarMatrizLapack(novoAlhs, alhs, Ap, Ai, LMstencilEq, neq, numCoefPorLinha, nbt, lda)

       implicit none

       real*8, intent(inout) :: novoAlhs(lda,neq)
       real*8, intent(in) :: alhs(*)
       integer, intent(in) :: LMstencilEq(neq,numCoefPorLinha)
       integer :: Ap(*), Ai(*)
       integer :: neq, numCoefPorLinha, nbt, lda
!
       integer :: cont, inicio, fim, linha, i, j, n
       integer :: maiorLM, posicao, coluna
       real*8 :: linhaMatriz(numCoefPorLinha)
       

       do i=1, neq
          cont=1
          inicio = ap(i)
          fim    = ap(i+1)-1
          linha=i
          linhaMatriz=0.0

         do j=inicio, fim

            coluna = ai(j)
            posicao=idiag(coluna)-(ai(j)-linha)
            if(ai(j)<i) posicao=idiag(i)-(linha-ai(j))
            linhaMatriz(cont)=alhs(posicao)
            cont=cont+1
         end do

         maiorLM=maxval(LMstencilEq(i,:),1)
         do n = 1, numCoefPorLinha
            if(LMstencilEq(i,n).ne.0) then
              linha = nbt - i +  LMstencilEq(i,n) 
              novoAlhs(linha,i) = novoAlhs(linha,i)+linhaMatriz(n)
            end if
         end do
      end do

      end subroutine montarMatrizLapack

      end subroutine solverLapack_RT
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
!=======================================================================
! 
#ifdef withHYPRE

      subroutine inicializarMPI(myid_, num_procs_, mpi_comm_)

      integer, intent(out) :: myid_, num_procs_
      integer (KIND=8)     :: mpi_comm_

      integer   ::  ierr  

      print *, " em subroutine inicializarMPI(myid_, num_procs_, mpi_comm_)"
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid_,      ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs_, ierr)

!   Convert the Fortran MPI communicator to a C version that can be used in hypre.
!   Uncomment the line below if you are using MPICH
      mpi_comm_ = MPI_COMM_WORLD

!   Uncomment the line below if you are using LAM-MPI
!     call HYPRE_MPI_Comm_f2c(mpi_comm, MPI_COMM_WORLD, ierr)
      print *, " em subroutine inicializarMPI(myid_, num_procs_, mpi_comm_)"
 
      end subroutine inicializarMPI

      subroutine gerarCoefSistemaAlg(alhs_, Ap_, Ai_, rhs_, x_, rows_, &
                                 neq_, nonZerosT_, Clower_, Cupper_) 

      implicit none
      double precision, allocatable :: alhs_(:)  
      double precision, allocatable :: x_(:), rhs_(:)
      integer,          allocatable :: Ap_(:), Ai_(:), rows_(:)  
      integer, intent(out)          :: neq_, nonzerosT_, Clower_, Cupper_


      integer                      :: Flower, Fupper
      integer                      ::    nnz, i, j, k
      integer                      :: ClowerLido, CupperLido
      character(LEN=*) , parameter :: nomeArquivoLHS        ="LHS.txt";
      character(LEN=*) , parameter :: nomeArquivoLHSstruct  ="LHSstruct.txt";
      character(LEN=*) , parameter :: nomeArquivoRHS        ="RHS.txt";
      character(LEN=*) , parameter :: nomeArquivoIG         ="initialGuess.txt";
      integer,           parameter :: arqHandleLHS=50, arqHandleLHSstruct=51
      integer,           parameter :: arqHandleRHS=52,        arqHandleIG=53
      integer                      :: valorDesprezadoI, valorDesprezadoF
      integer    local_size

!     Each processor knows only of its own rows - the range is denoted by lower
!     and upper.  Here we partition the rows. We account for the fact that
!     N may not divide evenly by the number of processors.
       open ( unit=arqHandleLHSstruct, file=nomeArquivoLHSstruct, status="old");
       open ( unit=arqHandleLHS,       file=nomeArquivoLHS,       status="old");
       open ( unit=arqHandleRHS,       file=nomeArquivoRHS,       status="old");
       open ( unit=arqHandleIG,        file=nomeArquivoIG,        status="old");

     read (arqHandleLHSstruct ,*) ClowerLido, CupperLido, nonZerosT_;
     neq_ =  CupperLido - ClowerLido + 1
     allocate(alhs_(nonZerosT_))
     allocate(Ap_(neq_+1)      )
     allocate(Ai_(nonZerosT_))
     allocate(rows_(neq_+1)      )
     allocate(rhs_(neq_+1)      )
     allocate(x_(neq_+1)      )


    Clower_    = ClowerLido;
    Cupper_    = CupperLido;
    local_size = Cupper_ - Clower_ + 1
    Flower    = Clower_+1; Fupper    = Cupper_+1;

     print*, " criando o LHS do sistema algebrico", Flower, Fupper

     print*, " lendo      o LHS do sistema algebrico"

     Ap_(1) = 1   
     do i = Flower, Fupper
        read(arqHandleLHSstruct,* ) j, nnz ;
        do k = 1, nnz
          read(arqHandleLHSstruct,*) Ai_(Ap_(i)+k-1)
        end do
        Ap_(i+1) = Ap_(i) + nnz
    enddo
    
    read (arqHandleLHS,*)  valorDesprezadoI, valorDesprezadoF
    write(*,*) " ler o arquivo com ALHS ", nonZerosT_
    do i = 1, nonZerosT_
         read (arqHandleLHS,*)  alhs_(i)
    enddo

    print*, " criando   o  RHS do sistema algebrico"
    print*, " lendo o      RHS do sistema algebrico", Flower, Fupper
    read( arqHandleRHS,*) ClowerLido, CupperLido;
    Flower=ClowerLido+1; Fupper=CupperLido+1;
    do k = Flower, Fupper;
     read (arqHandleRHS,*)  rhs_(k), rows_(k)
    end do
    print*, " fim da leitura", Flower, Fupper

      if (.true.) then
         print*, " lendo       o initial guest do sistema algebrico"
         read (arqHandleIG,*)  ClowerLido, CupperLido;
         write(*,*)  ClowerLido, CupperLido;
         Flower=ClowerLido+1; Fupper=CupperLido+1;
         do k = Flower, Fupper;
           read (arqHandleIG,*)  x_(k), rows_(k) 
           if(k<5) write (*,*)  x_(k), rows_(k) 
           if(k>Fupper-5) write (*,*)  x_(k), rows_(k) 
         end do
      else
         print*, " atribuindo 2 o initial guest do sistema algebrico"
            x_ = 0.0
      end if

      end subroutine gerarCoefSistemaAlg

      subroutine criarSistemaAlgHYPRE(A_, parcsr_A_, b_, par_b_, x_, par_u_, solver_,  &
                                      Clower_, Cupper_,  mpi_comm_)

      integer (kind=8)              :: A_, parcsr_A_, b_, par_b_, x_, par_u_, solver_
      integer, intent(in)           :: Clower_, Cupper_
      integer (KIND=8), intent(in)  :: mpi_comm_

      integer, parameter ::  HYPRE_PARCSR=5555
      integer :: ierr 
!     Create the matrix.

!     Note that this is a square matrix, so we indicate the row partition
!     size twice (since number of rows = number of cols)

      call HYPRE_IJMatrixCreate        (mpi_comm_, Clower_, Cupper_, Clower_, Cupper_, A_, ierr )
!     Choose a parallel csr format storage (see the User's Manual)
      call HYPRE_IJMatrixSetObjectType (A_, HYPRE_PARCSR, ierr)
!     Initialize before setting coefficients
      call HYPRE_IJMatrixInitialize    (A_, ierr)

!     Create the rhs and solution
      call HYPRE_IJVectorCreate        (mpi_comm_, Clower_, Cupper_, b_, ierr )
      call HYPRE_IJVectorSetObjectType (b_, HYPRE_PARCSR, ierr)
      call HYPRE_IJVectorInitialize    (b_, ierr)

      call HYPRE_IJVectorCreate        (mpi_comm_, Clower_, Cupper_, x_, ierr )
      call HYPRE_IJVectorSetObjectType (x_, HYPRE_PARCSR, ierr)
      call HYPRE_IJVectorInitialize    (x_, ierr)

      end subroutine criarSistemaAlgHYPRE

      subroutine atribuirValoresSistemaAlgHYPRE (alhs_ , Ap_, Ai_, rhs_, x_, rows_, &
                                   neq_, nonzerosT_, Clower_, Cupper_,  &
                                A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_,  &
                                 myid_, mpi_comm_) 

      implicit none
      integer, intent(in)           :: neq_, nonzerosT_, Clower_, Cupper_
      double precision :: alhs_(nonZerosT_)  
      integer          :: Ap_(neq_), Ai_(nonZerosT_)  
      double precision :: rhs_(neq_), x_(neq_)
      integer          :: rows_(neq_)  
      integer (kind=8)              :: A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_
      integer, intent(in)           :: myid_
      integer (KIND=8), intent(in)  :: mpi_comm_

      integer       :: nnz, i, j, k, ierr
      integer       :: local_size, Flower, Fupper
      integer       :: cols(50)
      double precision ::  values(50)
     
    Flower    = Clower_+1; Fupper    = Cupper_+1;
    print*, " atribuindo valores da matrix para o HYPRE_IJMatrix "
    do i = Flower, Fupper
       j = 1
       nnz=  Ap_(i+1) - Ap_(i) 
       do k = Ap_(i), Ap_(i+1)-1
            cols(j)   = Ai_(Ap_(i) +  (k - Ap_(i)))  
            values(j) = alhs_(k) 
            j = j + 1
       end do
       call HYPRE_IJMatrixSetValues(A_, 1, nnz, i-1, cols, values, ierr)
!     Note that here we are setting one row at a time, though
!     one could set all the rows together (see the User's Manual).
    enddo

!     Assemble after setting the coefficients
      call HYPRE_IJMatrixAssemble( A_, ierr)
!     Get parcsr matrix object
      call HYPRE_IJMatrixGetObject( A_, parcsr_A_, ierr)

      local_size = Cupper_ - Clower_ + 1

    print*, " atribuindo valores de RHS para o HYPRE_IJVector "
      call HYPRE_IJVectorSetValues (b_, local_size, rows_, rhs_, ierr )
      call HYPRE_IJVectorAssemble  (b_, ierr)
      call HYPRE_IJVectorGetObject (b_, par_b_, ierr)

      print*, " atribuindo valores de RHS para o HYPRE_IJVector "
      call HYPRE_IJVectorSetValues (u_, local_size, rows_, x_, ierr)
      call HYPRE_IJVectorAssemble  (u_, ierr)
      call HYPRE_IJVectorGetObject (u_, par_u_, ierr)

      end subroutine atribuirValoresSistemaAlgHYPRE
      !call resolverSistemaAlgHYPRE        (A, parcsr_A, b, par_b, u, par_u, solver, precond, &

      subroutine resolverSistemaAlgHYPRE  (A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_, &
                                          solver_id_, num_iterations_, final_res_norm_, &
                                          myid_, mpi_comm_) 

      integer (kind=8), intent(in)   :: A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_
      integer, intent(in)            :: solver_id_
      integer, intent(out)           :: num_iterations_
      double precision , intent(out) :: final_res_norm_
      integer, intent(in)            :: myid_
      integer (KIND=8), intent(in)   :: mpi_comm_

      integer (KIND=8) ::  precond
      integer    precond_id;
      integer :: ierr

!     AMG
      if ( solver_id_ .eq. 0 ) then

!        Create solver
         write(*,*) " BoomerAMGCreate , BD"
         call HYPRE_BoomerAMGCreate(solver_, ierr)


!        Set some parameters (See Reference Manual for more parameters)

!        print solve info + parameters 
         call HYPRE_BoomerAMGSetPrintLevel  (solver_, 3,      ierr)  
!        Falgout coarsening
         call HYPRE_BoomerAMGSetCoarsenType (solver_, 6,      ierr) 
!        G-S/Jacobi hybrid relaxation 
         call HYPRE_BoomerAMGSetRelaxType   (solver_, 3,      ierr)     
!        Sweeeps on each level
         call HYPRE_BoomerAMGSetNumSweeps   (solver_, 1,      ierr)  
!         maximum number of levels 
         call HYPRE_BoomerAMGSetMaxLevels   (solver_, 20,     ierr) 
!        conv. tolerance
         call HYPRE_BoomerAMGSetTol         (solver_, 1.0d-7, ierr)    

!        Now setup and solve!
         call HYPRE_BoomerAMGSetup( solver_, parcsr_A_, par_b_, par_u_, ierr )
         call HYPRE_BoomerAMGSolve( solver_, parcsr_A_, par_b_, par_u_, ierr )


!        Run info - needed logging turned on 
         call HYPRE_BoomerAMGGetNumIterations(solver_, num_iterations_, ierr)
         call HYPRE_BoomerAMGGetFinalReltvRes(solver_, final_res_norm_, ierr)

!        Destroy solver_
         call HYPRE_BoomerAMGDestroy( solver_, ierr )

!     PCG with AMG preconditioner
      elseif (solver_id_ == 1) then

        print* , " ... PCG with AMG preconditioner, BD";
     
!        Create solver_
         call HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, solver_, ierr)

!        Set some parameters (See Reference Manual for more parameters) 
         call HYPRE_ParCSRPCGSetMaxIter    (solver_, 5000,   ierr)
         call HYPRE_ParCSRPCGSetTol        (solver_, 1.0d-7, ierr)
         call HYPRE_ParCSRPCGSetTwoNorm    (solver_, 1,      ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (solver_, 2,      ierr)
         call HYPRE_ParCSRPCGSetLogging    (solver_, 1,      ierr)

!        Now set up the AMG preconditioner and specify any parameters

         call HYPRE_BoomerAMGCreate          (precond, ierr)

!        Set some parameters (See Reference Manual for more parameters)

!        print less solver_ info since a preconditioner
         call HYPRE_BoomerAMGSetPrintLevel   (precond, 1,     ierr); 
!        Falgout coarsening
         call HYPRE_BoomerAMGSetCoarsenType  (precond, 6,     ierr) 
!        SYMMETRIC G-S/Jacobi hybrid relaxation 
         call HYPRE_BoomerAMGSetRelaxType    (precond, 6,     ierr)     
!        Sweeeps on each level
         call HYPRE_BoomerAMGSetNumSweeps    (precond, 1,     ierr)  
!        conv. tolerance
         call HYPRE_BoomerAMGSetTol          (precond, 0.0d0, ierr)     
!        do only one iteration! 
         call HYPRE_BoomerAMGSetMaxIter      (precond, 1,     ierr)

!        set amg as the pcg preconditioner
         precond_id = 2
         call HYPRE_ParCSRPCGSetPrecond      (solver_, precond_id, precond, ierr)


!        Now setup and solve!
         call HYPRE_ParCSRPCGSetup(solver_, parcsr_A_, par_b_, par_u_, ierr)
         call HYPRE_ParCSRPCGSolve(solver_, parcsr_A_, par_b_, par_u_, ierr)


!        Run info - needed logging turned on 

        call HYPRE_ParCSRPCGGetNumIterations(solver_, num_iterations_, ierr)
        call HYPRE_ParCSRPCGGetFinalRelative(solver_, final_res_norm_, ierr)

!       Destroy precond and solver

        call HYPRE_BoomerAMGDestroy(precond, ierr )
        call HYPRE_ParCSRPCGDestroy(solver_, ierr)

!     PCG with ParaSails
      elseif (solver_id_ .eq. 2) then

!        Create solver
         call HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, solver_, ierr)

!        Set some parameters (See Reference Manual for more parameters) 
         call HYPRE_ParCSRPCGSetMaxIter(solver_, 1000, ierr)
         call HYPRE_ParCSRPCGSetTol(solver_, 1.0d-7, ierr)
         call HYPRE_ParCSRPCGSetTwoNorm(solver_, 1, ierr)
         call HYPRE_ParCSRPCGSetPrintLevel(solver_, 2, ierr)
         call HYPRE_ParCSRPCGSetLogging(solver_, 1, ierr)

!        Now set up the Parasails preconditioner and specify any parameters
         call HYPRE_ParaSailsCreate(MPI_COMM_WORLD, precond,ierr)
         call HYPRE_ParaSailsSetParams(precond, 0.1d0, 1, ierr)
         call HYPRE_ParaSailsSetFilter(precond, 0.05d0, ierr)
         call HYPRE_ParaSailsSetSym(precond, 1)
         call HYPRE_ParaSailsSetLogging(precond, 3, ierr)

!        set parsails as the pcg preconditioner
         precond_id = 4
         call HYPRE_ParCSRPCGSetPrecond(solver_, precond_id, precond, ierr)


!        Now setup and solve!
         call HYPRE_ParCSRPCGSetup(solver_, parcsr_A_, par_b_, par_u_, ierr)
         call HYPRE_ParCSRPCGSolve(solver_, parcsr_A_, par_b_, par_u_, ierr)


!        Run info - needed logging turned on 

        call HYPRE_ParCSRPCGGetNumIterations(solver_, num_iterations_, ierr)
        call HYPRE_ParCSRPCGGetFinalRelative(solver_, final_res_norm_, ierr)

!       Destroy precond and solver
        call HYPRE_ParaSailsDestroy(precond, ierr )
        call HYPRE_ParCSRPCGDestroy(solver_, ierr)

      else
         if (myid_ .eq. 0) then 
           print *,'Invalid solver id specified'
           stop
         endif  
      endif

      end subroutine resolverSistemaAlgHYPRE

      subroutine destruirSistemaAlgHYPRE   (A_, b_, x_)
      integer (kind=8), intent(in) :: A_, b_, x_
      integer   ::  ierr  
         call HYPRE_IJMatrixDestroy(A_, ierr)
         call HYPRE_IJVectorDestroy(b_, ierr)
         call HYPRE_IJVectorDestroy(x_, ierr)
      end subroutine destruirSistemaAlgHYPRE       

      subroutine escreverResultados             (x_, num_iterations_, final_res_norm_, elapsedT_, myId_, &
                            print_solution_)

      double precision, allocatable :: x_(:)
      integer, intent(in) :: myId_
      integer             :: print_solution_ 
!
      integer          num_iterations_
      double precision final_res_norm_, elapsedT_

      integer   ::  ierr  
      character (LEN=16)  :: SolverStatus="incomplete"

    SolverStatus="complete"
    if(SolverStatus=="complete" .and. myid_ .eq. 0) then
        print *
        print *, "Final Relative Residual Norm = ", final_res_norm_
        print *, "Iterations                   = ", num_iterations_
        print *, 'Elapsed real time            = ', elapsedT_
        print *
    endif

!     Print the solution
      if ( print_solution_ .ne. 0 ) then
         call HYPRE_IJVectorPrint( x_, "ij.out.x", ierr)
      endif
      end subroutine escreverResultados           

      subroutine finalizarMPI                 ()
      integer :: ierr
         call MPI_Finalize(ierr)
      end subroutine finalizarMPI   

#endif

!
!=======================================================================
!    
      subroutine solverHYPRE(alhs_, brhs_, Ap_, Ai_, neq_, nonzerosT_)

      implicit none 
!
      real*8  :: alhs_(*), brhs_(neq_)
      integer :: Ap_(*), Ai_(*)
      integer, intent(in) :: neq_, nonzerosT_
! 
      real*8, allocatable  :: x(:)
!
#ifdef withHYPRE

      integer          :: myid, num_procs
      integer (KIND=8) :: mpi_comm

      integer             :: solver_id, print_solution 
      integer (KIND=8)    :: A, parcsr_A, b, par_b, u, par_u, solver
      character (LEN=16)  :: SolverStatus="incomplete"
      integer             :: Clower, Cupper
      integer             :: rows(neq_)


      integer          :: t1, t2, clock_rate, clock_max
      integer          ::  num_iterations
      double precision :: final_res_norm
      double precision :: elapsedT

      allocate(x(neq_)); x=0.0d0

      write(*,*) ' +++  iniciando a execucao do prototipo para o HYPRE '
      solver_id       = 1;
      print_solution  = 1
 
!-----------------------------------------------------------------------
!     Initialize MPI
!-----------------------------------------------------------------------
!     call inicializarMPI                 (myid, num_procs, mpi_comm)

      call criarSistemaAlgHYPRE           (A, parcsr_A, b, par_b, u, par_u, solver,  &
                                             Clower, Cupper, mpi_comm)

      call atribuirValoresSistemaAlgHYPRE (alhs_ , Ap_, Ai_, brhs_, x, rows, &
                                          neq_, nonZerosT_, Clower, Cupper,  & 
                                          A, parcsr_A, b, par_b, u, par_u, solver,  &
                                          myid, mpi_comm) 

      print* , 'solver_id =', solver_id
      call system_clock                   (t1, clock_rate, clock_max)
      call resolverSistemaAlgHYPRE        (A, parcsr_A, b, par_b, u, par_u, solver, solver_id, &
                                              num_iterations, final_res_norm, myid, mpi_comm) 
      call system_clock                   (t2, clock_rate, clock_max)
      elapsedT  =  real ( t2 - t1 ) / real ( clock_rate )

!      call escreverResultados             (x, num_iterations, final_res_norm, elapsedT,&
!                                            myId, print_solution)

!     Clean up
      call destruirSistemaAlgHYPRE        (A, b, u)

!     Finalize MPI
!      call finalizarMPI                   ()


!       call solverUMFPackPPD(Ap, Ai, alhs, x, b, neq, nonzeros)

      brhs_ = x
      deallocate(x)
#endif

      end subroutine solverHYPRE


!
!=======================================================================
!    
      subroutine criarPonteirosMatEsparsa_CRS_MR(nsd, ndof, neq, numCoefPorLinha,  &
              conectsElem, listaDosElems, id, numConexoes, numConexoesPorElem,  &
              nonzeros, simetria)

      use mGlobaisEscalares, only: optSolver

      implicit none 
!
      integer, intent(in) :: nsd, neq, numConexoes, numConexoesPorElem,nonzeros,ndof
      logical, intent(in) :: simetria
      integer, intent(out) :: numCoefPorLinha
      integer :: conectsElem(numConexoesPorElem,*), listaDosElems(numConexoesPorElem,*)
      integer :: id(ndof,*)
!
       if(numCoefPorLinha==18) then
!
         if(.not.allocated(LMstencilEqGeo))allocate(LMstencilEqGeo(neq,numCoefPorLinha)); 
         LMstencilEqGeo=0

         call montarLmStencilGeo_CRS    (LMstencilEqGeo,listaDosElems, id, &
                   conectsElem, numCoefPorLinha, ndof, numConexoes, numConexoesPorElem, neq, simetria)

         allocate(ApGeo(neq+1));    ApGeo=0  
         call montarPonteiroAp_CRS(ApGeo, LMstencilEqGeo, numCoefPorLinha, neq, nonzeros)
!
         allocate(AiGeo(nonzeros)); AiGeo=0
         call montarPonteiroAi_CRS(AiGeo, LMstencilEqGeo, numCoefPorLinha, neq, nonzeros)
!
      else
!
         if(.not.allocated(LMstencilEqVel))allocate(LMstencilEqVel(neq,numCoefPorLinha))
         LMstencilEqVel=0

         call montarLmStencilVel_CRS    (LMstencilEqVel,listaDosElems, id, &
                   conectsElem, numCoefPorLinha, ndof, numConexoes, numConexoesPorElem, neq, simetria)

         allocate(ApVel(neq+1));    ApVel=0  
         call montarPonteiroAp_CRS(ApVel, LMstencilEqVel, numCoefPorLinha, neq, nonzeros)
!
         allocate(AiVel(nonzeros)); AiVel=0
         call montarPonteiroAi_CRS(AiVel, LMstencilEqVel, numCoefPorLinha, neq, nonzeros)


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

!            write(*,'(a,i4,a,18i4)'), "LmStencil", numEq, "->",LMstencilEq(numEq,1:numCoefPorLinha)
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
     end module mSolucoesExternas

