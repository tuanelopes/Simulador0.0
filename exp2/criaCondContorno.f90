program main

implicit none

   integer :: numnp, numelM, numelF, nsd, nenB, nenF
   real*8  :: numno
   integer :: IIN, IOUT
   character*100 :: NAMEIN_coord, NAMEOUT_coord
   character*100 :: NAMEIN_conecM, NAMEOUT_conecM
   character*100 :: NAMEIN_conecF, NAMEOUT_conecF
   character*100 :: NAMEOUT_ccc_dir, NAMEOUT_ccc_sup, NAMEOUT_vcc
   character*100 :: NAMEOUT_vccfd, NAMEOUT_vccfs, NAMEOUT_ccc
   
   
   character*100 :: NAMEOUT_ccc_4bordas, NAMEOUT_ccc_esq_dir, NAMEOUT_ccc_inf_sup
   character*100 :: NAMEOUT_vcc_FlDir_LinearInfSup, NAMEOUT_vcc_FlDir_NuloInfSup 
   character*100 :: NAMEOUT_vcc_FlSup_LinearEsqDir, NAMEOUT_vcc_FlSup_NuloEsqDir
   character*100 :: NAMEOUT_pressaoMedia
   
   
   character*20 :: valor
   character*100 :: expDirName
   
   real*8 :: xmin, xmax, ymin, ymax
   real*8 :: p1, p2
   
   
   real*8, allocatable :: x(:,:)
!  
   integer :: i, tipo
   
     
   !!!!!! PARAMETROS DE ENTRADA !!!!!! 
   
   numnp=9397
   nsd=2
   nenB=3
   nenF=2
   
   xmin=0.0
   xmax=200
   ymin=0.0
   ymax=9.60009

!    CALL GetArg(1, expDirName)
!    
!    CALL GetArg(2, valor)
!    read (valor,*) p1
!    
!    CALL GetArg(3, valor)
!    read (valor,*) p2
!    
!    CALL GetArg(4, valor)
!    read (valor,*) tipo

   
   NAMEIN_coord ='nodes_rsrv_outsburden.dat'
   NAMEOUT_coord='coordenadas.inc'
!    
!    NAMEIN_conecM  =trim(expDirName)//'PlaCONEC.DAT'
!    NAMEOUT_conecM =trim(expDirName)//'conectividades.inc'
!    
!    NAMEIN_conecF  =trim(expDirName)//'LinCONEC.DAT'
!    NAMEOUT_conecF =trim(expDirName)//'conectividadesFratura.inc'
!    
!    NAMEOUT_ccc_4bordas = trim(expDirName)//'codCondContorno_4bordas.inc'
!    NAMEOUT_ccc_esq_dir = trim(expDirName)//'codCondContorno_esq_dir.inc'
!    NAMEOUT_ccc_inf_sup = trim(expDirName)//'codCondContorno_inf_sup.inc'
!    
!    NAMEOUT_pressaoMedia= trim(expDirName)//'pressaoMedia.inc'
!    
! !    NAMEOUT_vccfd = trim(expDirName)//'valCondContornoFluxoDir.inc'
! !    NAMEOUT_vccfs = trim(expDirName)//'valCondContornoFluxoSup.inc'
! !    
!    
!    NAMEOUT_vcc_FlDir_LinearInfSup = trim(expDirName)//'valCondContornoFluxoDir_LinearInfSup.inc'
!    NAMEOUT_vcc_FlDir_NuloInfSup   = trim(expDirName)//'valCondContornoFluxoDir_NuloInfSup.inc'
!    NAMEOUT_vcc_FlSup_LinearEsqDir = trim(expDirName)//'valCondContornoFluxoSup_LinearEsqDir.inc'
!    NAMEOUT_vcc_FlSup_NuloEsqDir   = trim(expDirName)//'valCondContornoFluxoSup_NuloEsqDir.inc'
!    
   !!!!!! FIM PARAMETROS DE ENTRADA !!!!!! 
   
   allocate(x(nsd,numnp))
   
   call alterarArqCoordenadas   (x, nsd,numnp,NAMEIN_coord,NAMEOUT_coord)
 
!    call alterarArqConectividadesMatriz(nenB,numelM,NAMEIN_conecM,NAMEOUT_conecM)
!    
!    call alterarArqConectividadesFratura(nenF,numelF,NAMEIN_conecF,NAMEOUT_conecF)
! !  
!    if      (tipo .EQ. 1) then
!       call gerarCodigosCondicoesContornoBordasEsq_Dir(x, nsd, numnp, xmin, xmax, ymin, ymax, NAMEOUT_ccc_esq_dir)
!       call gerarValoresCondicoesContornoFluxoDireita_FlNuloInfSup(x,nsd,numnp, xmin, xmax, ymin, ymax, &
!       &                                                         p1, p2, NAMEOUT_vcc_FlDir_NuloInfSup)
!    else if (tipo .EQ. 2) then
!       call gerarCodigosCondicoesContornoBordasInf_Sup(x, nsd, numnp, xmin, xmax, ymin, ymax, NAMEOUT_ccc_inf_sup)   
!       call gerarValoresCondicoesContornoFluxoSuperior_FlNuloEsqDir(x,nsd,numnp, xmin, xmax, ymin, ymax, &
!      &                                                           p1, p2, NAMEOUT_vcc_FlSup_NuloEsqDir)
!    else if (tipo .EQ. 3) then
!       call gerarCodigosCondicoesContorno4Bordas      (x, nsd, numnp, xmin, xmax, ymin, ymax, NAMEOUT_ccc_4bordas)
!       call gerarValoresCondicoesContornoFluxoDireita_LinearInfSup (x,nsd,numnp, xmin, xmax, ymin, ymax, &
!      &                                                          p1, p2, NAMEOUT_vcc_FlDir_LinearInfSup)   
!    else if (tipo .EQ. 4) then
!       call gerarCodigosCondicoesContorno4Bordas      (x, nsd, numnp, xmin, xmax, ymin, ymax, NAMEOUT_ccc_4bordas)
!       call gerarValoresCondicoesContornoFluxoSuperior_LinearEsqDir(x,nsd,numnp, xmin, xmax, ymin, ymax, & 
!      &                                                             p1, p2, NAMEOUT_vcc_FlSup_LinearEsqDir)
!    endif      
!    
!    call calcularPressaoMedia(p1, p2, NAMEOUT_pressaoMedia)
   
   
end program



subroutine calcularPressaoMedia(p1, p2, NAMEOUT)

   implicit none
   
   real*8, intent(in) :: p1, p2
   character(len=100) :: NAMEOUT   
!   
   integer :: IOUT
   
   IOUT=22
   OPEN(UNIT=IOUT,FILE=NAMEOUT)
   write(IOUT,'(a)') "*pressao_media_pc{"
   write(IOUT,('(es15.7)')) (p2+p1)/2.0
   write(IOUT,*) "}"
   close(IOUT)

end subroutine
!
!**** new **********************************************************************
!
subroutine gerarCodigosCondicoesContorno4Bordas(x, nsd, numnp, xmin, xmax, ymin, ymax, NAMEOUT)

   implicit none

   integer, intent(in)    :: nsd, numnp
   real*8, intent(in)     :: xmin, xmax, ymin, ymax
   real*8 , intent(inout) :: x(nsd,numnp)
   
   integer :: IIN, IOUT
   character(len=100) :: NAMEOUT
   
   integer :: i

   IOUT=22
   OPEN(UNIT=IOUT,FILE=NAMEOUT)
   
      write(IOUT,'(a)') "*codigos_cond_contorno_potencial_pc{"
   do i=1, numnp
      if(x(1,i)==xmin.or.x(1,i)>=xmax) then
         write(IOUT,('(4i10)')) i,i,1,1
      else if (x(2,i)==ymin.or.x(2,i)>=ymax) then
         write(IOUT,('(4i10)')) i,i,1,1
      endif
   enddo

   write(IOUT,'(i10)') 0
   write(IOUT,*) "}"
   close(IOUT)

end subroutine


!
!**** new **********************************************************************
!
subroutine gerarCodigosCondicoesContornoBordasInf_Sup(x, nsd, numnp, xmin, xmax, ymin, ymax, NAMEOUT)

   implicit none

   integer, intent(in)    :: nsd, numnp
   real*8, intent(in)     :: xmin, xmax, ymin, ymax
   real*8 , intent(inout) :: x(nsd,numnp)
   
   integer :: IIN, IOUT
   character(len=100) :: NAMEOUT
   
   integer :: i

   IOUT=22
   OPEN(UNIT=IOUT,FILE=NAMEOUT)
   
      write(IOUT,'(a)') "*codigos_cond_contorno_potencial_pc{"
   do i=1, numnp
      if(x(2,i)==ymin.or.x(2,i)>=ymax) then
         write(IOUT,('(4i10)')) i,i,1,1
      endif
   enddo

   write(IOUT,'(i10)') 0
   write(IOUT,*) "}"
   close(IOUT)

end subroutine

!
!**** new **********************************************************************
!
subroutine gerarCodigosCondicoesContornoBordasEsq_Dir(x, nsd, numnp, xmin, xmax, ymin, ymax, NAMEOUT)

   implicit none

   integer, intent(in)    :: nsd, numnp
   real*8, intent(in)     :: xmin, xmax, ymin, ymax
   real*8 , intent(inout) :: x(nsd,numnp)
   
   integer :: IIN, IOUT
   character(len=100) :: NAMEOUT
   
   integer :: i

   IOUT=22
   OPEN(UNIT=IOUT,FILE=NAMEOUT)
   
      write(IOUT,'(a)') "*codigos_cond_contorno_potencial_pc{"
   do i=1, numnp
      if(x(1,i)==xmin.or.x(1,i)>=xmax) then
         write(IOUT,('(4i10)')) i,i,1,1
      endif
   enddo

   write(IOUT,'(i10)') 0
   write(IOUT,*) "}"
   close(IOUT)

end subroutine

!
!**** new **********************************************************************
!
subroutine gerarValoresCondicoesContornoFluxoDireita_LinearInfSup (x,nsd,numnp, xmin, xmax, ymin, ymax, p1, p2, NAMEOUT)
! subroutine gerarValoresCondicoesContornoFluxoDireita(x,nsd,numnp, xmin, xmax, ymin, ymax, p1, p2, NAMEOUT)

   implicit none

   integer, intent(in)  :: nsd, numnp
   real*8, intent(in)   :: xmin, xmax, ymin, ymax
   real*8 , intent(inout) :: x(nsd,numnp), p1, p2
   
   integer :: IIN, IOUT
   character(len=100) :: NAMEOUT
   real*8 :: pCalc
   
   integer :: i

   IOUT=22
   OPEN(UNIT=IOUT,FILE=NAMEOUT)
   
   write(IOUT,'(a)') "*valores_cond_contorno_potencial_pc{"
   do i=1, numnp
      if(x(1,i)==xmin) then                    ! borda esquerda
         write(IOUT,('(2i10,es15.7e1)')) i,0,p1
         write(1000,('(i10,3es15.7e1)')) i,x(:,i),p1
      else if(x(1,i)>=xmax) then               ! borda direita
         write(IOUT,('(2i10,es15.7e1)')) i,0,p2
         write(1000,('(i10,3es15.7e1)')) i,x(:,i),p2
      else if((x(2,i)==ymin) .or. (x(2,i)>=ymax)) then    ! bordas inferior e superior
         pCalc=p2*(x(1,i)-xmin)/(xmax-xmin) + p1*(xmax-x(1,i))/(xmax-xmin) 
         write(IOUT,('(2i10,es15.7e1)')) i,0, pCalc
         write(1000,('(i10,3es15.7e1)')) i,x(:,i),pCalc
         write(258,('(2es10.3)')) x(1,i), pCalc
         write(259,('(i10, f15.5, f15.5)')) i, x(1,i) , pCalc
      endif
!       endif
!       endif
   enddo

   write(IOUT,'(i10)') 0
   write(IOUT,*) "}"
   close(IOUT)
   
end subroutine

!
!**** new **********************************************************************

subroutine gerarValoresCondicoesContornoFluxoSuperior_LinearEsqDir(x,nsd,numnp, xmin, xmax, ymin, ymax, p1, p2, NAMEOUT)

   implicit none

   integer, intent(in)  :: nsd, numnp
   real*8, intent(in)   :: xmin, xmax, ymin, ymax
   real*8 , intent(inout) :: x(nsd,numnp), p1, p2
   
   integer :: IIN, IOUT
   character(len=100) :: NAMEOUT
   real*8 :: pCalc
   
   integer :: i

   IOUT=22
   OPEN(UNIT=IOUT,FILE=NAMEOUT)
   
   write(IOUT,'(a)') "*valores_cond_contorno_potencial_pc{"
   do i=1, numnp
      if((x(1,i)==xmin) .or. (x(1,i)>=xmax)) then            ! bordas esquerda e direita
         pCalc=p2*(x(2,i)-ymin)/(ymax-ymin) + p1*(ymax-x(2,i))/(ymax-ymin) 
         write(IOUT,('(2i10,es15.7e1)')) i,0, pCalc
!                write(259,('(2es10.3)')) x(2,i), pCalc
      else if(x(2,i)>=ymax) then                 ! borda superior          
         write(IOUT,('(2i10,es15.7e1)')) i,0,p2
      else if(x(2,i)==ymin) then                 ! borda inferior
         write(IOUT,('(2i10,es15.7e1)')) i,0,p1
      endif
!          endif
!       endif
   enddo

   write(IOUT,'(i10)') 0
   write(IOUT,*) "}"
   close(IOUT)
   

end subroutine

!
!**** new **********************************************************************
!
subroutine gerarValoresCondicoesContornoFluxoDireita_FlNuloInfSup (x,nsd,numnp, xmin, xmax, ymin, ymax, p1, p2, NAMEOUT)
! subroutine gerarValoresCondicoesContornoFluxoDireita(x,nsd,numnp, xmin, xmax, ymin, ymax, p1, p2, NAMEOUT)

   implicit none

   integer, intent(in)  :: nsd, numnp
   real*8, intent(in)   :: xmin, xmax, ymin, ymax
   real*8 , intent(inout) :: x(nsd,numnp), p1, p2
   
   integer :: IIN, IOUT
   character(len=100) :: NAMEOUT
   real*8 :: pCalc
   
   integer :: i

   IOUT=22
   OPEN(UNIT=IOUT,FILE=NAMEOUT)
   
   write(IOUT,'(a)') "*valores_cond_contorno_potencial_pc{"
   do i=1, numnp
      if(x(1,i)==xmin) then                    ! borda esquerda
         write(IOUT,('(2i10,es15.7e1)')) i,0,p1
      else if(x(1,i)>=xmax) then               ! borda direita
         write(IOUT,('(2i10,es15.7e1)')) i,0,p2
      endif
   enddo

   write(IOUT,'(i10)') 0
   write(IOUT,*) "}"
   close(IOUT)
   
end subroutine

!
!**** new **********************************************************************

subroutine gerarValoresCondicoesContornoFluxoSuperior_FlNuloEsqDir(x,nsd,numnp, xmin, xmax, ymin, ymax, p1, p2, NAMEOUT)

   implicit none

   integer, intent(in)  :: nsd, numnp
   real*8, intent(in)   :: xmin, xmax, ymin, ymax
   real*8 , intent(inout) :: x(nsd,numnp), p1, p2
   
   integer :: IIN, IOUT
   character(len=100) :: NAMEOUT
   real*8 :: pCalc
   
   integer :: i

   IOUT=22
   OPEN(UNIT=IOUT,FILE=NAMEOUT)
   
   write(IOUT,'(a)') "*valores_cond_contorno_potencial_pc{"
   do i=1, numnp
      if(x(2,i)>=ymax) then                 ! borda superior          
         write(IOUT,('(2i10,es15.7e1)')) i,0,p2
      else if(x(2,i)==ymin) then                 ! borda inferior
         write(IOUT,('(2i10,es15.7e1)')) i,0,p1
      endif
   enddo

   write(IOUT,'(i10)') 0
   write(IOUT,*) "}"
   close(IOUT)  

end subroutine

!
!**** new **********************************************************************
!
subroutine alterarArqCoordenadas(x,nsd,numnp,NAMEIN, NAMEOUT)
   implicit none

   integer, intent(in)    :: nsd, numnp
   real*8 , intent(inout) :: x(nsd,numnp)
  
   real*8  :: Rnumel
   integer :: IIN, IOUT
   real*8  :: numno
   character(len=100) :: NAMEIN, NAMEOUT
   character*20 :: titulo
!  
   integer :: i 
   
   !COORDENADAS
 
   IIN=21
   OPEN(UNIT=IIN,FILE=NAMEIN,STATUS='OLD')
   
   IOUT=22
   OPEN(UNIT=IOUT,FILE=NAMEOUT)
   
!    read(IIN,'(a,f15.0)') titulo, Rnumel
!    read(IIN,*) titulo
   
   write(IOUT,'(a)') "*coordenadas_nodais_pc{"
   do i=1, 5207
      read (IIN,*) x(1:nsd,i)
      write(IOUT,'(i10,10x,2es15.7)') i, x(1:nsd,i)
  end do
   write(IOUT,'(i10)') 0
   write(IOUT,*) "}"

   
   close(IIN)
   close(IOUT)
 
end subroutine
!
!**** new **********************************************************************
!
subroutine alterarArqConectividadesMatriz(nen,numel,NAMEIN, NAMEOUT)

   implicit none

   integer :: nen, numel
   integer :: IIN, IOUT
   real*8  :: conects(nen), Rnumel, nel, tipo
   character*100 :: NAMEIN, NAMEOUT
   character*20 :: titulo
!  
   integer :: i 
   
   !CONECTIVIDADES
   
   IIN=21
   OPEN(UNIT=IIN,FILE=NAMEIN,STATUS='OLD')
   
   IOUT=22
   OPEN(UNIT=IOUT,FILE=NAMEOUT)
   
   read(IIN,'(a,f15.0)') titulo, Rnumel
   read(IIN,*) titulo
   
   numel=int(Rnumel)
   write(IOUT,'(a)') "*conectividades_nodais_pc{"
   do i=1, numel
      read (IIN,*) nel, tipo, conects(1:nen)
      write(IOUT,'(6i10)')  i, int(tipo), int(conects(1:nen))
   end do
   write(IOUT,'(i10)') 0
   write(IOUT,*) "}" 
   
   close(IIN)
   close(IOUT)
   
end subroutine   

!
!**** new **********************************************************************
!
subroutine alterarArqConectividadesFratura(nen,numel,NAMEIN, NAMEOUT)

   implicit none

   integer :: nen, numel
   integer :: IIN, IOUT
   real*8  :: conects(nen), Rnumel, nel, tipo
   character*100 :: NAMEIN, NAMEOUT
   character*20 :: titulo
!  
   integer :: i 
   
   !CONECTIVIDADES
   
   IIN=21
   OPEN(UNIT=IIN,FILE=NAMEIN,STATUS='OLD')
   
   IOUT=22
   OPEN(UNIT=IOUT,FILE=NAMEOUT)
   
   read(IIN,'(a,f15.0)') titulo, Rnumel
   read(IIN,*) titulo
   
   numel=int(Rnumel)
   write(IOUT,'(a)') "*conectividades_nodais_fratura_pc{"
   do i=1, numel
      read (IIN,*) nel, tipo, conects(1:nen)
      write(IOUT,'(6i10)')  i, int(tipo), int(conects(1:nen))
   end do
   write(IOUT,'(i10)') 0
   write(IOUT,*) "}" 
   
   close(IIN)
   close(IOUT)
   
end subroutine
