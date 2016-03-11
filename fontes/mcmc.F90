!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MMCMC
  INTEGER :: IFMCMC,IFPRODF,IFDISF,IFCONC,IFPRESF
  integer :: IFLAG_PRODF,IFLAG_DISF,IFLAG_CONC,IFLAG_PRESF
  REAL(8) :: TPRT_PRODF,DTPRT_PRODF,TPRT_DISF,DTPRT_DISF
  REAL(8) :: TPRT_PRESF,DTPRT_PRESF,TPRT_CONC,DTPRT_CONC
  REAL(8), DIMENSION(2,100) :: PCONDC,PCONDD,PCONDP,PCONDPR
  INTEGER, DIMENSION(100)   :: ELEM_CONDC,ELEM_CONDD,ELEM_CONDP,ELEM_CONDPR
  CHARACTER(LEN=128) :: DISF_IN,DISF_OUT,PRODF_IN,PRODF_OUT
  CHARACTER(LEN=128) :: PRESF_IN,PRESF_OUT,CONC_OUT,CONC_IN
  INTEGER :: MCMCCONC,MCMCPRODF,MCMCDISF,NCONDC,NCONDD,NCONDP
  INTEGER :: MCMCPRESF,NCONDPR,NPRESPRODWELL
  REAL(8) :: XLEFT,XRIGHT,YTOP,YBOTTOM
  REAL(8), ALLOCATABLE :: PRESPRODWELL(:)
  REAL(8) :: TPRESPRODWELL,PRESMEDIAINICIAL,PRESPROD
  REAL(8) :: ZFUNDOPOCCO
END MODULE MMCMC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE escreverArqParaview_vetor(arquivo,campo,passo,fname,nen, &
     NRESERV,nprint,LABEL,dim)
!
    USE mMalha,            only: x,nsd,numel,numnp
    USE mPropGeoFisica,    only: GEOFORM
    use mMalha,            only: conecNodaisElem, conecLadaisElem
!
    IMPLICIT NONE
!
    INTEGER :: ARQUIVO,ISTAT,d,i,n,nen,dim,j
    REAL(8), DIMENSION(dim,*) :: CAMPO
    CHARACTER(LEN=128)    :: fname,NAME
    CHARACTER(LEN=10)     :: TEMP
    REAL(8)               :: COORDZ = 0.0,PASSO
    integer               :: tipo=1
    CHARACTER(len=4)      :: EXT,LABEL
    CHARACTER(len=5)      :: C
    INTEGER               :: VARIOSARQ,ZERO,NPRINT
    LOGICAL               :: NRESERV
!
    VARIOSARQ = 1
    ZERO = 0
!
    IF(VARIOSARQ.EQ.1)THEN
       WRITE(C,300)nprint
       C=ADJUSTL(C)
       EXT='.vtk'
       NAME=ADJUSTL(TRIM(fname))//TRIM(C)//TRIM(EXT)
       OPEN(UNIT=ARQUIVO,FILE=name,STATUS='UNKNOWN',&
            ACTION='READWRITE',IOSTAT=ISTAT)
    ELSE
       WRITE(C,300)ZERO
       C=ADJUSTL(C)
       EXT='.vtk'
       NAME=ADJUSTL(TRIM(fname))//TRIM(C)//TRIM(EXT)
       OPEN(UNIT=ARQUIVO,FILE=name,STATUS='UNKNOWN',&
            ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
    END IF
    WRITE(*,*)'ARQUIVO DE IMPRESSAO:',NAME
!
    IF(ISTAT.NE.0)THEN
       WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',fname
       STOP
    END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(VARIOSARQ.EQ.1)THEN
       write(arquivo,'(a)')'# vtk DataFile Version 3.0'
       write(arquivo,'(a)')'vtk output'
       write(arquivo,'(a)')'ASCII'
       write(arquivo,'(a)')'DATASET UNSTRUCTURED_GRID'
       write(arquivo,'(a,i10,a)')'POINTS', numnp,' float '
!    write(arquivo,'(a,i10,a)')'POINTS', numnpReserv,' float '
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Escreve as coordenadas nodais
!
       if(nsd==2) then 
          do i=1,numnp
             write(arquivo,'(3(1x, 1pe15.8))') (x(d,i),d=1,nsd), coordZ 
          end do
       end if
!
       if(nsd==3) then
          do i=1,numnp
             write(arquivo,'(3(1x, 1pe15.8))') (x(d,i),d=1,nsd)
          end do
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Escreve as conectividades
       write(arquivo,'(a,i10,i10)')'CELLS', numel , (nen+1) * numel
       if(nsd==2) then
          do  n=1,numel
             write(arquivo,'(i10,9(2x,i10))') nen, (conecNodaisElem(i,n)-1, i = 1, nen) 
          end do
       end if
!
       if(nsd==3) then
          do  n=1,numel
             write(arquivo,'(i10,18(2x,i10))') nen, (conecNodaisElem(i,n)-1, i = 1, nen) 
          end do
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Escreve o tipo de celula
       write(arquivo,'(a,i10)')'CELL_TYPES ', numel
!
       if(nsd==2) then
          do  i =1,numel
             write(arquivo,'(a)') '9'!trim(adjustl(tipo))
          end do
       end if
!
       if(nsd==3) then
          do  i =1,numel
             write(arquivo,'(a)') '12'!trim(adjustl(tipo))
          end do
       end if
!
       if(tipo==1) write(arquivo,'(a,i10)')'POINT_DATA ', numnp

!    if(tipo==2) write(arquivo,'(a,i10)')'POINT_DATA',  numel*ndofP
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(TEMP,100)PASSO
    TEMP=ADJUSTL(TRIM(TEMP))
    WRITE(*,*)'TEMPO da IMPRESSAO:',TEMP
    IF(VARIOSARQ.EQ.1)THEN
       write(arquivo,'(3a,i10)')'SCALARS ', 'DISP' , ' float ',dim
    ELSE
       write(arquivo,'(3a,i10)')'SCALARS ', 't='//TRIM(TEMP)//'' , ' float ',dim
    END IF
    write(arquivo,'(2a)')'LOOKUP_TABLE ','default'
!
    DO I=1,NUMNP
       WRITE(ARQUIVO,*)(CAMPO(J,I),J=1,DIM)
    ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CLOSE(ARQUIVO) 
100 FORMAT(F10.5)
200 FORMAT(E15.7)
300 FORMAT(I5)
!
  END SUBROUTINE escreverArqParaview_vetor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE escreverArqParaview_escalar(arquivo,campo,passo,fname,nen, &
     NRESERV,nprint,LABEL,dim)
!
    USE mMalha,            only: x,nsd,numel,numnp,numelReserv,numnpReserv
    USE mPropGeoFisica,    only: GEOFORM,SRW,PHIINICIAL,PERMINICIAL
    use mMalha,            only: conecNodaisElem, conecLadaisElem
    USE mMCMC,             only: PRESMEDIAINICIAL,PRESPROD
!
    IMPLICIT NONE
!
    INTEGER :: ARQUIVO,ISTAT,d,i,n,nen,nprint,NPOSN,dim
    REAL(8), DIMENSION(dim,*) :: CAMPO
    CHARACTER(LEN=128)    :: fname,NAME
    CHARACTER(LEN=21)     :: TEMPO
    REAL(8)               :: COORDZ = 0.0,PASSO,PROP,MINIMO
    integer               :: tipo=1
    INTEGER               :: NUMNPLOCAL,NUMELLOCAL
    LOGICAL               :: NRESERV
    CHARACTER(len=4)      :: EXT,LABEL
    CHARACTER(len=5)      :: C
    INTEGER               :: VARIOSARQ,ZERO
    INTEGER, DIMENSION(numel) :: ELEMRESERV
!
    MINIMO = 1E30
    NUMNPLOCAL = numnp
    NUMELLOCAL = NUMEL
    VARIOSARQ = 1
    IF(VARIOSARQ.EQ.1)THEN
       WRITE(C,300)nprint
       C=ADJUSTL(C)
       EXT='.vtk'
       NAME=ADJUSTL(TRIM(fname))//TRIM(C)//TRIM(EXT)
    ELSE
       WRITE(C,300)ZERO
       C=ADJUSTL(C)
       EXT='.vtk'
       NAME=ADJUSTL(TRIM(fname))//TRIM(C)//TRIM(EXT)
    END IF
    WRITE(*,*)'ARQUIVO DE IMPRESSAO:',NAME
!
    IF(VARIOSARQ.EQ.1)THEN
       OPEN(UNIT=ARQUIVO,FILE=name,STATUS='UNKNOWN',&
            ACTION='WRITE',IOSTAT=ISTAT)
    ELSE
       OPEN(UNIT=ARQUIVO,FILE=name,STATUS='UNKNOWN',&
            ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
    END IF
    IF(ISTAT.NE.0)THEN
       WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',fname
       STOP
    END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    IF(VARIOSARQ.EQ.1)THEN
       write(arquivo,'(a)')'# vtk DataFile Version 3.0'
       write(arquivo,'(a)')'vtk output'
       write(arquivo,'(a)')'ASCII'
       write(arquivo,'(a)')'DATASET UNSTRUCTURED_GRID'
       write(arquivo,'(a,i10,a)')'POINTS', NUMNPLOCAL,' float '
       !    write(arquivo,'(a,i10,a)')'POINTS', numnpReserv,' float '
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Escreve as coordenadas nodais
!
       if(nsd==2) then 
          do i=1,NUMNPLOCAL
             write(arquivo,'(3(1x, 1pe15.8))') (x(d,i),d=1,nsd), coordZ 
          end do
       end if
!
       if(nsd==3) then
          do i=1,NUMNPLOCAL
             write(arquivo,'(3(1x, 1pe15.8))') (x(d,i),d=1,nsd)
          end do
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Escreve as conectividades
       write(arquivo,'(a,i10,i10)')'CELLS', NUMELLOCAL , (nen+1) * NUMELLOCAL
       if(nsd==2) then
          do  n=1,NUMELLOCAL
             write(arquivo,'(i10,9(2x,i10))') nen, (conecNodaisElem(i,n)-1, i = 1, nen) 
          end do
       end if
!
       if(nsd==3) then
          do  n=1,NUMELLOCAL
             write(arquivo,'(i10,18(2x,i10))') nen, (conecNodaisElem(i,n)-1, i = 1, nen) 
          end do
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Escreve o tipo de celula
       write(arquivo,'(a,i10)')'CELL_TYPES ', NUMELLOCAL
!
       if(nsd==2) then
          do  i =1,NUMELLOCAL
             write(arquivo,'(a)') '9'!trim(adjustl(tipo))
          end do
       end if
!
       if(nsd==3) then
          do  i =1,NUMELLOCAL
             write(arquivo,'(a)') '12'!trim(adjustl(tipo))
          end do
       end if
!
       if(tipo==1) write(arquivo,'(a,i10)')'CELL_DATA ', NUMELLOCAL
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(TEMPO,'(F20.4)')PASSO
    TEMPO=ADJUSTL(TRIM(TEMPO))
    WRITE(*,*)'TEMPO da IMPRESSAO:',TEMPO
    IF(VARIOSARQ.EQ.1)THEN
       write(arquivo,'(3a)')'SCALARS ', ' '//TRIM(LABEL)//'' , ' float '
    ELSE
       write(arquivo,'(3a)')'SCALARS ', 't='//TRIM(TEMPO)//'' , ' float '
    ENDIF
    write(arquivo,'(2a)')'LOOKUP_TABLE ','default'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(LABEL.EQ.' SAT')MINIMO = SRW
    if(LABEL.EQ.'PRES')MINIMO = PRESPROD!PRESMEDIAINICIAL
    IF(LABEL.EQ.'PORE')MINIMO = PHIINICIAL
    IF(LABEL.EQ.'PERM')MINIMO = PERMINICIAL
    IF(LABEL.EQ.'YOUN')MINIMO = 5e19
!
    DO I=1,NUMELLOCAL
       IF(NRESERV)THEN
          PROP = CAMPO(dim,I)
       ELSE
          IF(GEOFORM(I).EQ.'RESERVATORIO')THEN
             PROP=CAMPO(dim,I)
          ELSE
             PROP=MINIMO
          END IF
       END IF
       WRITE(ARQUIVO,*)PROP
    ENDDO
!
    CLOSE(ARQUIVO) 
100 FORMAT(F10.5)
200 FORMAT(E15.7)
300 FORMAT(I5)
!
END SUBROUTINE escreverArqParaview_escalar
!=======================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE MONITORMASSATOTAL(U,POROS,T,iifmass_out)
      use mMalha,            only: numelReserv,numel
      use mPropGeoFisica,    only: hy,hx
      use mPropGeoFisica,    only: rhow
      use mPropGeoFisica,    only: GEOFORM
!
      IMPLICIT NONE
!
      INTEGER :: N,INFILE,NEL
      REAL(8) :: T,VOL
      REAL(8), DIMENSION(numelReserv) :: U,POROS
      INTEGER :: ISTAT,I,ZERO=0D0
      CHARACTER(len=128) :: NAME,FILE_IN
      CHARACTER(len=4)   :: EXT
      CHARACTER(len=5)   :: C
      CHARACTER(len=128) :: iifmass_out
!
      INFILE = 170
      WRITE(C,113)ZERO
      C=ADJUSTL(C)
      FILE_IN = IIFMASS_OUT
      EXT='.dat'
      NAME=TRIM(FILE_IN)//TRIM('_')//TRIM(C)//TRIM(EXT)
      NAME=ADJUSTL(TRIM(NAME))
!
      OPEN(UNIT=INFILE,FILE=NAME,STATUS='UNKNOWN',&
           ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
      IF(ISTAT.NE.0)THEN
         WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',IIFMASS_OUT
         STOP
      END IF
!
      VOL=0.0D0
      DO NEL=1,numel
         IF(GEOFORM(NEL).EQ.'RESERVATORIO')THEN
            VOL = VOL+U(NEL)*POROS(NEL)
         END IF
      END DO
      write(*,*)'##### MONITORANDO MASSA TOTAL DE AGUA #####'
      WRITE(INFILE,400)T,VOL*RHOW*HX*HY
!
      CLOSE(INFILE)
!
113   FORMAT(I5)
400   FORMAT(100(e15.8,2x))
    END SUBROUTINE MONITORMASSATOTAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ENCONTRARLIMITESRESERVATORIO(POINTINF,POINTSUP)
!
    USE mMalha, only: x,nsd,numel,numnp, conecNodaisElem,NEN
    USE mMalha, only: numelReserv
    USE mPropGeoFisica, only: GEOFORM
    IMPLICIT NONE
!
    INTEGER :: J,I,K,NO,NELEMI,NELEMF
    REAL(8),DIMENSION(NSD) :: POINTINF,POINTSUP
    REAL(8) :: AUX,AUXMAX,AUXMIN
    REAL(8) :: DIMRESERVX,DIMRESERVY

!
    AUXMAX=-1E20
    AUXMIN=1E20
!
!    DO J=1,NUMEL
!       IF(GEOFORM(J).EQ.'RESERVATORIO')THEN
!          DO I=1,NEN
!             NO=conecNodaisElem(I,J)-1
!             IF(NSD.EQ.2)THEN
!                AUX=SQRT(X(1,NO)*X(1,NO)+X(2,NO)*X(2,NO))
!                IF(AUX.GT.AUXMAX)THEN
!                   AUXMAX=AUX
!                   NELEMF=J
!                   POINTSUP(1)=X(1,NO)
!                   POINTSUP(2)=X(2,NO)
!                END IF
!                IF(AUX.LT.AUXMIN)THEN
!                   AUXMIN=AUX
!                   NELEMI=J
!                   POINTINF(1)=X(1,NO)
!                   POINTINF(2)=X(2,NO)
!                END IF
!             END IF
!          END DO
!       END IF
!    END DO
!ACHAR AS COORDENADAS INICIAIS E FINAIS DO RESERVATORIO
    NELEMI=1
    NELEMF=numelReserv
!
    NO=conecNodaisElem(1,NELEMI)
    POINTINF(1)=X(1,NO)
    POINTINF(2)=X(2,NO)
    NO=conecNodaisElem(3,NELEMF)
    POINTSUP(1)=X(1,NO)
    POINTSUP(2)=X(2,NO)
!
    DIMRESERVX=POINTSUP(1)-POINTINF(1)
    DIMRESERVY=POINTSUP(2)-POINTINF(2)
      WRITE(*,*)'COORDENADA INICIAL DO RESERVATORIO:',POINTINF
      WRITE(*,*)'COORDENADA FINAL DO RESERVATORIO..:',POINTSUP
      WRITE(*,*)'DIMENSOES DO RESERVATORIO.........:',DIMRESERVX,DIMRESERVY
!
    END SUBROUTINE ENCONTRARLIMITESRESERVATORIO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE INITMCMC
    USE MMCMC
    USE mGlobaisEscalares, only: tt,tzero,NVEL

!
    IMPLICIT NONE
    CHARACTER(LEN=128) :: FNOME,FLAG,FLAG1
    INTEGER :: ISTAT
    REAL(8) :: TEMPO, TOL=1E-6
!    
    character(len=50) :: keyword_name
    integer :: ierr
!

    keyword_name='mcmc_concentracao'
    call readRandMCMC(keyword_name,IFLAG_CONC, CONC_IN, CONC_OUT, MCMCCONC, ierr)
    WRITE(*,*)'#####################'
    WRITE(*,*)'ARQUIVO MCMC CONC......................:',ADJUSTL(TRIM(CONC_OUT))
    WRITE(*,*)'NUMERO DE IMPRESSOES:',MCMCCONC
    WRITE(*,*)'ARQUIVO COM POSICOES DE MONITORAMENTO..:',ADJUSTL(TRIM(CONC_IN))   
    IF(IFLAG_CONC==1) CALL COND_POINTS(CONC_IN,NCONDC,PCONDC)
    WRITE(*,*)'#####################'
!
    keyword_name='mcmc_producao'
    call readRandMCMC(keyword_name,IFLAG_PRODF, PRODF_IN, PRODF_OUT, MCMCPRODF, ierr)
    WRITE(*,*)'#####################'
    WRITE(*,*)'ARQUIVO MCMC PROD......................:',ADJUSTL(TRIM(PRODF_OUT))
    WRITE(*,*)'NUMERO DE IMPRESSOES:',MCMCPRODF
    WRITE(*,*)'ARQUIVO COM POSICOES DE MONITORAMENTO..:',ADJUSTL(TRIM(PRODF_IN))
    IF(IFLAG_PRODF==1) CALL COND_POINTS(PRODF_IN,NCONDP,PCONDP)
    WRITE(*,*)'#####################'
!
    keyword_name='mcmc_deslocamentos'
    call readRandMCMC(keyword_name,IFLAG_DISF, DISF_IN, DISF_OUT, MCMCDISF, ierr)
    WRITE(*,*)'#####################'
    WRITE(*,*)'ARQUIVO MCMC DISLOC....................:',ADJUSTL(TRIM(DISF_OUT))
    WRITE(*,*)'NUMERO DE IMPRESSOES:',MCMCDISF
    WRITE(*,*)'ARQUIVO COM POSICOES DE MONITORAMENTO..:',ADJUSTL(TRIM(DISF_IN))
    IF(IFLAG_DISF==1) CALL COND_POINTS(DISF_IN,NCONDD,PCONDD)
    WRITE(*,*)'#####################'
!
    keyword_name='mcmc_pressao'
    call readRandMCMC(keyword_name,IFLAG_PRESF, PRESF_IN, PRESF_OUT, MCMCPRESF, ierr)
    WRITE(*,*)'#####################'
    WRITE(*,*)'ARQUIVO MCMC PROD.....................:',ADJUSTL(TRIM(PRESF_OUT))
    WRITE(*,*)'NUMERO DE IMPRESSOES:',MCMCPRESF
    WRITE(*,*)'ARQUIVO COM POSICOES DE MONITORAMENTO..:',ADJUSTL(TRIM(PRESF_IN))
    IF(IFLAG_PRESF==1) CALL COND_POINTS(PRESF_IN,NCONDPR,PCONDPR)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TEMPO=(tt-tzero)/REAL(NVEL)
    DTPRT_CONC = (tt-tzero)/REAL(MCMCCONC)
    IF(MCMCCONC.GE.NVEL)THEN
       MCMCCONC = NVEL
       DTPRT_CONC = (tt-tzero)/REAL(MCMCCONC)
    ELSE
       DO WHILE(DMOD(DTPRT_CONC,TEMPO).GT.TOL)
          MCMCCONC = MCMCCONC+1
          DTPRT_CONC = (tt-tzero)/REAL(MCMCCONC)
       END DO
    END IF
    TPRT_CONC = DTPRT_CONC
    write(*,*)'NUMERO de IMPRESSOES CONC',MCMCCONC
    write(*,*)'Delta t IMPRESSAO        ',DTPRT_CONC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DTPRT_PRODF= (tt-tzero)/REAL(MCMCPRODF)
    IF(MCMCPRODF.GE.NVEL)THEN
       MCMCPRODF= NVEL
       DTPRT_PRODF= (tt-tzero)/REAL(MCMCPRODF)
    ELSE
       DO WHILE(DMOD(DTPRT_PRODF,TEMPO).GT.TOL)
          MCMCPRODF= MCMCPRODF+1
          DTPRT_PRODF= (tt-tzero)/REAL(MCMCPRODF)
       END DO
    END IF
    TPRT_PRODF=DTPRT_PRODF  
    write(*,*)'NUMERO de IMPRESSOES PROD',MCMCPRODF
    write(*,*)'Delta t IMPRESSAO        ',DTPRT_PRODF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DTPRT_DISF= (tt-tzero)/REAL(MCMCDISF)
    IF(MCMCDISF.GE.NVEL)THEN
       MCMCDISF= NVEL
       DTPRT_DISF= (tt-tzero)/REAL(MCMCDISF)
    ELSE
       DO WHILE(DMOD(DTPRT_DISF,TEMPO).GT.TOL)
          MCMCDISF= MCMCDISF+1
          DTPRT_DISF= (tt-tzero)/REAL(MCMCDISF)
       END DO
    END IF
    TPRT_DISF = DTPRT_DISF
    write(*,*)'NUMERO de IMPRESSOES DISP',MCMCDISF
    write(*,*)'Delta t IMPRESSAO        ',DTPRT_DISF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DTPRT_PRESF= (tt-tzero)/REAL(MCMCPRESF)
    IF(MCMCPRESF.GE.NVEL)THEN
       MCMCPRESF= NVEL
       DTPRT_PRESF= (tt-tzero)/REAL(MCMCPRESF)
    ELSE
       DO WHILE(DMOD(DTPRT_PRESF,TEMPO).GT.TOL)
          MCMCPRESF= MCMCPRESF+1
          DTPRT_PRESF= (tt-tzero)/REAL(MCMCPRESF)
       END DO
    END IF
    TPRT_PRESF=DTPRT_PRESF  
    write(*,*)'NUMERO de IMPRESSOES PRES',MCMCPRESF
    write(*,*)'Delta t IMPRESSAO        ',DTPRT_PRESF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(PRESPRODWELL(NCONDP))
    PRESPRODWELL = 0.0D0
!
    CALL MONITORPOSITION(PCONDC, NCONDC, ELEM_CONDC)
    CALL MONITORPOSITION(PCONDP, NCONDP, ELEM_CONDP)
    CALL MONITORPOSITION(PCONDD, NCONDD, ELEM_CONDD)
    CALL MONITORPOSITION(PCONDPR,NCONDPR,ELEM_CONDPR)
!
    contains
                             
      subroutine readRandMCMC(keyword_name,IFLAG, ArqIN, ArqOUT, MCMC, ierr)
!
        use mInputReader,    only: file_lines, findKeyword      
        
        implicit none
        
        integer keyword_line, ierr
        character(len=50) keyword_name
        integer :: IFLAG, MCMC
        character(len=128) :: ArqIN, ArqOUT

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
            ierr=1
            return
        end if
        read(file_lines(keyword_line:),'(a)') IFLAG
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),'(a)') ArqOut
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),*) MCMC
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),'(a)') ArqIn

    end subroutine 
    
  END SUBROUTINE INITMCMC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MONITORPOSITION(PCOND,N,ELEM_COND)
    use mMalha,         only: x,NUMEL,conecNodaisElem
    IMPLICIT NONE
!
    INTEGER :: N,J,I,NO,NEL
    REAL(8), DIMENSION(2,100) :: PCOND
    REAL(8) :: X1,X2,Y1,Y2
    INTEGER, DIMENSION(100)   :: ELEM_COND
!
    IF(N.GT.100)THEN
       WRITE(*,*)'NUMERO EXCESSIVO DE PONTOS'
       WRITE(*,*)'MODIFICAR EM MODULE MMCMC '
       STOP
    END IF
!
    WRITE(*,*)'########################################'
    WRITE(*,*)'### MONITORAMENTO ######################'
    WRITE(*,*)'### NUMERO DE PONTOS:',N
    WRITE(*,*)'########################################'
    DO 200 I=1,N
       DO NEL=1,NUMEL
          X1=x(1,conecNodaisElem(1,NEL))
          X2=x(1,conecNodaisElem(2,NEL))
          Y1=x(2,conecNodaisElem(1,NEL))
          Y2=x(2,conecNodaisElem(4,NEL))
          IF(PCOND(1,I).LT.X2)THEN
             IF(PCOND(1,I).GE.X1)THEN
                IF(PCOND(2,I).LT.Y2)THEN
                   IF(PCOND(2,I).GE.Y1)THEN
                      ELEM_COND(I)=NEL
                      GO TO 200
                   END IF
                END IF
             END IF
          END IF
       END DO
200    CONTINUE
!
  END SUBROUTINE MONITORPOSITION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MONITOR(FNOME,PCOND,N,ELEM_COND,U,T)
    IMPLICIT NONE
!
    INTEGER :: N,INFILE
    REAL(8) :: T
    REAL(8), DIMENSION(*) :: U
    REAL(8), DIMENSION(2,100) :: PCOND
    CHARACTER(LEN=128) :: FNOME
    INTEGER :: ISTAT,I,ZERO=0D0
    INTEGER, DIMENSION(100)   :: ELEM_COND
    CHARACTER(len=128) :: NAME,FILE_IN
    CHARACTER(len=4)   :: EXT
    CHARACTER(len=5)   :: C
!
    INFILE = 163
    WRITE(C,113)ZERO
    C=ADJUSTL(C)
    FILE_IN = FNOME
    EXT='.dat'
    NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
    NAME=ADJUSTL(TRIM(NAME))
!
    OPEN(UNIT=INFILE,FILE=NAME,STATUS='UNKNOWN',&
         ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
    IF(ISTAT.NE.0)THEN
       WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FNOME
       STOP
    END IF
!
    write(*,*)'##### MONITORANDO #####'
    WRITE(INFILE,400)T,(U(ELEM_COND(I)),I=1,N)
!
    CLOSE(INFILE)
!
113 FORMAT(I5)
400 FORMAT(100(e15.8,2x))
  END SUBROUTINE MONITOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MONITORV(FNOME,PCOND,N,ELEM_COND,U,T)
    use mGlobaisEscalares, only : ndofD
    USE mMalha,             only: x,numnp,conecNodaisElem
    IMPLICIT NONE
!
    INTEGER :: N,INFILE,NEL,J
    REAL(8) :: T
    REAL(8), DIMENSION(ndofD,numnp) :: U
    REAL(8), DIMENSION(2,100) :: PCOND
    CHARACTER(LEN=128) :: FNOME
    INTEGER :: ISTAT,I,ZERO=0D0
    INTEGER, DIMENSION(100)   :: ELEM_COND
    CHARACTER(len=128) :: NAME,FILE_IN
    CHARACTER(len=4)   :: EXT
    CHARACTER(len=5)   :: C
    REAL(8) :: X1,X2,Y1,Y2
!
    INFILE = 164
    WRITE(C,113)ZERO
    C=ADJUSTL(C)
    FILE_IN = FNOME
    EXT='.dat'
    NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
    NAME=ADJUSTL(TRIM(NAME))
!
    OPEN(UNIT=INFILE,FILE=NAME,STATUS='UNKNOWN',&
         ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
    IF(ISTAT.NE.0)THEN
       WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FNOME
       STOP
    END IF
!
    write(*,*)'##### MONITORANDO #####'
    WRITE(INFILE,400)T,(U(2,conecNodaisElem(3,ELEM_COND(NEL))),NEL=1,N)
!
    CLOSE(INFILE)
!
113 FORMAT(I5)
400 FORMAT(100(e15.8,2x))
  END SUBROUTINE MONITORV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IREADSTATM(flag,ifile,N)
!
      implicit none
      integer :: ifile,NTESTE
      LOGICAL :: N
      character(len=128) :: flag,flag1
!
      read(ifile,"(a)") flag1
      if(trim(flag).eq.trim(flag1)) then
      read(ifile,*) NTESTE
      IF(NTESTE.EQ.1)THEN
         N=.TRUE.
      ELSE
         N=.FALSE.
      END IF
      else
      write(*,*) "Erro na leitura de ", flag
      stop
      end if
!
    END SUBROUTINE IREADSTATM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE COND_POINTS(file_in,N,PCOND)
!
      implicit none
!
      integer, intent(out)           :: N
      integer                        :: idata,i,j,ISTAT
      character(len=128), intent(in) :: FILE_IN
      real(8)                        :: TOL
      REAL(8), intent(out), DIMENSION(2,100) :: PCOND
!
      TOL=1.e-6
!
      idata=168
      OPEN(UNIT=IDATA,FILE=FILE_IN,STATUS='UNKNOWN',&
           ACTION='READ',IOSTAT=ISTAT)
      IF(ISTAT.NE.0)THEN
         WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_IN
         STOP
      ELSE
         WRITE(*,*)'READING FILE:',ADJUSTL(TRIM(FILE_IN))
      END IF
!
      read(idata,*)N
      write(*,101)N
!
      do i=1,N
      read(idata,*)(pcond(j,i),j=1,2)
      write(*,201)i,pcond(1,i),pcond(2,i)
      enddo
!
      do i=1,N
         do j=1,2
            pcond(j,i)=pcond(j,i)+TOL
         enddo
      enddo
      CLOSE(IDATA)
!
100   format(I5)
101   format('Numero de Pontos de Controle:',i5,/)
200   format(3F12.7)
201   format('Coordenada do Ponto',i5,' (',F10.4,',',F10.4,')')
!
    end SUBROUTINE COND_POINTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE MAPEIARESERVATORIO(permk,perm,conecNodaisElem,x,xlx,xly,xlz, &
                                  nelemx,nelemy,nelemz,nelem,nen,nsd)
!
!     MAPEA O CAMPO DE PERMEABILIDADES NO DOMINIO DE SIMULACAO
!
      use mMalha,            only: numel
      use mGlobaisEscalares, only: geomech
      use mPropGeoFisica,    only: nelxReserv, nelyReserv, nelzReserv, &
                                   dimx, dimy, dimz, hy, hx, hz
!
      implicit none
!
      integer :: nelemx,nelemy,nelemz,nelem,nen,nsd,NOI,NOF
      integer,dimension(nen,*)  :: conecNodaisElem
      real(8),dimension(nsd,*)  :: x
      real(8),dimension(*)      :: perm,permk
!      
      integer :: nel,i,no,j,k,ncont,NBLOCKX,NBLOCKY
      real(8) :: xi,xf,yi,yf,zi,zf,xx,yy,zz,xlx,xly,xlz,TOL
      real(8) :: cdx,cdy,cdz,aux,aux1,aux2,X1,X2,Y1,Y2
      REAL(8),DIMENSION(NSD) :: POINTINF,POINTSUP,CENTRO
      REAL(8) :: XPOS,YPOS
!
      TOL=1E-6
      XPOS=0.0D0
      YPOS=0.0D0
!
      if(nsd==2) then
         nelzReserv=1
         nelemz=1
      endif
!
!     definicacao do tamanlhos dos elementos
!     Campo de permeabilidades
!
      cdx=xlx/nelemx
      cdy=xly/nelemy
      if(nsd==3)cdz=xlz/nelemz
!
!     Impressao dos parametros
!
      IF(nsd==2)THEN
         write(*,800)dimx,dimy,nelxReserv,nelyReserv
         write(*,890)xlx,xly,nelemx,nelemy
         write(*,891)hx,hy,cdx,cdy
      ELSE
         write(*,900)dimx,dimy,dimz,nelxReserv,nelyReserv,nelzReserv
         write(*,990)xlx,xly,xlz,nelemx,nelemy,nelemz
         write(*,991)hx,hy,hz,cdx,cdy,cdz
      END IF
!
!     verificacao da consistencia
!
      aux = (hx-cdx)
      aux1= (hy-cdy)
      if(nsd==3) aux2= (hz-cdz)
!
      if(nsd==2) then
         if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,110)
            stop
         end if
      else
         if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,110)
            stop
         end if
      end if
!
      aux = (dimx-xlx)
      aux1= (dimy-xly)
      if(nsd==3) aux2= (dimz-xlz)
      tol = 1.0e-08
!
      if(nsd==2) then
         if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,111)
            stop
         end if
      else
         if(aux.gt.tol.or.aux1.gt.tol.or.aux2.gt.tol)then
            write(*,111)
            stop
         end if
      endif
      CALL ENCONTRARLIMITESRESERVATORIO(POINTINF,POINTSUP)
!
!     loop nos elementos
!             
      nelem=nelxReserv*nelyReserv
!
      DO NEL=1,NUMEL
         X1=X(1,conecNodaisElem(1,NEL))+TOL
         X2=X(1,conecNodaisElem(3,NEL))-TOL
         Y1=X(2,conecNodaisElem(1,NEL))+TOL
         Y2=X(2,conecNodaisElem(3,NEL))-TOL
         IF(X1.GE.POINTINF(1))THEN
            IF(X2.LE.POINTSUP(1))THEN
              IF(Y1.GE.POINTINF(2))THEN
                 IF(Y2.LE.POINTSUP(2))THEN
                    CENTRO(1)=(X2+X1)*0.5D0
                    CENTRO(2)=(Y1+Y2)*0.5D0
                    I=0
                    DO NBLOCKY=1,nelemy
                       YPOS = POINTINF(2)+REAL(NBLOCKY-1)*CDY
                       DO NBLOCKX=1,nelemx
                          I=I+1
                          XPOS = POINTINF(1)+REAL(NBLOCKX-1)*CDX
                          IF((CENTRO(1).GE.XPOS).AND.(CENTRO(1).LE.(XPOS+CDX)))THEN
                             IF((CENTRO(2).GE.YPOS).AND.(CENTRO(2).LE.(YPOS+CDY)))THEN
                                PERMK(NEL) = PERM(I)
                                exit
                             END IF
                          END IF
                       END DO
                    END DO
                 END IF
              END IF
            END IF
         END IF
      END DO
!      STOP
!
!
 800  format(/,                            &
       '##############################',/, &
       'TAMANHO DO RESERVATORIO:',/,'Lx=', &
       f10.5,2x,'Ly=',f10.5,/,             &
       'MALHA:',/,'nx=',i7,' ny=',i7,/,    &
       '##############################')
!
 890  format(/,                              &
       '##############################',/,   &
       'TAMANHO DO CAMPO ALEATORIO:',/,'Lx=',&
       f10.5,2x,'Ly=',f10.5,/,               &
       'MALHA:',/,'nx=',i7,' ny=',i7,/,      &
       '##############################')
!
 891  format(/,                            &
       '##############################',/, &
       'Tamanho do elemento:',/,'hx=',     &
       f10.5,2x,'hy=',f10.5,/, &
       '##############################',/, &
       'Tamanho do bloco geologico:',/,'gx=', &
       f10.5,2x,'gy=',f10.5,/,&
       '##############################',/)

 900  format(/,                            &
       '##############################',/, &
       'TAMANHO DO RESERVATORIO:',/,'Lx=', &
       f10.5,2x,'Ly=',f10.5,2x,'Lz=',f10.5,/,&
       'MALHA:',/,'nx=',i7,' ny=',i7,' nz=',i7,/,    &
       '##############################')
!
 990  format(/,                              &
       '##############################',/,   &
       'TAMANHO DO CAMPO ALEATORIO:',/,'Lx=',&
       f10.5,2x,'Ly=',f10.5,2x,'Lz=',f10.5,/,&
       'MALHA:',/,'nx=',i7,' ny=',i7,' nz=',i7,/,&
       '##############################')
!
 991  format(/,                            &
       '##############################',/, &
       'Tamanho do elemento:',/,'hx=',     &
       f10.5,2x,'hy=',f10.5,2x,'hz=',f10.5,/, &
       '##############################',/, &
       'Tamanho do bloco geologico:',/,'gx=', &
       f10.5,2x,'gy=',f10.5,2x,'gz=',f10.5,/,&
       '##############################',/)

 110  format('####################################',/, &
             'INCONSISTENCIA NO TAMANHO DOS BLOCOS',/, &
             'GEOLOGICOS: MALHA COMPUTACIONAL     ',/, &
             'MAIS GROSSEIRA QUE A MALHA GEOLOGICA',/, &
             '####################################',/)

 111  format('####################################',/, &
             'INCONSISTENCIA NO TAMANHO DOS       ',/, &
             'DOMINIOS: O DOMINIO GEOLOGICO DEVE  ',/, &
             'SER MAIOR OU IGUAL AO DOMINIO       ',/, &
             '####################################',/)
!
      return
!
    END SUBROUTINE MAPEIARESERVATORIO
!
!=======================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE MONITORPROD(FNOME,PCOND,N,ELEM_COND,VEL,U,T)
      use mMalha,            only: numLadosReserv,numel
      use mMalha,            only: conecLadaisElem
      use mGlobaisEscalares, only: ndofV
      use mPropGeoFisica,    only: hy,xlw,xlo,XLT
      use MMCMC,             only: dtprt_prodF
!
      IMPLICIT NONE
!
      INTEGER :: N,INFILE,NEL,LADO
      REAL(8) :: T,PROD,OIL_MOBI,ST
      REAL(8), DIMENSION(numel) :: U
      REAL(8), DIMENSION(ndofV,numLadosReserv) :: VEL
      REAL(8), DIMENSION(2,100) :: PCOND
      CHARACTER(LEN=128) :: FNOME
      INTEGER :: ISTAT,I,ZERO=0D0
      INTEGER, DIMENSION(100)   :: ELEM_COND
      CHARACTER(len=128) :: NAME,FILE_IN
      CHARACTER(len=4)   :: EXT
      CHARACTER(len=5)   :: C
!
      INFILE = 169
      WRITE(C,113)ZERO
      C=ADJUSTL(C)
      FILE_IN = FNOME
      EXT='.dat'
      NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
      NAME=ADJUSTL(TRIM(NAME))
!
      OPEN(UNIT=INFILE,FILE=NAME,STATUS='UNKNOWN',&
           ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
      IF(ISTAT.NE.0)THEN
         WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FNOME
         STOP
      END IF
!
      PROD=0.0D0
      DO NEL=1,N
         ST   = U(ELEM_COND(NEL))
         LADO = conecLadaisElem(2,ELEM_COND(NEL))
         OIL_MOBI = XLO(ST)/XLT(ST)
         PROD = PROD + VEL(1,LADO)*OIL_MOBI
      END DO
      write(*,*)'##### MONITORANDO PRODUCAO #####'
      WRITE(INFILE,400)T,PROD*HY!*dtprt_prodF
!
      CLOSE(INFILE)
!
113   FORMAT(I5)
400   FORMAT(100(e15.8,2x))
    END SUBROUTINE MONITORPROD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE MAPEIADOMINIO(permk,perm,conecNodaisElem,x,xlx,xly,xlz, &
                                  nelemx,nelemy,nelemz,nelem,nen,nsd,YX1,YX2)
!
!     MAPEA O CAMPO DE PERMEABILIDADES NO DOMINIO DE SIMULACAO
!
      use mMalha,            only: numel,numnp
      use mGlobaisEscalares, only: geomech
      use mPropGeoFisica,    only: nelxReserv, nelyReserv, nelzReserv, &
                                   dimx, dimy, dimz, hy, hx, hz
      use MMCMC,             only: XLEFT,XRIGHT,YTOP,YBOTTOM
!
      implicit none
!
      integer :: nelemx,nelemy,nelemz,nelem,nen,nsd,NOI,NOF
      integer,dimension(nen,*)  :: conecNodaisElem
      real(8),dimension(nsd,*)  :: x
      real(8),dimension(*)      :: perm,permk
!      
      integer :: nel,i,no,j,k,ncont,NBLOCK
      real(8) :: xi,xf,yi,yf,zi,zf,xx,yy,zz,xlx,xly,xlz,TOL
      real(8) :: cdx,cdy,cdz,aux,aux1,aux2,X1,X2,Y1,Y2
      REAL(8),DIMENSION(NSD) :: POINTINF,POINTSUP,CENTRO
      REAL(8) :: XPOS,YPOS,AREA,YX1,YX2
!
! PONTO INICIAL SERA A POSICAO (XLEFT,YTOP)
!      WRITE(*,*)'()()()MAPEANDO()()()'
!      WRITE(*,*) 'XLEFT = ', XLEFT,'XRIGHT  = ',XRIGHT
!      WRITE(*,*) 'YTOP  = ',  YTOP,'YBOTTOM = ',YBOTTOM 
!
      TOL=1E-6
      XPOS=0.0D0
      YPOS=0.0D0
!
      if(nsd==2) then
         nelzReserv=1
         nelemz=1
      endif
!
!     definicacao do tamanlhos dos elementos
!     Campo de permeabilidades
!
      cdx=xlx/nelemx
      cdy=xly/nelemy
      if(nsd==3)cdz=xlz/nelemz
!
!     Impressao dos parametros
!
      IF(nsd==2)THEN
!         write(*,800)dimx,dimy,nelxReserv,nelyReserv
!         write(*,890)xlx,xly,nelemx,nelemy
         write(*,891)hx,hy,cdx,cdy
      ELSE
         write(*,900)dimx,dimy,dimz,nelxReserv,nelyReserv,nelzReserv
         write(*,990)xlx,xly,xlz,nelemx,nelemy,nelemz
         write(*,991)hx,hy,hz,cdx,cdy,cdz
      END IF
!
!     verificacao da consistencia
!
      AUX = -10.0
      AUX1= -10.0
      DO NEL=1,NUMEL
         X1=X(1,conecNodaisElem(1,NEL))
         X2=X(1,conecNodaisElem(3,NEL))
         Y1=X(2,conecNodaisElem(1,NEL))
         Y2=X(2,conecNodaisElem(3,NEL))
         IF(AUX.LT.ABS(X2-X1))AUX=ABS(X2-X1)
         IF(AUX1.LT.ABS(Y2-Y1))AUX1=ABS(Y2-Y1)
      END DO
      WRITE(*,*)'###################################'
      WRITE(*,*)'###################################'
      WRITE(*,*)'## MAPEANDO MODULO DE YOUNG #######'
      WRITE(*,*)'###################################'
      WRITE(*,*)'### PONTO INICIAL DE MAPEAMENTO ###'
      WRITE(*,*)'###################################'
      WRITE(*,*)'POSICAO X:',YX1
      WRITE(*,*)'POSICAO Y:',YX2
      WRITE(*,*)'###################################'
      WRITE(*,*)'DIMENSAO MAXIMA DE UM ELEMENTO EM X'
      WRITE(*,*)AUX
      WRITE(*,*)'DIMENSAO MAXIMA DE UM ELEMENTO EM Y'
      WRITE(*,*)AUX1
      WRITE(*,*)'###################################'
      WRITE(*,*)'###################################'
!
      aux = (AUX-cdx)
      aux1= (AUX1-cdy)
      if(nsd==3) aux2= (hz-cdz)
!
      if(nsd==2) then
         if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,110)
            stop
         end if
      else
         if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,110)
            stop
         end if
      end if
!
      aux = (dimx-xlx)
      aux1= (dimy-xly)
      if(nsd==3) aux2= (dimz-xlz)
      tol = 1.0e-08
!
      if(nsd==2) then
         if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,111)
            stop
         end if
      else
         if(aux.gt.tol.or.aux1.gt.tol.or.aux2.gt.tol)then
            write(*,111)
            stop
         end if
      endif
!      CALL ENCONTRARLIMITESRESERVATORIO(POINTINF,POINTSUP)
      POINTINF(1)=XLEFT
      POINTINF(2)=YBOTTOM
      POINTSUP(1)=XRIGHT
      POINTSUP(2)=YTOP
!
!     loop nos elementos
!             
      nelem=nelemx*nelemy
!
      DO NEL=1,NUMEL
         AREA=(X(1,conecNodaisElem(4,NEL))*X(2,conecNodaisElem(1,NEL))- &
                 X(1,conecNodaisElem(1,NEL))*X(2,conecNodaisElem(4,NEL)))
         DO J=1,NEN-1
            AREA = AREA+(X(1,conecNodaisElem(J,NEL))*X(2,conecNodaisElem(J+1,NEL))- &
                 X(1,conecNodaisElem(J+1,NEL))*X(2,conecNodaisElem(J,NEL)))
         END DO
         AREA=AREA*0.5
         AREA=1.0D0/(6.0*AREA)
         AUX=(X(1,conecNodaisElem(nen,NEL))*X(2,conecNodaisElem(1,NEL))-  &
                 X(1,conecNodaisElem(1,NEL))*X(2,conecNodaisElem(nen,NEL)))
         CENTRO(1)=(X(1,conecNodaisElem(nen,NEL))+X(1,conecNodaisElem(1,NEL)))*AUX
         CENTRO(2)=(X(2,conecNodaisElem(nen,NEL))+X(2,conecNodaisElem(1,NEL)))*AUX
         DO J=1,NEN-1
            AUX = (X(1,conecNodaisElem(J,NEL))*X(2,conecNodaisElem(J+1,NEL))-  &
                 X(1,conecNodaisElem(J+1,NEL))*X(2,conecNodaisElem(J,NEL)))
            CENTRO(1) = CENTRO(1)+ &
                 (X(1,conecNodaisElem(J,NEL))+X(1,conecNodaisElem(J+1,NEL)))*AUX
            aux1 = (X(1,conecNodaisElem(J,NEL))+X(1,conecNodaisElem(J+1,NEL)))
            CENTRO(2) = CENTRO(2)+ &
                 (X(2,conecNodaisElem(J,NEL))+X(2,conecNodaisElem(J+1,NEL)))*AUX
         END DO
         CENTRO=CENTRO*AREA
         NBLOCK=0
         DO K=1,NELEMY
            YPOS=YX2+REAL(K-1)*CDY
            DO J=1,NELEMX
               NBLOCK=NBLOCK+1
               XPOS=YX1+REAL(J-1)*CDX
               IF(CENTRO(1).GE.XPOS)THEN
                  IF(CENTRO(1).LE.XPOS+CDX)THEN
                     IF(CENTRO(2).GE.YPOS)THEN
                        IF(CENTRO(2).LE.YPOS+CDY)THEN
                           PERMK(NEL) = PERM(NBLOCK)
                           exit
                        END IF
                     END IF
                  END IF
               END IF
            END DO
         END DO
      END DO
!
800   format(/,                            &
           '##############################',/, &
           'TAMANHO DO RESERVATORIO:',/,'Lx=', &
           f10.5,2x,'Ly=',f10.5,/,             &
           'MALHA:',/,'nx=',i7,' ny=',i7,/,    &
           '##############################')
!
890   format(/,                              &
           '##############################',/,   &
           'TAMANHO DO CAMPO ALEATORIO:',/,'Lx=',&
           f10.5,2x,'Ly=',f10.5,/,               &
           'MALHA:',/,'nx=',i7,' ny=',i7,/,      &
           '##############################')
!
891   format(/,                            &
           '##############################',/, &
           'Tamanho do elemento:',/,'hx=',     &
           f10.5,2x,'hy=',f10.5,/, &
           '##############################',/, &
           'Tamanho do bloco geologico:',/,'gx=', &
           f10.5,2x,'gy=',f10.5,/,&
           '##############################',/)

900   format(/,                            &
           '##############################',/, &
           'TAMANHO DO RESERVATORIO:',/,'Lx=', &
           f10.5,2x,'Ly=',f10.5,2x,'Lz=',f10.5,/,&
           'MALHA:',/,'nx=',i7,' ny=',i7,' nz=',i7,/,    &
           '##############################')
!
990   format(/,                              &
           '##############################',/,   &
           'TAMANHO DO CAMPO ALEATORIO:',/,'Lx=',&
           f10.5,2x,'Ly=',f10.5,2x,'Lz=',f10.5,/,&
           'MALHA:',/,'nx=',i7,' ny=',i7,' nz=',i7,/,&
           '##############################')
!
991   format(/,                            &
           '##############################',/, &
           'Tamanho do elemento:',/,'hx=',     &
           f10.5,2x,'hy=',f10.5,2x,'hz=',f10.5,/, &
           '##############################',/, &
           'Tamanho do bloco geologico:',/,'gx=', &
           f10.5,2x,'gy=',f10.5,2x,'gz=',f10.5,/,&
           '##############################',/)
!
100   FORMAT('POSICAO X:',F10.5,/,'POSICAO Y:',F10.5)
!
110   format('####################################',/, &
           'INCONSISTENCIA NO TAMANHO DOS BLOCOS',/, &
           'GEOLOGICOS: MALHA COMPUTACIONAL     ',/, &
           'MAIS GROSSEIRA QUE A MALHA GEOLOGICA',/, &
           '####################################',/)
!
111   format('####################################',/, &
           'INCONSISTENCIA NO TAMANHO DOS       ',/, &
           'DOMINIOS: O DOMINIO GEOLOGICO DEVE  ',/, &
           'SER MAIOR OU IGUAL AO DOMINIO       ',/, &
           '####################################',/)
!
      return
!
    END SUBROUTINE MAPEIADOMINIO
!
!==========================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===========================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
