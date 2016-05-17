opcaoA=${1:-"1"}
opcaoB=${2:-"1"}
opcaoC=${3:-"1"}

     listaSolvers=(zero Gauss Pardiso HYPRE P+H)
listaCompiladores=(zero gfortran ifort) # pgf90)
     listaOPTIMIZ=(zero "-g -O0" "-O4" "-fast")

 solver=${listaSolvers[opcaoA]}
     FC=${listaCompiladores[opcaoB]}
OPTIMIZ=${listaOPTIMIZ[opcaoC]}
maquina=$(hostname)

     lS=(z G P H PH)
     lC=(z G I)  # P)
     lO=(z 0 4 f)
sufixoExec="${lS[opcaoA]}${lC[opcaoB]}${lO[opcaoC]}"

LIBPARDISO=MKL

echo "escolhas: $opcaoA	 $opcaoB	$opcaoC"
echo "escolhas: $solver	 $FC 		$OPTIMIZ ...  $maquina"
echo "sufixo  : $sufixoExec"

dirFontes=fontes
dirBin=bin

OUTROSF="-DmostrarTempos -Ddebug"
OUTROSF="-DmostrarTempos        "

     listaFontes=(${dirFontes}/variaveisGlobais.F90  \
\
        ${dirFontes}/malha.F90 \
        ${dirFontes}/mInputReader.F90 \
        ${dirFontes}/utilitarios.F90 \
        ${dirFontes}/leituraEscrita.F90 \
        ${dirFontes}/funcoesDeForma.F90 \

        ${dirFontes}/propGeoFisica.F90 \
\
        ${dirFontes}/solverGaussSkyline.F90 \
        ${dirFontes}/solverHypre.F90 \
        ${dirFontes}/solverPardiso.F90 \
\
        ${dirFontes}/leituraEscritaSimHidroGeoMec.F90 \
        ${dirFontes}/hidrodinamica.F90 \
        ${dirFontes}/geomecanica.F90 \
        ${dirFontes}/transporte.F90 \
        ${dirFontes}/mcmc.F90 \
\
        ${dirFontes}/driverGeocreep.F90 )


     i=0
     for f in ${listaFontes[*]} 
       do
           listaObjetos[i]=${dirBin}/$(basename ${f/F90/o/})
           ((i=i+1))
       done

if [ "$FC" = "ifort" ]; then 
   ARGINC="-w -module include";
   COMP="-openmp  -DwithOMP"
   LOMP="-openmp"
fi

if [ "$FC" = "gfortran" ]; then
   ARGINC="-J include";
   COMP="-fopenmp  -DwithOMP"
   LOMP="-fopenmp"
fi

case $maquina in
"petropolis.apos3.lncc.br")
    PARDISO_DIR="/usr/local/pardiso"
    HYPRE_DIR="/usr/local/hypre2.9" # em petropolis
;;

"bidu-debian" )
    HYPRE_DIR="/usr/local/hypre-2.10.1/"
    PARDISO_DIR="/usr/local/lib/pardiso5.0/"
    export PARDISO_LIC_PATH=/usr/local/lib/pardiso5.0/
    export  LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/pardiso5.0/
;;
# /hpc/pardiso5.0
#Architecture X86-64, 64-bit, icc/ifort 13.01 libpardiso500-INTEL1301-X86-64.so
#Architecture X86-64, 64-bit, gcc/gfortran 4.7.2 libpardiso500-GNU472-X86-64.so

#Architecture X86-64, 64-bit, icc/ifort 13.01     Linux MPI libpardiso500-MPI-INTEL1301-X86-64.so
#Architecture X86-64, 64-bit, mpicc/mpif77 4.7.2     Linux MPI libpardiso500-MPI-GNU472-X86-64.so
"altix-xe.hpc.lncc.br")
    PARDISO_DIR="/hpc/pardiso5.0"
   if [ "$FC" = "gfortran" ]; then
      comando="source /hpc/modulos/bash/gcc-4.7.sh";       echo +++ $comando; eval $comando
      comando="source /hpc/modulos/bash/hypre-2.9.0b.sh";  echo +++ $comando; eval $comando
      comando="source /hpc/modulos/bash/openmpi-gcc44-1.4.1.sh ";   echo +++ $comando; eval $comando
      HYPRE_DIR="/hpc/hypre-2.9.0b"         # altix-xe, gcc
   fi
   if [ "$FC" = "ifort" ]; then
      comando="source /hpc/modulos/bash/intel-cluster_studio_xe_2013.sh"; echo $comando; eval $comando
      comando="source /hpc/modulos/bash/hypre-2.9.0b-intel.sh";           echo $comando; eval $comando
      HYPRE_DIR="/hpc/hypre-2.9.0b-intel"   # altix-xe, intel
   fi
;;

"no41.hpc.lncc.br")
     LAPACK_DIR=/hpc/lapack-3.6.0-k20
     PARDISO_DIR="/hpc/pardiso5.0"
  if [ "$FC" = "gfortran" ]; then
     comando="source /hpc/modulos/bash/gcc-4.7.sh";          echo $comando; eval $comando
     comando="source /hpc/modulos/bash/hypre-2.9.0b-k20.sh"; echo $comando; eval $comando
     comando="source /hpc/modulos/bash/libblas-k20.sh"     ; echo $comando; eval $comando
     HYPRE_DIR="/hpc/hypre-2.9.0b-babel-k20"  # k20, gcc 
  fi
  if [ "$FC" = "ifort" ]; then
     comando="source /hpc/modulos/bash/intel-cluster_studio_xe_2013.sh"; echo $comando; eval $comando
     comando="source /hpc/modulos/bash/hypre-2.9.0b-intel-k20.sh";       echo $comando; eval $comando
     HYPRE_DIR="/hpc/hypre-2.9.0b-intel-k20"  # k20, intel 
  fi
;;
*)
  echo "..... erro .... em rodarExperimento.sh" 
  exit;;
esac

case $solver in
"Gauss" )
   ppSolver="-DSkyline"; LIBS="" ;  
;;

"Pardiso")
  if [ "$FC" = "gfortran" ]; then
    FC=gfortran
    comando="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/hpc/pardiso5.0/"; eval $comando
    PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-MPI-GNU472-X86-64 -lblas -llapack -fopenmp -lpthread -lm" 
    PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-GNU472-X86-64 -L${LAPACK_DIR}     -lblas -llapack -fopenmp -lpthread -lm" 
  fi
  if [ "$FC" = "ifort" ]; then
    FC=ifort
    if [ "$LIBPARDISO" = "MKL" ]; then
       PARDISOLNK="-mkl"; 
    else
       PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-MPI-INTEL1301-X86-64 -lblas -llapack -fopenmp -lpthread -lm" 
       PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-INTEL1301-X86-64  -L${LAPACK_DIR}    -lblas -llapack -fopenmp -lpthread -lm" 
    fi
  fi
   comando="LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PARDISO_DIR}"; eval $comando
   ppSolver="-DwithPardiso -Dwithcrs -DSkyline  $COMP ";  
   LIBS="${PARDISOLNK}"; 
;;

"HYPRE")
  if [ "$FC" = "gfortran" ]; then
    FC=mpif90
  fi
  if [ "$FC" = "ifort" ]; then
    FC=mpiifort
  fi
  ppSolver="-DwithHYPRE -Dwithcrs $COMP"; 
  HYPRELNK="-L${HYPRE_DIR}/lib -lHYPRE " 
  HYPRE_INC="${HYPRE_DIR}/include/"
  LIBS="${HYPRELNK}"; 
  comando="LD_LIBRARY_PATH=$LD_LIBRARY_PATH::${HYPRE_DIR}"; eval $comando
;;

"P+H")
  ppSolver="-DwithHYPRE -DwithPardiso -Dwithcrs $COMP";
  HYPRELNK="-L${HYPRE_DIR}/lib -lHYPRE " 
  HYPRE_INC="${HYPRE_DIR}/include/"
  if [ "$FC" = "gfortran" ]; then
    FC=mpif90
    PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-GNU472-X86-64 -L${LAPACK_DIR} -lblas -llapack -fopenmp -lpthread -lm"
    LIBS="${HYPRELNK} ${PARDISOLNK} "; 
  fi
  if [ "$FC" = "ifort" ]; then
    FC=mpiifort
    PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-INTEL1301-X86-64     -lblas -llapack -fopenmp -lpthread -lm"
    LIBS="${HYPRELNK} -mkl"; 
  fi

  if [ "$LIBPARDISO" = "MKL" ]; then
       PARDISOLNK="-mkl";
  else
       PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-MPI-INTEL1301-X86-64 -lblas -llapack -fopenmp -lpthread -lm"
  fi
   ppSolver="-DwithHYPRE -DwithPardiso -Dwithcrs -DSkyline $COMP ";
   comando="LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PARDISO_DIR}:${HYPRE_DIR}"; eval $comando
;;
*)
  echo "...erro na escolha do solver....em rodarExperimento.sh" 
  exit;;
esac

#echo +++  maquina FC solver ppSolver ARGINC COMP LOMP LIBS
comando="($FC --version)"; echo $comando; eval $comando

if [ $# == 0 ]; then
nomeExecutavel=simulador${simulador}.exe
else 
nomeExecutavel=simulador${simulador}${sufixoExec}.exe
fi
#echo "+++ digite Enter para gerar o executavel: $nomeExecutavel "
#read
comando="rm -rf $dirBin/$nomeExecutavel"
echo $comando
eval $comando


for i in $(seq 1 ${#listaFontes[*]})
 do
   FFLAGS="${OPTIMIZ}  ${OUTROSF}"
   FFLAGS="${OPTIMIZ}  ${OUTROSF} ${COMP}"
   comando="${FC} -c ${ARGINC} ${FFLAGS} ${listaFontes[i-1]} -o ${listaObjetos[i-1]}  " ; 

if [ "${listaFontes[i-1]}" = "${dirFontes}/solverHypre.F90" ]; then
   comando="${FC} -c ${ARGINC} ${ppSolver} ${listaFontes[i-1]} -o ${listaObjetos[i-1]}" ; 
fi

if [ ".${listaFontes[i-1]}." = ".${dirFontes}/solverPardiso.F90." ]; then
   comando="${FC} -c ${ARGINC} ${ppSolver} ${listaFontes[i-1]} -o ${listaObjetos[i-1]}" ; 
fi

if [ "${listaFontes[i-1]}" = "${dirFontes}/driverGeocreep.F90" ]; then
   comando="${FC} -c ${ARGINC} ${ppSolver}  ${FFLAGS} ${listaFontes[i-1]} -o ${listaObjetos[i-1]}" ; 
fi
   echo $comando;  eval $comando
done

LFLAGS="${LOMP} ${OPTIMIZ}"
comando="${FC} ${LFLAGS} -o ${dirBin}/${nomeExecutavel} ${listaObjetos[*]} ${LIBS} ";
echo $comando; eval $comando


echo +++
echo +++ executavel criado
ls -ltr ${dirBin}/${nomeExecutavel}
echo +++

#rm $dirBin/*.o include/*

#echo +++
#comando="./rodarExperimento.sh  $opcaoA $opcaoB $opcaoC exp05x02_DS"
#echo $comando; 
#echo +++
#eval $comando
 

