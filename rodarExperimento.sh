dirExp=${1:-"exp05x02/"}
opcaoA=${2:-"1"}
opcaoB=${3:-"1"}
opcaoC=${4:-"1"}
numThreads=${5:-"1"}
sufixoTela=${6:-""};

dirBin=/prj/prjedlg/bidu/BTsimuladorMEFacademico/bin
dirBin=/prj/prjedlg/bidu/simuladorMEFacademicoBT-Patricia/bin
dirBin="$(pwd)/bin"

#if [ "$opcaoA" = "2" ]; then
#   opcaoB="2"
#fi

     listaSolvers=(zero Gauss Pardiso HYPRE PH PardisoMKL)
listaCompiladores=(zero gfortran ifort) # pgf90)
     listaOPTIMIZ=(zero "-g -O0" "-O4" "-fast")

 solver=${listaSolvers[opcaoA]}
     FC=${listaCompiladores[opcaoB]}
OPTIMIZ=${listaOPTIMIZ[opcaoC]}
maquina=$(hostname)
#maquina="altix-xe.hpc.lncc.br"

     lS=(z G P H PH P-mkl)
     lC=(z G I)  # P)
     lO=(z 0 4 f)

sufixoExec="${lS[opcaoA]}${lC[opcaoB]}${lO[opcaoC]}"

LIBPARDISO=""
#LIBPARDISO=MKL; sufixoExec="${sufixoExec}${LIBPARDISO}"

#echo "escolhas: $opcaoA	 $opcaoB	$opcaoC"
#echo "escolhas: $solver	 $FC 		$OPTIMIZ ...  $maquina"
#echo "sufixo  : $sufixoExec"

case $maquina in
 "petropolis")
   PARDISO_DIR="/usr/local/pardiso"
   HYPRE_DIR="/usr/local/hypre2.9" # em petropolis
;;
 "bidu-debian")
   PARDISO_DIR="/usr/local/lib/pardiso5.0"
   HYPRE_DIR="/usr/local/hypre-2.10.1/"
   export PARDISO_LIC_PATH=$PARDISO_DIR
   export  LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PARDISO_DIR
;;
# /hpc/pardiso5.0
#Architecture X86-64, 64-bit, icc/ifort 13.01 libpardiso500-INTEL1301-X86-64.so
#Architecture X86-64, 64-bit, gcc/gfortran 4.7.2 libpardiso500-GNU472-X86-64.so

#Architecture X86-64, 64-bit, icc/ifort 13.01     Linux MPI libpardiso500-MPI-INTEL1301-X86-64.so
#Architecture X86-64, 64-bit, mpicc/mpif77 4.7.2     Linux MPI libpardiso500-MPI-GNU472-X86-64.so
"altix-xe.hpc.lncc.br")
   PARDISO_DIR="/hpc/pardiso5.0"
   if [ "$FC" = "gfortran" ]; then
      comando="source /hpc/modulos/bash/gcc-4.7.sh";              echo +++ $comando; eval $comando
      comando="source /hpc/modulos/bash/hypre-2.9.0b.sh";         echo +++ $comando; eval $comando
      comando="source /hpc/modulos/bash/openmpi-gcc44-1.4.1.sh "; echo +++ $comando; eval $comando
      comando="source /hpc/modulos/bash/libblas.sh";              echo $comando; eval $comando
      HYPRE_DIR="/hpc/hypre-2.9.0b"         # altix-xe, gcc
      PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-MPI-GNU472-X86-64 -lblas -llapack -fopenmp -lpthread -lm" 
      PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-GNU472-X86-64     -lblas -llapack -fopenmp -lpthread -lm" 
   fi
   if [ "$FC" = "ifort" ]; then
      comando="source /hpc/modulos/bash/intel-cluster_studio_xe_2013.sh"; echo $comando; eval $comando
      comando="source /hpc/modulos/bash/hypre-2.9.0b-intel.sh";           echo $comando; eval $comando
      HYPRE_DIR="/hpc/hypre-2.9.0b-intel"   # altix-xe, intel
      PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-MPI-INTEL1301-X86-64 -lblas -llapack -fopenmp -lpthread -lm" 
      PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-INTEL1301-X86-64     -lblas -llapack -fopenmp -lpthread -lm" 
    if [ "$LIBPARDISO" = "MKL" ]; then
      PARDISOLNK="-mkl"; 
      sufixoExec="${sufixoExec}${PARDISOLNK}"
    fi
   fi
;;
 "no41.hpc.lncc.br"|no42.hpc.lncc.br|no43.hpc.lncc.br|no44.hpc.lncc.br)
  PARDISO_DIR="/hpc/pardiso5.0"
  if [ "$FC" = "gfortran" ]; then
     comando="source /hpc/modulos/bash/gcc-4.7.sh";          echo $comando; eval $comando
     comando="source /hpc/modulos/bash/hypre-2.9.0b-k20.sh"; echo $comando; eval $comando
     comando="source /hpc/modulos/bash/libblas-k20.sh";      echo $comando; eval $comando
  fi
  if [ "$FC" = "ifort" ]; then
     comando="source /hpc/modulos/bash/intel-cluster_studio_xe_2013.sh"; echo $comando; eval $comando
     comando="source /hpc/modulos/bash/hypre-2.9.0b-intel-k20.sh";       echo $comando; eval $comando
      fi
;;
 *)
 echo "..... erro .... em rodarExperimento.sh" 
 exit;;
esac

comando="LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PARDISO_DIR}:${HYPRE_DIR}";echo $comando; eval $comando

if [ $# == 1 ]; then
 nomeExecutavel=simulador${simulador}.exe
else
 nomeExecutavel=simulador${simulador}${sufixoExec}.exe
fi

arqTela="tela_${sufixoExec}_${numThreads}Thr${sufixoTela}";
echo +++ |tee $dirExp/${arqTela}.txt

echo +++ |tee -a $dirExp/${arqTela}.txt 
date     |tee -a $dirExp/${arqTela}.txt 
echo +++ |tee -a $dirExp/${arqTela}.txt 
echo +++  RODAR EXPERIMENTO COM AS SEGUINTES ESCOLHAS: |tee -a $dirExp/${arqTela}.txt
echo +++ |tee -a $dirExp/${arqTela}.txt 
echo +++  maquina:    $maquina                  |tee -a $dirExp/${arqTela}.txt
echo +++  compilador: $FC                       |tee -a $dirExp/${arqTela}.txt
echo +++  solver:     $solver                   |tee -a $dirExp/${arqTela}.txt
echo +++  executavel: ${dirBin}/$nomeExecutavel |tee -a $dirExp/${arqTela}.txt
echo +++  diretorio:  ${dirExp}                 |tee -a $dirExp/${arqTela}.txt

echo +++ |tee -a $dirExp/${arqTela}.txt
echo +++ |tee -a $dirExp/${arqTela}.txt
echo "+++ digite ENTER para executar o simulador "  |tee -a $dirExp/${arqTela}.txt
echo +++ |tee -a $dirExp/${arqTela}.txt
comando="(export OMP_NUM_THREADS=1;time mpirun -np 2 ${dirBin}/${nomeExecutavel}  |tee -a ${arqTela}.txt)";
comando="(export OMP_NUM_THREADS=1;time ${dirBin}/${nomeExecutavel} <resp |tee -a ${arqTela}.txt)";
comando="(export OMP_NUM_THREADS=1;time mpirun -np 1 ${dirBin}/${nomeExecutavel} <resp  |tee -a ${arqTela}.txt)";
comando="(cd $dirExp; export OMP_NUM_THREADS=1; time mpiexec -np 1 ${dirBin}/${nomeExecutavel} <resp|tee -a ${arqTela}.txt)";
comando="(export OMP_NUM_THREADS=$numThreads; cd $dirExp ; time ${dirBin}/${nomeExecutavel}   |tee -a ${arqTela}.txt)";
echo $comando |tee -a $dirExp/${arqTela}.txt
#read
formato="\"\t%E real,\t%U user,\t%S sys\" "
#echo $formato
#time -f $formato eval $comando 2> tempoMedido.txt
dataInicio=$(date)
eval $comando &> tempoMedido.txt
dataFinal=$(date)
echo +++ inicio da simulacao: $dataInicio |tee -a $dirExp/${arqTela}.txt
echo +++ fim .. da simulacao: $dataFinal  |tee -a $dirExp/${arqTela}.txt
cat tempoMedido.txt; rm tempoMedido.txt

