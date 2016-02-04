#set -x
#SCRIPT PARA AUXILIAR NA EXECUCAO DE PROGRAMAS 
# DE SIMULACAO DE ESCOAMENTOS EM RESERVATORIOS 
#
#
#AUTORES: bidu@lncc.br e tuane@lncc.br
# 06.2013
#
# foma de utilizacao . rodarSimulador.sh exp10x10x30

#### EXPERIMENTOS ####

# diretorio do experimento:
# dirPadrao ou lido de primeiro argumento da linha de comando
dirExp="$(pwd)"/${1:-${dirPadrao}} 

#rm -rf ${dirExp}/out/*.vtk
rm -rf ${dirExp}/initial_*
#rm -rf ${dirExp}/out/*.vtk
rm -rf ${dirExp}/out/*.dat
rm -rf ${dirExp}/conc/*amostra_0.dat
rm -rf ${dirExp}/desl/*amostra_0.dat
rm -rf ${dirExp}/pres/*amostra_0.dat
rm -rf ${dirExp}/prod/*amostra_0.dat
#rm -rf ${dirExp}/conc/*0.dat
#rm -rf ${dirExp}/desl/*0.dat
#rm -rf ${dirExp}/pres/*0.dat
#rm -rf ${dirExp}/prod/*0.dat
rm -rf ${dirExp}/fort.*
arqTela=${dirExp}/tela.txt

#### DEFINICAO DO NUMERO DE THREADS (OPENMP) E PROCESSOS (MPI) ####
ntPadrao=1
numThreads=${2:-${ntPadrao}}
export OMP_NUM_THREADS=${numThreads}


npPadrao=2
#numProcs=${3:-${npPadrao}}
numProcs=${npPadrao}
export NP=${numProcs} # atribuir valor padrao aa variavel NP

#ARGSP=" "
#echo "${ARGSP}"

#### DEFINICAO DO EXECUTAVEL ####
DIRBIN="/prj/prjedlg/tuane/Siger/simuladorGeocreep/bin"
DIRBIN=${3:-"./bin/"}


LOCAL=$(pwd)
NOMEEXECUTAVEL=simulador.exe

DIRBIN="/prj/prjedlg/bidu/Siger/simuladorGeocreep/ultimosPassos/bin"
DIRBIN="${LOCAL}/bin"

EXECUTAVEL=${DIRBIN}/${NOMEEXECUTAVEL}

#### definicao do comando a ser executado
comando="(export  OMP_NUM_THREADS=${numThreads} ; cd ${dirExp}; time  ${EXECUTAVEL})"

if [ -e ${EXECUTAVEL} ] 
then
  printf "\n diretorio do experimento.: %s\n" ${dirExp}  
  printf "\n nome do executavel.......: %s\n" ${EXECUTAVEL} 
  printf "\n numero de threads .......: %d\n" ${OMP_NUM_THREADS}
  printf "\n numero de processos......: %d\n" ${NP}
  printf "\n comando .................: %s\n" "${comando}"
  eval ${comando}  |tee  ${arqTela}
else
  printf "\n EXECUTAVEL NAO ENCONTRADO \n"
  printf "\n comando .................: %s\n" "${comando}"
fi
