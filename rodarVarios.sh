 echo Nova simulacao ...................
 echo "solvers escolhidos no arquivo: exp120x16CREEP/SimulatorParam.inc"
 echo 
 sed '/solver_hidro/ {n; d}' exp120x16CREEP/SimulatorParam.inc -i       # remove linha com o solver para Hidromecannica
 sed '/solver_hidro/ a pardiso' exp120x16CREEP/SimulatorParam.inc -i    # inclui linha com o solver para Hidromecannica
 sed '/solver_geo/ {n; d}' exp120x16CREEP/SimulatorParam.inc -i         # remove linha com o solver para GEomecannica
 sed '/solver_geo/ a pardiso' exp120x16CREEP/SimulatorParam.inc -i      # inclui linha com o solver para GEomecannica
grep solver_  exp120x16CREEP/SimulatorParam.inc -A 2
comando="./rodarExperimento.sh exp120x16CREEP 1 4 1 2 PP |grep equacoes:"
echo $comando;
 echo "solvers escolhidos e lidos pelo simulador "
eval $comando;

 echo Nova simulacao ...................
 echo "solvers escolhidos no arquivo: exp120x16CREEP/SimulatorParam.inc"
 echo 
 sed '/solver_geo/ {n; d}' exp120x16CREEP/SimulatorParam.inc -i
 sed '/solver_geo/ a pardiso' exp120x16CREEP/SimulatorParam.inc -i
grep solver_  exp120x16CREEP/SimulatorParam.inc -A 2
comando="./rodarExperimento.sh exp120x16CREEP 1 4 1 2 GP |grep equacoes:"
#echo $comando;
#echo "solvers escolhidos e lidos pelo simulador "
#eval $comando;


 echo Nova simulacao ...................
 echo "solvers escolhidos no arquivo: exp120x16CREEP/SimulatorParam.inc"
 echo 
 sed '/solver_geo/ {n; d}' exp120x16CREEP/SimulatorParam.inc -i
 sed '/solver_geo/ a hypre' exp120x16CREEP/SimulatorParam.inc -i
grep solver_  exp120x16CREEP/SimulatorParam.inc -A 2
comando="./rodarExperimento.sh exp120x16CREEP 1 4 1 2 GH |grep equacoes:"
echo $comando;
 echo "solvers escolhidos e lidos pelo simulador "
eval $comando;


 echo Nova simulacao ...................
 echo "solvers escolhidos no arquivo: exp120x16CREEP/SimulatorParam.inc"
 echo 
 sed '/solver_hidro/ {n; d}' exp120x16CREEP/SimulatorParam.inc -i
 sed '/solver_hidro/ a pardiso' exp120x16CREEP/SimulatorParam.inc -i
grep solver_  exp120x16CREEP/SimulatorParam.inc -A 2
comando="./rodarExperimento.sh exp120x16CREEP 1 4 1 2 PH |grep equacoes:"
echo $comando;
 echo "solvers escolhidos e lidos pelo simulador "
eval $comando;

 sed '/solver_hidro/ {n; d}'    exp120x16CREEP/SimulatorParam.inc -i # remove linha com o solver para Hidromecannica
 sed '/solver_hidro/ a skyline' exp120x16CREEP/SimulatorParam.inc -i # inclui linha com o solver para Hidromecannica
 sed '/solver_geo/   {n; d}'    exp120x16CREEP/SimulatorParam.inc -i # remove linha com o solver para GEomecannica
 sed '/solver_geo/   a skyline' exp120x16CREEP/SimulatorParam.inc -i # inclui linha com o solver para GEomecannica

echo ..........
echo ...verificacao comparando a soma do  numero de passos do transposte
echo ...valor esperado 336
# for f in $(find  exp120x16CREEP/tela_*4Th*); 
 for f in $(ls -tr  exp120x16CREEP/tela_* |tail -4); 
  do
   echo arquivo tela: $f
   grep "Solver para a solucao"  $f -a
#   awk '/passos/{s=s+$4} END{print "num. passos do transporte=", s;} ' $f;
   grep passos $f -a | awk '/passos/{s=s+$4} END{print " soma do num. passos do transporte=", s;}'
 done

