#PBS -l walltime=240:30:00,nodes=2:ppn=12

cd /unicluster/home/matvey.kraposhin/run/Unicfdlab/github-salinityBoussinesqPimpleFoam/Tutorials-2.3.0/case-5


rm -rf log.new

mpirun -np 6 -npernode 3 -machinefile $PBS_NODEFILE salinityBoussinesqPimpleEqDyMFoam -parallel | tee -a log.new



