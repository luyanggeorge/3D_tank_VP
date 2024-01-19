#$ -cwd -V
#$ -l h_rt=40:00:00
#$ -l node_type=40core-192G
#$ -pe ib 32
#$ -l h_vmem=12G
#$ -m be

module swap openmpi mvapich2
module add apptainer

mkdir -p /nobackup/$USER/.cache
mpiexec -n 32 singularity exec --env 'PATH=/home/firedrake/firedrake/bin:$PATH' -B /run -B /nobackup -B /nobackup/$USER/.cache:/home/firedrake/firedrake/.cache /home/home02/mmyl/firedrake_latest.sif python3 3DtankVP_TC4snapshots.py