#$ -cwd -V
#$ -l h_rt=40:00:00
#$ -l node_type=40core-768G
#$ -pe smp 40
#$ -l h_vmem=16G
#$ -l disk=10G
#$ -m be

module swap openmpi mvapich2
module add apptainer

mkdir $TMPDIR/.cache
mpiexec -n 32 singularity exec --env 'PATH=/home/firedrake/firedrake/bin:$PATH' -B /run -B /nobackup -B $TMPDIR/.cache:/home/firedrake/firedrake/.cache /home/home02/mmyl/firedrake_latest.sif python3 3DtankVP_TC4snapshots.py
