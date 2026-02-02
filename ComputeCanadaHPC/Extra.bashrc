#To request interactive nodes on HPC machines
alias JOB1='srun -q primary -N 1 -n 1 --mem=4G -t 2-0  --pty bash'
alias JOB4='srun -q primary -N 1 -n 4 --mem=16G -t 2-0  --pty bash'
alias JOB8='srun -q primary -N 1 -n 8 --mem=8G -t 2-0  --pty bash'
alias JOB12='srun -q primary -N 1 -n 12 --mem=12G -t 2-0  --pty bash'
alias JOB16='srun -q primary -N 1 -n 16 --mem=16G -t 2-0  --pty bash'
alias JOB24='srun -q primary -N 1 -n 24 --mem=24G -t 2-0  --pty bash'
alias JOB32='srun -q primary -N 1 -n 32 --mem=32G -t 2-0  --pty bash'

alias JOB1E='srun -q express -N 1 -n 1 --mem=8G  -t 4-0  --pty bash'
alias JOB4E='srun -q express -N 1 -n 4 --mem=16G -t 4-0  --pty bash'
alias JOB8E='srun -q express -N 1 -n 8 --mem=32G -t 4-0  --pty bash'

alias JOB8L='srun -q primary -N 1 -n 8 --mem=16G -t 4-0  --pty bash'
alias JOB16L='srun -q primary -N 1 -n 16 --mem=16G -t 4-0  --pty bash'
alias JOB32L='srun -q primary -N 1 -n 32 --mem=32G -t 4-0  --pty bash'
