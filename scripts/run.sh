#!/bin/bash -l
## Nazwa zlecenia
#SBATCH -J 2DBosons_g=GG-N=NN
## Liczba alokowanych węzłów
#SBATCH -N 1
## Liczba zadań per węzeł (domyślnie jest to liczba alokowanych rdzeni na węźle)
#SBATCH --ntasks-per-node=1
## Ilość pamięci przypadającej na jeden rdzeń obliczeniowy (domyślnie 5GB na rdzeń)
#SBATCH --mem-per-cpu=7500MB
## Maksymalny czas trwania zlecenia (format HH:MM:SS)
#SBATCH --time=24:00:00 
## Nazwa grantu do rozliczenia zużycia zasobów
#SBATCH -A plgquantmol8-cpu
## Specyfikacja partycji
#SBATCH -p plgrid
## Plik ze standardowym wyjściem
#SBATCH --output="/net/people/plgrid/plgsuchorowski/2DBosons/tdl_run_g=GG-N=NN/output.out"
## Plik ze standardowym wyjściem błędów
#SBATCH --error="/net/people/plgrid/plgsuchorowski/2DBosons/tdl_run_g=GG-N=NN/error.err"



## Przejście do katalogu, z którego wywołany został sbatch
cd $SLURM_SUBMIT_DIR

srun /bin/hostname

## Dodaj biblioteki

# module load intel
#module load intel/16.0.1
module load gcc
module load cmake
module load eigen
module load fftw

# skompiluj program
./build.sh

# uruchom program
echo "Start: " $(date)
./simulation.sh
echo "Koniec: " $(date)
