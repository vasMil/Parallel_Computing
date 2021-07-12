CC=gcc
CFLAGS=-Wall -fopenmp -O3

#TODO: add the following implementations: multistart_hooke_omp multistart_hooke_omp_tasks multistart_hooke_mpi multistart_hooke_mpi_omp

all: benchmark multistart_hooke_seq multistart_hooke_omp multistart_hooke_omp_tasks multistart_hooke_mpi multistart_hooke_mpi_omp

hooke: multistart_hooke_seq.c
	$(CC) $(CFLAGS) -o multistart_hooke_seq multistart_hooke_seq.c -lm

multistart_hooke_omp: multistart_hooke_omp.c
	$(CC) $(CFLAGS) -o multistart_hooke_omp multistart_hooke_omp.c -lm

multistart_hooke_omp_tasks: multistart_hooke_omp_tasks.c
	$(CC) $(CFLAGS) -o multistart_hooke_omp_tasks multistart_hooke_omp_tasks.c -lm

multistart_hooke_mpi: multistart_hooke_mpi.c
	mpicc $(CFLAGS) -o multistart_hooke_mpi multistart_hooke_mpi.c -lm

multistart_hooke_mpi_omp: multistart_hooke_mpi_omp.c
	mpicc $(CFLAGS) -o multistart_hooke_mpi_omp multistart_hooke_mpi_omp.c -lm

benchmark: benchmark.c
	gcc -o benchmark benchmark.c

clean:
	rm -f multistart_hooke_seq multistart_hooke_omp multistart_hooke_omp_tasks multistart_hooke_mpi multistart_hooke_mpi_omp benchmark
