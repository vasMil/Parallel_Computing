#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <string.h>

#define NUM_THR 50

void benchmark_me(char str[]){
    int i;
    for(i = 1; i <= NUM_THR; i++){
        if(strstr(str, "mpi") != NULL){
            if(i <= 6){
                char start[200] = "mpirun -n ";
                char int_char[5];
                sprintf(int_char, " %d", i);
                strcat(start, int_char);
                strcat(start, " ./");
                strcat(start, str);
                if(strstr(str, "omp") != NULL){
                    char num_thr[5];
                    char temp[sizeof(start)];
                    strcpy(temp, start);
                    for(int j = 1; j <= 8; j++){
                        strcpy(start, temp);
                        sprintf(num_thr, " %d", j);
                        strcat(start, num_thr);
                        system(start);
                        sleep(20);
                    }
                }
                else{
                    system(start);
                    sleep(20);
                }
            }
        }
        else{
            if( (i <= 26) || (i == 30) || (i == 40) || (i == 50) ){
                char start[50] = "./";
                char int_char[5];
                strcat(start, str);
                sprintf(int_char, " %d", i);
                strcat(start, int_char);
                system(start);
                sleep(20);
            }
        }
    }
}

int main(int argc, char *argv[]){
    int i;
    char exec_command[100] = "python3 graph_benchmark.py ";
    if(argc == 2){
        if(!strcmp(argv[1], "all")){ //just run benchmark to all (no plot)
            printf("I will benchmark all parallel codes\n");
            benchmark_me("multistart_hooke_omp");
            sleep(20);
            benchmark_me("multistart_hooke_omp_tasks");
            sleep(20);
            benchmark_me("multistart_hooke_mpi");
            sleep(20);
            benchmark_me("multistart_hooke_mpi_omp");
        }
        else{//only benchmark the specified file
            if(!strcmp(argv[1], "multistart_hooke_omp_tasks") || !strcmp(argv[1], "multistart_hooke_omp") 
                || !strcmp(argv[1], "multistart_hooke_mpi") || !strcmp(argv[1], "multistart_hooke_mpi_omp")){

                printf("I will benchmark %s\n",argv[1]);
                benchmark_me(argv[1]);
            }
            else{
                printf("This is not a valid input!\n");
                goto options; //Bad programming
            }
        }
    }
    else if(argc == 3){
        if(!strcmp(argv[2], "plot")){
            if(!strcmp(argv[1], "all")){
                printf("I will benchmark all and then plot them\n");
                benchmark_me("multistart_hooke_omp");
                sleep(20);
                benchmark_me("multistart_hooke_omp_tasks");
                sleep(20);
                benchmark_me("multistart_hooke_mpi");
                sleep(20);
                benchmark_me("multistart_hooke_mpi_omp");
                chdir("Results");
                strcat(exec_command, argv[1]);
                system(exec_command);
                chdir("..");
            }
            else{ //benchmark and plot the specified file
                if(!strcmp(argv[1], "multistart_hooke_omp_tasks") || !strcmp(argv[1], "multistart_hooke_omp") 
                    || !strcmp(argv[1], "multistart_hooke_mpi") || !strcmp(argv[1], "multistart_hooke_mpi_omp")){

                    printf("I will benchmark and plot %s\n", argv[1]);
                    benchmark_me(argv[1]);
                    chdir("Results");
                    strcat(exec_command, argv[1]);
                    system(exec_command);
                    chdir("..");
                }
                else{
                    printf("This is not a valid input!\n");
                    goto options;
                }
            }
        }
        else{
            goto options; //Bad programming
        }
    }
    else if(argc == 4 && !strcmp(argv[3], "only")){
        printf("I will plot %s\n", argv[1]);
        chdir("Results");
        strcat(exec_command, argv[1]);
        system(exec_command);
        chdir("..");
    }
    else{
        options:
        printf("\n");
        printf("\t\t\t\t   ---------\n");
        printf("\t\t\t\t   |OPTIONS|\n");
        printf("\t\t\t\t   ---------\n\n\n");
        printf("./benchmark filename [plot] [only]\n\n\n");
        printf("filename:\nmultistart_hooke_omp | multistart_hooke_omp_tasks | multistart_hooke_mpi_tasks | multistart_hooke_mpi_omp_tasks\n\n");
        printf("plot:\nIncluding this argument will to plot the measurements (saved in .txt files, in Results directory)\n\n");
        printf("only:\nIncluding this argument will result to only plot the current measurements (will not run the benchmark)\n\n");

        printf("\t\t\t\t    ------\n");
        printf("\t\t\t\t    |INFO|\n");
        printf("\t\t\t\t    ------\n\n\n");
        printf("Benchmark will run tests on the specified codes for number of threads 1 - 26, 30, 40, 50!\n");
        printf("If it is an mpi code it will run it for number of communicators 1 - 6!\n");
        printf("If you want to benchmark for a specific amount of threads you may run the executable directly, ");
        printf("with the first argument being the number of threads.\n"
                "(The measurement will be saved in the .txt file with the results)\n\n");
        exit(0);
    }
}