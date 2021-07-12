import sys
import os
from struct import *
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import time

#change the style of the plot (there are a lot more)
plt.style.use("fivethirtyeight")

#Will use matplotlib in order to plot the average execution times
#of the parallel codes.
#The codes were run with 1 to 20 threads.

#Store values of the 

#Calculate the average, using the .txt files
#Pretty sure would be better of using csv...

#Returns the position of val, in array arr
#If val not in arr returns -1
def valuePos(arr, val):
    if len(arr) == 0:
        return -1
    else:
        pos = 0
        for x in arr:
            if int(x) == int(val):
                return pos
            else:
                pos += 1
    return -1


#Scans the whole document and returns an array, of size equal to the amount of measurements
#that have different number of threads
#Each position contains the number of threads
#for mpi_omp it refers to size not num_threads
def get_xaxis(openf):
    xaxis = []
    lines = openf.readlines()
    for line in lines:
        #get the number of threads (concerning the current line)
        num_threads = line.split(" ")[2]
        num_threads = num_threads.split("\n")[0]
        pos = 0
        if  valuePos(xaxis,num_threads) == -1:
            xaxis.append(int(num_threads))
    xaxis.sort()
    return xaxis

#Calculates the average values of the measurements in openf
#xaxis MUST be an array returned by get_xaxis(openf)
def calc_aver(openf, xaxis):
    yaxis = [0 for i in range(len(xaxis))]
    counter = [0 for i in range(len(xaxis))]
    openf.seek(0)
    lines = openf.readlines()
    for line in lines:
        words = line.split(" ")
        pos = valuePos(xaxis,int(words[2].split("\n")[0]))
        yaxis[pos] += float(words[0])
        counter[pos] += 1

    #Round to 3 decimal points since each sample has a precision of 3 decimal points
    #Pretty sure need to round before each addition to avoid cumulative errors, but oh well
    for i in range(len(yaxis)):
        yaxis[i] = round(round(yaxis[i], 3) / counter[i], 3)
    return yaxis


#Plots the measurements contained in openf
#Calls get_xaxis() and get_yaxis()
def quick_plot(program_name, openf):
    #configure the window of the plot I am creating
    plt.figure(program_name, (9,7))
    #Get x,y values
    xaxis = get_xaxis(openf)
    yaxis = calc_aver(openf, xaxis)

    #Check if dimensions are ok (DELETEME)
    if(len(xaxis) != len(yaxis)):
        print("Lengths do not match")
        print("length of xaxis: ", len(xaxis))
        print("length of yaxis: ", len(yaxis))

    #create a simple line plot
    plt.plot(xaxis, yaxis, marker="o", label="our code")
    xl = "Number of threads"
    if program_name.find("mpi") != -1:
        xl = "MPI_COMM_WORLD size"
    #add labels
    plt.title(f"Average running times of parallel code ({program_name})")
    plt.xlabel(xl)
    plt.ylabel("Running time (sec)")

    #X - axis ticks must only be integers (19.5 threads cannot exist)
    #plt.xticks( range(0,21,2) )

    #TODO: add optimal line - must be linear + legend
    #Do I need to study the roofline model? -> probably

    #Adding a grid
    plt.grid(True, linestyle='--')
    #Show and save the plot
    plt.savefig("fig_"+program_name)
    plt.show()
    return plt


#ploting the speedups
def speedUp_plot(program_name, openf):
    #configure the window of the plot I am creating
    plt.figure(program_name, (9,7))
    #Get x,y values
    xaxis = get_xaxis(openf)
    y = calc_aver(openf, xaxis)
    su = [0 for i in range(len(xaxis))]
    yaxis = [0 for i in range(len(xaxis))]
    #speedup =  sequential time / parallel time
    for i in range(len(xaxis)):
        su[i] = y[0] / y[i]

    yaxis = su
    #Check if dimensions are ok (DELETEME)
    if(len(xaxis) != len(yaxis)):
        print("Lengths do not match")
        print("length of xaxis: ", len(xaxis))
        print("length of yaxis: ", len(yaxis))

    #create a simple line plot
    plt.plot(xaxis, yaxis, marker="o", label="our code")
    xl = "Number of threads"
    if program_name.find("mpi") != -1:
        xl = "MPI_COMM_WORLD size"
    #add labels
    plt.title(f"Parallel code speedup ({program_name})")
    plt.xlabel(xl)
    plt.ylabel("Speed up")

    #X - axis ticks must only be integers (19.5 threads cannot exist)
    #plt.xticks( range(0,21,2) )

    #TODO: add optimal line - must be linear + legend
    #Do I need to study the roofline model? -> probably

    #Adding a grid
    plt.grid(True, linestyle='--')
    #Show and save the plot
    plt.savefig("fig_"+program_name+"_speedup")
    plt.show()
    return plt


def hybrid_plot(program_name, openf):
    plt.figure(figsize=(15,8), dpi=80)

    #yaxis
    yaxis = []
    labels = []
    lines = openf.readlines()
    for line in lines:
        # get the number of running (concerning the current line)
        time = line.split(" ")[0]
        pos = 0
        yaxis.append(float(time))
    #yaxis.sort()
    #return yaxis
    
    #xaxis
    x = list(range(0,48))
    for i in range(1, 7):
        for j in range(1, 9):
            lab_str = "(" + str(i) + "," + str(j) + ")"
            labels.append(lab_str)

    plt.plot(labels, yaxis)
    plt.xticks(x, labels, rotation = 'vertical')
    plt.title(f"Running times of parallel code ({program_name})")
    plt.xlabel("(size, threads)")
    plt.ylabel("Running time (sec)")

    plt.grid(True, linestyle='--')
    #Show and save the plot
    plt.savefig("fig_"+program_name)
    plt.show()
    return plt

def hybrid_speedUp_plot(program_name, openf):
    plt.figure(figsize=(15,8), dpi=80)

    y = []
    labels = []
    lines = openf.readlines()
    for line in lines:
        #get the number of running (concerning the current line)
        time = line.split(" ")[0]
        pos = 0
        y.append(float(time))
    
    su = [0 for i in range(len(y))]
    
    #speedup =  sequential time / parallel time
    for i in range(0, 48):
        su[i] = y[0] / y[i]

    y = su
    x = list(range(0,48))
    for i in range(1, 7):
        for j in range(1, 9):
            lab_str = "(" + str(i) + "," + str(j) + ")"
            labels.append(lab_str)

    plt.plot(labels, y)
    plt.xticks(x, labels, rotation = 'vertical')
    plt.title(f"Speedup of parallel code ({program_name})")
    plt.xlabel("(size, threads)")
    plt.ylabel("Speed Up")
    plt.grid(True, linestyle='--')

    #Show and save the plot
    plt.savefig("fig_"+program_name+"_speedup")
    plt.show()
    return plt

#Get the part of the filename that differentiates the files
def get_dif_fname(exec_filename):
    spltf = exec_filename.split("_")
    if len(spltf) == 3:
        return spltf[2]
    else:
        return spltf[2]+"_"+spltf[3]


#TODO: make if-else safe
if(sys.argv[1] == "all"): #plot all
    n = 4
    pids = list(np.zeros(n))
    for i in range(n):
        pids[i] = os.fork()
        if pids[i] < 0:
            print("Error on fork!")
            exit(1)
        elif pids[i] == 0:
            if i == 0:
                os.system("python3 graph_benchmark.py multistart_hooke_omp")
            elif i == 1:
                os.system("python3 graph_benchmark.py multistart_hooke_omp_tasks")
            elif i == 2:
                os.system("python3 graph_benchmark.py multistart_hooke_mpi")
            elif i == 3:
                os.system("python3 graph_benchmark.py multistart_hooke_mpi_omp")
            exit(0)
    
    while n > 0:
        status = os.waitpid(pids[n-1], 0)
        n-=1
else:
    dif = get_dif_fname(sys.argv[1])
    #open the correct res txt file, using the differentiating part of the filename provided
    f = open("res_" + dif + ".txt", "r")
    if dif.find("omp") != -1 and dif.find("mpi") != -1:
        hybrid_plot(program_name = dif, openf=f)
        f.seek(0)
        hybrid_speedUp_plot(program_name = dif, openf=f)
    else:
        quick_plot(program_name = dif, openf=f)
        f.seek(0)
        speedUp_plot(program_name = dif, openf=f)
    f.close()