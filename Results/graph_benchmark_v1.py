import sys
import os
from struct import *
from matplotlib import pyplot as plt
import numpy as np

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


#Get the part of the filename that differentiates the files
def get_dif_fname(exec_filename):
    spltf = exec_filename.split("_")
    if len(spltf) == 3:
        return spltf[2]
    else:
        return spltf[2]+"_"+spltf[3]


#TODO: make if-else safe
if(sys.argv[1] == "all"): #plot all
    #TODO: fork 3 times and wait for all pids
    pid = os.fork()
    if(pid == 0): #child
        os.system("python3 graph_benchmark.py multistart_hooke_omp")
    else: #parent
        os.system("python3 graph_benchmark.py multistart_hooke_omp_tasks")
        os.waitpid(pid,0)
else:
    dif = get_dif_fname(sys.argv[1])
    #open the correct res txt file, using the differentiating part of the filename provided
    f = open("res_" + dif + ".txt", "r")
    quick_plot(program_name = dif, openf=f)
    f.close()