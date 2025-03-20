#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon March 14 2025

Subject :
Testing computation time of solution for multilayering in matter

@author: Alexandre Pellet, INRAE-IGE
"""
import subprocess
from collections import deque
from time import time
import datetime

#
cmd_getMPM=["cp mpm.cpp ../../src/mpm.cpp"]

# WITHOUT MODIFICATIONS
cmd_getUNIMAT=[
    """
    cp tools_UNIMAT.hpp ../../src/tools.hpp ;
    cd ../../build ;
    cmake .. ;
    make -j 15
    """
]
cmd_run_UNIMAT=[
    """
    cd ../../build/src ;
    ./mpm 'testing_UNIMAT'
    """
]
# WITH MODIFICATIONS
cmd_getMULTIMAT=[
    """
    cp tools_MULTIMAT.hpp ../../src/tools.hpp ;
    cd ../../build ;
    cmake .. ;
    make -j 15
    """
]
cmd_run_MULTIMAT=[
    """
    cd ../../build/src ;
    ./mpm 'testing_MULTIMAT'
    """
]
# ===================================================
# ||           --- Running the tests ---           ||
# ===================================================
date = datetime.datetime.now()
NKeepLines = 5
output_buffer_0 = deque(maxlen=NKeepLines)
output_buffer_1 = deque(maxlen=NKeepLines)
subprocess.run(cmd_getMPM, shell=True)

# UNIMATERIAL
subprocess.run(cmd_getUNIMAT, shell=True)
deb = time()
with subprocess.Popen(cmd_run_UNIMAT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True) as process:
    for line in process.stdout:
        print(line, end='')
        output_buffer_0.append(line)  # Add line to the buffer
    process.wait()
Dt0= time() - deb

# Yes modif
subprocess.run(cmd_getMULTIMAT, shell=True)
deb = time()
with subprocess.Popen(cmd_run_MULTIMAT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True) as process:
    for line in process.stdout:
        print(line, end='')
        output_buffer_1.append(line)  # Add line to the buffer
    process.wait()
Dt1 = time() - deb

# SHOW RESULTS
print("==========================================")
print("|| Without MULTIMATERIAL | %5.0f s      ||"%Dt0)
print("|| --------------------- | -------------||")
print("|| With MULTIMATERIAL    | %5.0f s      ||"%Dt1)
print("==========================================")

# print output info into file
export_file = "MULTIMATERIAL_RESULTS.txt"
with open(export_file, 'a') as file:
    file.write("\n")
    file.write("Test at time : %s\n" % date)
    file.write("=========================================================================\n")
    file.write("|| Without MULTIMATERIAL | Measured simulation time : %5.0f s          ||\n"%Dt0)
    for line in output_buffer_0:
        file.write("||                       | %-43s ||\n"%line[:-1])
    file.write("|| --------------------- | ------------------------------------------- ||\n")
    file.write("|| With MULTIMATERIAL    | Measured simulation time : %5.0f s          ||\n"%Dt1)
    for line in output_buffer_1:
        file.write("||                       | %-43s ||\n"%line[:-1])
    file.write("=========================================================================\n")


