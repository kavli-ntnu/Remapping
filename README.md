### **Realignment of grid modules \& place cell remapping**

#### 

This repository contains analysis code and processed data for the figures associated with the following paper: "Functional independence of entorhinal grid cell modules enables remapping in hippocampal place cells" by Lykken, Kanter, Nagelhus, Carpenter, Guardamagna, Moser \& Moser (2025). The manuscript is available on bioRxiv. DOI: https://www.biorxiv.org/content/10.1101/2025.09.24.677985v1



The repository contains:

-Code to plot the main figures, supplementary figures, and supplementary tables of the paper.

-Pre-processed data required to generate those figures.



Note that the code \& data provided here are intended only for reproducing the paper figures. The full dataset of neural recordings will be available at the time of publication.



The folder structure is as follows:

-data: pre-processed data needed to reproduce paper figures; organized by figure number: folder for figure/table (figure # or ed figure #; e.g., figure 1) > subfolder for figure panel (e.g., 1b) > script to produce figure panel (e.g., fig1b.m)

-helpers: extra code dependencies



**-------------------------------------------------------------------------------------------------------------------------------------------------------**



**System requirements** 



* System requirements for MATLAB 2020a (Windows):



Operating systems:

Windows 10 (version 1709 or higher)

Windows 7 Service Pack 1

Windows Server 2019

Windows Server 2016



Processors:

Minimum: any Intel or AMD x86-64 processor

Recommended: Any Intel or AMD x86-64 processors with four logical cores and AVX2 instruction set support



Disk:

Minimum: 3 GB of HDD space for MATLAB only; 5-8 GB for a typical installation

Recommended: An SSD is recommended

A full installation of all MathWorks products may take up to 31 GB of disk space



RAM:

Minimum: 4 GB

Recommended: 8 GB



Graphics: 

No specific graphics card is required.

Hardware accelerated graphics card supporting OpenGL 3.3 with 1GB GPU memory is recommended.

GPU acceleration using the Parallel Computing Toolbox requires a GPU that supports CUDA or newer.



* Software has been tested on:

MATLAB 2020a



* Non-standard hardware is not required.



**Installation guide**



* Instructions for MATLAB installation can be found here:

https://se.mathworks.com/help/install/ug/install-products-with-internet-connection.html



* Typical installation time: ~1-3 hours, but may vary depending on hardware, network delays, etc.



**Instructions for use**



Copy or clone repository to your local machine. (https://github.com/kavli-ntnu/Remapping)



To get started, change the MATLAB directory to the root directory of this library (where this README file is located).



To generate the figures, choose a subdirectory within the folder, e.g., figure 1.

Inside the folder are subfolders with code and pre-processed data that will generate all panels for each figure (very short run time).



Other scripts within subfolders include code and pre-processed data for main analysis procedures (e.g., generation of module crosscorrelograms). Longer run time (6-8 hours) is expected.





