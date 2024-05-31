## Project Overview
This project is Carole Schöpfer's Bachelor thesis conducted under the supervision of Loïc Marrec and Claudia Bank in the Division of Theoretical Ecology and Evolution (THEE) at the Institute of Ecology and Evolution, University of Bern during the spring semester of 2024. 
The primary focus was on evolutionary rescue through different modes of migration (asymmetric and genotype-dependent) in spatially structured populations facing environmental deterioration. The simulations aimed to investigate the impact of these migration modes on rescue probability across various total migration rates.

The folders "Asymmetric Migration" and "Genotype-dependent Migration" contain MATLAB code for simulations and plotting, as well as raw data and corresponding plots used in the thesis. 
The folder "Visualization of Trajectories" contains R-code for generating and plotting single simulation trajectories, as well as plots of trajectories for different parameter values.


## Installation of Required Environment and Provided Code
For the simulations with different parameters I used the software MATLAB (23.2.0.2485118 (R2023b) Update 6). To run the code in MATLAB, a licences for this software is necessary. If your campus does not provide free access to MATAB, a licence can be bought. https://ch.mathworks.com/products/matlab/student.html
MATLAB can be used directly in a browser or as a desktop application. I used the MATLAB desktop environment to perform all simulations and plotting.
On request I can also provide the code as an R-Script.

1. Install MATLAB: [Download MATLAB](https://ch.mathworks.com/products/matlab/student.html)

2. Save the provided MATLAB code or raw data on your local machine.

3. Double-click on any provided MATLAB file to automatically start the software, or open MATLAB and access the folder containing the MATLAB files.

4. Note that MATLAB functions must be stored in the same folder as the code that calls them. (e.g. Function for asymmetric migration: "Gillespie_fct_deme" and the code recalling that function: "automated_analysis_deme")

5. Follow the instructions in the ## Usage section.


For R visualization of trajectories, R version 4.3.1 (2023-06-16 ucrt) was used. To install R:
1. [Download R](https://www.r-project.org)

2. Install R Studio: [Download R Studio](https://posit.co/download/rstudio-desktop/)

3. Install required packages ggplot2 and dplyr by running the following commands in the R console: 
   ```R
   install.packages("ggplot2")
   install.packages("dplyr")


## Usage and Configuration

# For MATLAB simulations:
For single simulations with rescue probability as output, open "testing_gillespie_deme.m" or "testing_gillespie_geno.m" in MATLAB.
Define the parameters as desired and let it run (green triangle in toolstrip). You will find the result in the Workspace environment.

For simulations with varying parameter values, open "automated_analysis_deme.m" or "automated_analysis_geno.m" in MATLAB. Define fixed parameters and desired ranges for testing parameters (total migration rate (T_mig) and the ratios of migration rates (alpha or beta)).
Note that depending on the chosen parameter values and parameter space the simulations may take a while. The code prints the progress of each parameter combination, so you can keep track on how far it has progressed.

# For R visualization:

Open "Analysis of Trajectories_genotype.R" or "Analysis of Trajectories_asymmetric.R" in R Studio.
For better overview, click on the small triangle in the first line of code at "Gillespie_fct_deme" to fold in the function.
At line 154 you can define the parameter values as desired. Adapt the titel of the plots accordingly in line 220. Specify the path and folder for the plots to be saved in line 225.
Run the whole code by selecting all (ctrl + a) or just the section you want to run (e.g. line 183 - 226) and press ctrl + Enter.


## Credits
This code is based on and extends the code used by Loïc Marrec in his paper "Quantifying the Impact of Genotype-Dependent Gene Flow on Mutation Fixation in Subdivided Populations" (2023). Hardware support was provided by the Division of Theoretical Ecology and Evolution (THEE) at the University of Bern.
Assistence in coding, especially translating R code into MATLAB code and vice versa, was given by AI (provided by OpenAI's ChatGPT).

## Contact Information
For inquiries, please contact Carole Schöpfer at carole.schoepfer@students.unibe.ch.