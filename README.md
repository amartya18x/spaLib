#SPALIB#
##Compile the OenBLAS Library.
Go into the folders
Do a Make
##Compile the SuperLU Library.
Read the README in the SuperLU directory and compile.
It is currently not used but we aim to use it soon.
##Compile the Armadillo Library.
Go the folder and read the README
##Go in the individual folders of the GradeS, OMP, OMPR and COSAMP.
Open the Makefile in there
Make the appropriate changes to the library paths
Type make
##To do an experiemnt
Create a folder in ~/data for each dataset
Create a file X that will store the measurement matrix
Create a file theta that will store the original signal.
Create a file y that will store the transformed signal.
Run ./<executable> ~/data/X ~/data/y ~/data/theta
##Generate Data for testing
You need to have MATLAB installed for this.
Go to IHT_Code
Edit the genAllData.m ro include the particular parameters you want to generate the data for.
Type matlab -r "run genAllData.m"
##Sanity Check experiment
Complete all the installations and making of the code.
###This is the output I have got on my computer.
COSAMP:../data/Samples_p700_s70_e0.12_o1_Ceps0.1/X:700:459:1.529507:17.262565:0.857143
COSAMP:../data/Samples_p700_s70_e0.15_o1_Ceps0.1/X:700:459:3.444802:18.845160:0.828571
COSAMP:../data/Samples_p700_s70_e0.1_o1_Ceps0.1/X:700:459:1.763942:13.697971:0.885714
COSAMP:../data/Samples_p700_s70_e0.2_o1_Ceps0.1/X:700:459:1.643382:23.422408:0.757143
GRADES:../data/Samples_p700_s70_e0.12_o1_Ceps0.1/X:700:459:0.008804:4.636733
GRADES:../data/Samples_p700_s70_e0.15_o1_Ceps0.1/X:700:459:0.011711:18.545763
GRADES:../data/Samples_p700_s70_e0.1_o1_Ceps0.1/X:700:459:0.011234:4.557762
GRADES:../data/Samples_p700_s70_e0.2_o1_Ceps0.1/X:700:459:0.012382:19.461333
OMPR:../data/Samples_p700_s70_e0.12_o1_Ceps0.1/X:700:459:1.684092:2.352459
OMPR:../data/Samples_p700_s70_e0.15_o1_Ceps0.1/X:700:459:1.725066:3.587574
OMPR:../data/Samples_p700_s70_e0.1_o1_Ceps0.1/X:700:459:0.672750:48.637012
OMPR:../data/Samples_p700_s70_e0.2_o1_Ceps0.1/X:700:459:2.931193:14.421245
OMP:../data/Samples_p700_s70_e0.12_o1_Ceps0.1/X:700:459:7.249897:2.726920
OMP:../data/Samples_p700_s70_e0.15_o1_Ceps0.1/X:700:459:6.080593:3.757105
OMP:../data/Samples_p700_s70_e0.1_o1_Ceps0.1/X:700:459:7.172752:3.456214
OMP:../data/Samples_p700_s70_e0.2_o1_Ceps0.1/X:700:459:6.726980:25.066389
