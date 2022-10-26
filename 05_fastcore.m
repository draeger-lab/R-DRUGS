initCobraToolbox(false)
changeCobraSolver('ibm_cplex', 'LP');

 
files = dir('CoreReactsInfluenza/*.csv');
N=length(files);
Data=cell(1,N);
 

for i=1:N
Data{1,i} = csvread(strcat('Data/CoreReacts/', files(i).name));
end

load('Models/recon2_Influenza_consistent.mat')

for i=1:N
filename = (strcat('CoreModelsInfluenza/',files(i).name))
csvwrite(filename, fastcore(Data{1,i}, consistRecon2, 1e-4))
end 