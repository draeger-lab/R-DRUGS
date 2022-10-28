initCobraToolbox(false);
changeCobraSolver('ibm_cplex', 'LP');
load('recon2_Influenza.mat');
epsilon=1e-4;
new_model=fastcc(Recon2, epsilon);
csvwrite('csvs/consistent_Influenza.csv', new_model); 
