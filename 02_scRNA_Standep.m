function [] = scRNA_Standep()
%%This script uses StanDep (https://github.com/LewisLabUCSD/StanDep) to transform
%%normalized scRNA-Seq expression data into core reactions used by FastCore. To this end
%%the normalized expression matrix ("PBMC_expression.tsv"), the IDs of genes as annotated in Recon 2.2 (HGNC-IDs in "PBMC_genes.tsv")
%%and the cell types ("PBMC_conditions.tsv") as well as the metabolic model ("recon22.mat") are required.
			    
initCobraToolbox(0)


%%Load data
expression = dlmread('/scRNA/STANDEP/Input/PBMC_expression.tsv');
fid = fopen('/scRNA/STANDEP/Input/PBMC_genes.tsv');
genes = textscan(fid,'%s');
fclose(fid)
fid = fopen('/scRNA/STANDEP/Input/PBMC_conditions.tsv');
tissues = textscan(fid,'%s');
fclose(fid)

%%Load metabolic model
load('/scRNA/STANDEP/Model/recon22.mat');

rnaData.gene=genes{1};
rnaData.valuebyTissue=expression';
rnaData.Tissue=tissues{1};

%Create StanDep data structure
modelData = getModelData(rnaData,model);

%Obtain enzyme types
spec = getSpecialistEnzymes(model);  
prom = getPromEnzymes(model);

%calculate enzyme expression
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);

%Further parameters
minLog=log10(min(enzymeData.value(enzymeData.value>0)))
maxLog=log10(max(enzymeData.value(enzymeData.value>0)))
edgeX=minLog:((maxLog-minLog)/7):maxLog;
%distMethod = 'euclidean'; % distance method
distMethod = 'chi2dist'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering

%Number of clusters
k=40;

%calculate clusters of enzyme expression
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,k,distMethod,linkageMethod);

%calculate core reaction sets
coreRxnMat = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],false,0,[1 1]); 

%store coreRxnMat
csvwrite('/scRNA/STANDEP/Output/PBMC_coreRxns.csv',coreRxnMat);
writetable(cell2table(model.rxns),'/scRNA/STANDEP/Output/PBMC_reactions_recon22.csv');

end
