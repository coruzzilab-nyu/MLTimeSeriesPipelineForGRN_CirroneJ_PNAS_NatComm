
function ParseRootShootTimeSeries2DFGInput()

%Quant 1e2 universe
univ_type = 'Quant_df2_1e2'
%TF_list:
tf_list_path = 'Quant.df2.1e2.TFs.ids';
%Gene list:
gene_list_path = 'Quant.df2.1e2.ids';


%0, 5, 10, 15, 20, 30, 45, 60, 90, 120
num_timepoints = 10;

deltaT =[5,5,5,5,10,15,15,30,30];

nComb = 9;



%Path dataset:
dataset_path = 'RootNormalizedCounts_woheader.txt';

dataset = importdata(dataset_path);

tf_list = importdata(tf_list_path);

gene_list = importdata(gene_list_path);


num_nonTFs = length(gene_list) - length(tf_list);

nonTFNames = cell(1, num_nonTFs);

j = 1;
for i = 1 : length(gene_list)

    gene = gene_list(i);
    if sum(ismember(tf_list, gene)) == 0
        nonTFNames(j) = gene;
        j = j+ 1;
    end
        
end

gene_list = horzcat(tf_list, nonTFNames);

geneNames = gene_list;
TFNames = tf_list;


name_input_ds = strcat('Arabidopsis_n',num2str(length(geneNames)),'g',num2str(length(TFNames)),'TF_',univ_type,'_root');


universe = zeros(length(geneNames), length(dataset.data(1,:)));

for i = 1 : length(geneNames)
    gene = geneNames(i);
    indx = find(strcmp(dataset.textdata,gene));
    universe(i,:) = dataset.data(indx, :);
end


dataKCl1 = zeros(length(geneNames), num_timepoints);
dataKCl2 = zeros(length(geneNames), num_timepoints);
dataKCl3 = zeros(length(geneNames), num_timepoints);

dataKNO31 = zeros(length(geneNames), num_timepoints);
dataKNO32 = zeros(length(geneNames), num_timepoints);
dataKNO33 = zeros(length(geneNames), num_timepoints);

%Also from universe matrix derive the 3 replicates

for i = 1 : length(geneNames)
    dataKCl1(i,:) = [universe(i,1),universe(i,4),universe(i,7),universe(i,10),universe(i,13),universe(i,16),universe(i,19),universe(i,22),universe(i,25),universe(i,28)];
    dataKCl2(i,:) = [universe(i,2),universe(i,5),universe(i,8),universe(i,11),universe(i,14),universe(i,17),universe(i,20),universe(i,23),universe(i,26),universe(i,29)];
    dataKCl3(i,:) = [universe(i,3),universe(i,6),universe(i,9),universe(i,12),universe(i,15),universe(i,18),universe(i,21),universe(i,24),universe(i,27),universe(i,30)];
    dataKNO31(i,:) = [universe(i,1), universe(i,31),universe(i,34),universe(i,37),universe(i,40),universe(i,43),universe(i,46),universe(i,49),universe(i,52),universe(i,55)];
    dataKNO32(i,:) = [universe(i,2), universe(i,32),universe(i,35),universe(i,38),universe(i,41),universe(i,44),universe(i,47),universe(i,50),universe(i,53),universe(i,56)];
    dataKNO33(i,:) = [universe(i,3), universe(i,33),universe(i,36),universe(i,39),universe(i,42),universe(i,45),universe(i,48),universe(i,51),universe(i,54),universe(i,57)];
end

KNO3=cell(9,1);
for k = 1:nComb
  KNO3{k} = zeros(length(geneNames), num_timepoints);
end

[nComb, KNO3] = GRN_CreateCombinations(geneNames,'', dataKNO31, dataKNO32, dataKNO33);


%vars: KNO3 (9x1 cell); TFNames(1xnum cell); dataKNO31(numgenesx10 double)
%dataKNO32(numgenesx10 double); dataKNO33(numgenesx10 double); 
%dataKnown(1x9cell)
%deltaT 5     5     5     5    10    15  15    30    30
%geneNames(1xnumgenes cell)
%ncomb 9
%nonTFNames(1xnumNonTFs cell)


dataKnown = cell(1, nComb);

for i = 1 : nComb
    dataKnown{1,i} = repmat(true,length(geneNames),num_timepoints);
end


save(strcat(name_input_ds,'.mat'), 'deltaT', 'nComb','geneNames', 'KNO3', 'TFNames', 'nonTFNames', 'dataKnown', 'dataKNO31', 'dataKNO32', 'dataKNO33');

end








