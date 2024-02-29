clearvars
clc

load CO_to_PHB_model
%% getting data obtained with eQuilibrator

%%% Get the metabolite concentrations calculated by MDF analysis
file_name = 'compound_concentrations.txt';
delimiter = "	";
concentrations_table = readtable(file_name, 'Delimiter', delimiter);
%re-ordering the metabolites in the same order of the stoichiometric model
[~, idx] = ismember(concentrations_table.Compound, model.mets);
[~, sortorder] = sort(idx);
concentrations_table = concentrations_table(sortorder,:);
list_mets = table2array(concentrations_table(:,1));
concs_MDF = table2array(concentrations_table(:,2));

%%% Get the standard deltaG of the reactions participating in the network
file_name = 'standard_dGs.txt';
delimiter = "	";
dGs_table = readtable(file_name, 'Delimiter', delimiter);
%re-ordering the reactions in the same order of the stoichiometric model
[~, idx] = ismember(dGs_table.reactions, model.rxns);
[~, sortorder] = sort(idx);
dGs_table = dGs_table(sortorder,:);
list_rxns = table2array(dGs_table(:,1));
list_dGs = table2array(dGs_table(:,2));
R = 0.008314; T = 30 + 273.15;
%calculating equilibrium constants from standard dGs 
list_Keqs = exp((list_dGs.*-1)/(R*T));

%% Reading enzyme properties
protein_table = readtable('enzyme_properties.txt');

%% Creating the matrices to be filled
KMS=zeros(length(model.mets),length(model.rxns));
stoich_coef_substrates=zeros(length(model.mets),length(model.rxns));
KMP=zeros(length(model.mets),length(model.rxns));
stoich_coef_products=zeros(length(model.mets),length(model.rxns));
kcatR = zeros(length(model.rxns),1);
kcatF = zeros(length(model.rxns),1);

%% Making KM equal to MDF concentrations and calculating kcatR

for r=1:length(model.rxns)

%     printRxnFormula(model,model.rxns(r))

    for m=1:length(model.mets)

        if model.S(m,r)<0;

            KMS(m,r)=(concs_MDF(m)/10)*1000;%key step to decide if saturated or sub-saturated
            stoich_coef_substrates(m,r)=abs(model.S(m,r));

            if strcmp(model.metNames{m}, 'H2O')
                KMS(m, r) = 1;
                stoich_coef_substrates(m, r) = abs(model.S(m, r));
            end

        end

        if model.S(m,r)>0;
            KMP(m,r)=concs_MDF(m)*1000;
            stoich_coef_products(m,r)=abs(model.S(m,r));

            if strcmp(model.metNames{m}, 'H2O')
                KMP(m, r) = 1;
                stoich_coef_products(m, r) = abs(model.S(m, r));
            end
        end
    
    end
  
  %for each reaction, get the KM values from the matrix KMS 
  A = KMS(:,r);
  valuesKMS = A(A ~= 0);
  C = stoich_coef_substrates(:,r);
  valuescoef_substrates = C(C ~= 0);
  %Each KMS value is powered to its corresponding stoichiometric coefficient
  adjusted_substrates = valuesKMS .^ valuescoef_substrates;

  B = KMP(:,r);
  valuesKMP = B(B ~= 0);
  order2 = stoich_coef_products(:,r);
  valuescoef_products = order2(order2 ~= 0);
  %Each KMP value is powered to its corresponding stoichiometric coefficient
  adjusted_products = valuesKMP .^ valuescoef_products;

% Specify the reaction name you're interested in
target_reaction_name = model.rxns(r);  % Replace with the desired reaction name

% Find the row index for the specified reaction name
row_index = find(strcmp(protein_table.Reaction_Name, target_reaction_name));

% Get the corresponding kcat_meta_data value
%kcatF(r) = protein_table.kcat_meta_data(row_index);
% Using average kcat for central metabolic pathways in prokaryotes
kcatF(r) = 36.72;
% kcat reverse are calculated using the Haldane relationship 
kcatR(r)=(kcatF(r)*prod(adjusted_products))/(list_Keqs(r)*prod(adjusted_substrates));
       
end

%% exporting the consistent kinetic parameters in the required format
combined_matrix = KMS+KMP;

[metabolites,reactions]=size(combined_matrix);
chain = [];
chain2 = [];
KMs = cell(reactions,1);

for g=1:reactions

    for u=1:metabolites
            if combined_matrix(u,g)~0;
                chain = strcat(model.mets(u),':',num2str(combined_matrix(u,g)));
                chain2 = strcat(chain2,{' '},chain);
            end
            
    end
    KMs(g) = chain2;
    chain2 = [];
end

all_KMS = cell2table(KMs);

T1 = table(model.rxns,kcatF,kcatR);
T1 = [T1, all_KMS];

T1.Properties.VariableNames = {'Reaction Name' 'kcatf (1/s)' 'kcatr (1/s)' 'KM (mM)'};
[~, idx] = ismember(T1.("Reaction Name"), order);
[~, sortorder] = sort(idx);
T1 = T1(sortorder,:);
writetable(T1,'consistent_kinetic_parameters.txt','Delimiter',"\t");

[~, idx] = ismember(protein_table.("Reaction_Name"), order);
[~, sortorder] = sort(idx);
protein_table = protein_table(sortorder,:);

T2 = T0.("Reaction Formula");
%T2.Properties.VariableNames = {'Reaction Formula'}; 
T3 = table(fluxes');
T3.Properties.VariableNames = {'Relative Flux'}; 
T4 = table(protein_table.MW);
T4.Properties.VariableNames = {'MWe(Da)'};

T5 = [T2,T3,T1,T4];
T5.Properties.VariableNames{'Var1'} = 'Reaction Formula';
T5.Properties.VariableNames{'kcatf (1/s)'} = 'kcrf(1/s)';
T5.Properties.VariableNames{'kcatr (1/s)'} = 'kcrr(1/s)';
T5.Properties.VariableNames{'KM (mM)'} = 'kM(mM)';

% Save the table as a CSV file
writetable(T5,'CO_to_PHB_ECM.csv');