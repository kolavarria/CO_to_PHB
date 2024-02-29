clc
clearvars

%% Creating the metabolic model and performing FBA
model = createModel();

model=addReaction(model,'CODH_pump','CO + ADP + Pi -> H2 + CO2 + ATP');
model=addReaction(model,'Hydrogenase','H2 + NAD -> NADH');
model=addReaction(model,'PPiR5PK','Ri5P + PPi -> RiBP + Pi');
model=addReaction(model,'RiBPiso','RiBP -> RuBP');
model=addReaction(model,'RbuCO','RuBP + CO2 + H2O -> P3G + P3G');
model=addReaction(model,'PGK','P3G + ATP <=> BPG + ADP');
model=addReaction(model,'GAPDH','BPG + NADH <=> G3P + NAD + Pi');
model=addReaction(model,'TPI','G3P <=> DHAP');
model=addReaction(model,'FBPald','DHAP + G3P <=> FBP');
model=addReaction(model,'PPiPFK','FBP + Pi <=> F6P + PPi');
model=addReaction(model,'TALA','F6P + E4P <=> G3P + S7P');
model=addReaction(model,'Tkt','G3P + S7P <=> Ri5P + X5P');
model=addReaction(model,'RibI','Ru5P <=> Ri5P');
model=addReaction(model,'RibE','X5P <=> Ru5P');
model = addReaction(model,'PKT','F6P + Pi <=> E4P + AcP + H2O');
model = addReaction(model,'PTA','AcP + CoA <=> AcCoA + Pi');
model=addReaction(model,'Thio','AcCoA + AcCoA <=> AACoA + CoA');
model=addReaction(model,'AAR','AACoA + NADH <=> HBCoA + NAD');
model=addReaction(model,'PHBsyn','HBCoA + HB <=> PHB + CoA');

model=addReaction(model,'ATPase','ATP + H2O -> ADP + Pi');

model = addExchangeRxn(model,{'CO2'},0,50);
model = addExchangeRxn(model,{'H2O'},-50,50);
model = addExchangeRxn(model,{'HB'},-50,0);
model = addExchangeRxn(model,{'PHB'},0,50);
model = addExchangeRxn(model,{'CO'},-9,0);
clc;
model = changeObjective(model,'ATPase');

FBAsolution = optimizeCbModel(model,'max');
disp ('------- q-rates -------------')
printFluxVector(model, FBAsolution.x,1,1);

%% Creating the CSV file for MDF calculations using eQuilibrator

ActiveRxnsFormula = printRxnFormula(model,model.rxns,false);
% Define the string to be replaced and the replacement string
stringToReplace = {'  '};
replacementString = {' '};

% Iterate through the array of strings and replace the target string
for i = 1:numel(ActiveRxnsFormula)
    ActiveRxnsFormula{i} = strrep(ActiveRxnsFormula{i}, stringToReplace, replacementString);
end

% Define the string to be replaced and the replacement string
stringToReplace = {'->'};
replacementString = {'<=>'};

% Iterate through the array of strings and replace the target string
for i = 1:numel(ActiveRxnsFormula)
    ActiveRxnsFormula{i} = strrep(ActiveRxnsFormula{i}, stringToReplace, replacementString);
end

T0 = table(ActiveRxnsFormula,FBAsolution.x,model.rxns,'VariableNames',{'Reaction Formula' 'Relative Flux' 'Reaction Name'});

list_zero_flux = {};

for i=1:length(model.rxns)
    if FBAsolution.x(i)==0
       list_zero_flux = [list_zero_flux, model.rxns(i)];
    end
end

% Identify rows to remove
rows_to_remove = ismember(T0.("Reaction Name"), list_zero_flux);

% Remove rows from the table
T0(rows_to_remove, :) = [];
writetable(T0, 'CO_to_PHB_model2.csv');

exchange_rxns = findExcRxns(model,1,1);

list_exchange = {};

for i=1:length(model.rxns)
    if exchange_rxns(i)==1
       list_exchange = [list_exchange, model.rxns(i)];
    end
end

rxns_to_remove = [list_exchange,list_zero_flux,'ATPase'];

% Identify rows to remove
rows_to_remove = ismember(T0.("Reaction Name"), rxns_to_remove);

% Remove rows from the table
T0(rows_to_remove, :) = [];
writetable(T0, 'CO_to_PHB2.csv');

model = removeRxns(model,rxns_to_remove);    
order = model.rxnNames; % the desired order

rxnID = findRxnIDs(model,model.rxns);

    for i=1:length(rxnID)    
        fluxes(i)=FBAsolution.x(rxnID(i));
    end

save CO_to_PHB_model2