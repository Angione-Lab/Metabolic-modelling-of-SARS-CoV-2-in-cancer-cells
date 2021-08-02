
labels = {'RNAseq', 'proteomics'};
fluxFoldChange_pathways_all = [];

% repeat for results obtained from both datasets
for i = 1:2
    data_label = labels{i};
    load(['simulation_results_' data_label '.mat']);

    %% Fix pathways
    model.subSystems(strcmp(model.subSystems, '')) = {'Unassigned'};
    model.subSystems(strcmp(model.subSystems, 'beta-Alanine metabolism')) = {'Beta-Alanine metabolism'};
    model.subSystems(7786:end) = {'Protein translation and secretion'};
    proteins = {'P60568','P01589','P14784','P31785','P05231','P08887','P40189','P13232','P16871','P01579','P01584','P18510',...
                'P05112','P05113','P10145','P15248','P22301','P29459','P29460','P35225','P40933','Q96F46','Q9Y2Y9','P09919',...
                'P02778','P13500','P10147','P01375','P15692','P01127','P13236','P04141','P09038','P51671','P78556','Q969D9',...
               'Q06830','P32119','Q16552','Q96PD4','Q8NAC3','P01009','P02763','P00738','P01011','P02765','P02647','P02787',...
               'P02787','P02741'};
    for i = 1:length(proteins)
        a = contains(model.rxns, proteins{i});
        model.subSystems(a) = {[proteins{i} ' translation and secretion']};
    end


    %% Calculate reaction activity fold changes
    % set to 0 values below cplex feasibility tolerance
    maxFlux(maxFlux < 1e-09) = 0;
    minFlux(minFlux < 1e-09) = 0;
    % split reversible reactions
    revIDs = model.lb < 0;
    rxns = [model.rxns; strcat(model.rxns(revIDs), '_b')];
    rxnNames = [model.rxnNames; strcat(model.rxnNames(revIDs), ' (backward)')];
    subSystems = [model.subSystems; model.subSystems(revIDs)];
    formulas = printRxnFormula(model, model.rxns, false);
    formulas = strrep(formulas, '[]', '');
    formulas = [formulas; formulas(revIDs)];
    maxFlux = [maxFlux; -minFlux(revIDs, :)];
    maxFlux(maxFlux < 0) = 0;

    % re-sort reactions
    [rxns, sort_idx] = sort(rxns);
    rxnNames = rxnNames(sort_idx);
    subSystems = subSystems(sort_idx);
    formulas = formulas(sort_idx);
    maxFlux = maxFlux(sort_idx, :);

    % % get solution span for each flux
    fluxFoldChange = maxFlux(:, 1:3) ./ repmat(maxFlux(:, 4), 1, 3);
    % set to 1 the NaNs due to 0/0 divisions
    fluxFoldChange(isnan(fluxFoldChange)) = 1;


    %% Export reaction data
    t1 = [{'Reaction ID', 'Reaction name', 'Formula', 'Pathway', 'Max flux 24h', 'Max flux 48h', 'Max flux 72h', 'Max flux control', 'FC 24h', 'FC 48h', 'FC 72h'}; ...
        rxns, rxnNames, formulas, subSystems, num2cell(maxFlux), num2cell(fluxFoldChange)];
    xlswrite(['covid19_FVA_results_' data_label], t1, 'reactions');

    % sort reaction fold changes
    t3 = t1(2:end, [1,2,3,4,9]);
    [~, a] = sort(cell2mat(t3(:, end)), 'descend');
    t3 = [{'Reaction ID', 'Reaction name', 'Formula', 'Pathway', 'FC 24h'}; t3(a, :)];
    xlswrite(['covid19_FVA_results_' data_label], t3, 'reactions24h');
    t3 = t1(2:end, [1,2,3,4,10]);
    [~, a] = sort(cell2mat(t3(:, end)), 'descend');
    t3 = [{'Reaction ID', 'Reaction name', 'Formula', 'Pathway', 'FC 48h'}; t3(a, :)];
    xlswrite(['covid19_FVA_results_' data_label], t3, 'reactions48h');
    t3 = t1(2:end, [1,2,3,4,11]);
    [~, a] = sort(cell2mat(t3(:, end)), 'descend');
    t3 = [{'Reaction ID', 'Reaction name', 'Formula', 'Pathway', 'FC 72h'}; t3(a, :)];
    xlswrite(['covid19_FVA_results_' data_label], t3, 'reactions72h');


    %% Pathway statistics
    [fluxFoldChange_pathways, fluxFoldChange_pathways_sem, fluxFoldChange_pathways_std] = grpstats(fluxFoldChange, subSystems, {'mean', 'sem', 'std'});
    uniqueSubSystems = unique(subSystems, 'stable');
    [uniqueSubSystems, sort_idx] = sort(uniqueSubSystems);
    fluxFoldChange_pathways = fluxFoldChange_pathways(sort_idx, :);
    fluxFoldChange_pathways_sem = fluxFoldChange_pathways_sem(sort_idx, :);
    fluxFoldChange_pathways_std = fluxFoldChange_pathways_std(sort_idx, :);
    t2 = [{'Pathway', 'FC mean 24h', 'FC mean 48h', 'FC mean 72h', 'FC sem 24h', 'FC sem 48h', 'FC sem 72h', 'FC std 24h', 'FC std 48h', 'FC std 72h'}; ...
        uniqueSubSystems, num2cell(fluxFoldChange_pathways), num2cell(fluxFoldChange_pathways_sem), num2cell(fluxFoldChange_pathways_std)];
    xlswrite(['covid19_FVA_results_' data_label], t2, 'pathways');

    % sort pathway fold changes
    t3 = t2(2:end, [1,2,5,8]);
    [~, a] = sort(cell2mat(t3(:, 2)), 'descend');
    t3 = [{'Pathway', 'FC mean 24h', 'FC sem 24h', 'FC std 24h'}; t3(a, :)];
    xlswrite(['covid19_FVA_results_' data_label], t3, 'pathways24h');
    t3 = t2(2:end, [1,3,6,9]);
    [~, a] = sort(cell2mat(t3(:, 2)), 'descend');
    t3 = [{'Pathway', 'FC mean 48h', 'FC sem 48h', 'FC std 48h'}; t3(a, :)];
    xlswrite(['covid19_FVA_results_' data_label], t3, 'pathways48h');
    t3 = t2(2:end, [1,4,7,10]);
    [~, a] = sort(cell2mat(t3(:, 2)), 'descend');
    t3 = [{'Pathway', 'FC mean 72h', 'FC sem 72h', 'FC std 72h'}; t3(a, :)];
    xlswrite(['covid19_FVA_results_' data_label], t3, 'pathways72h');


    %% Pathway enrichment
    a = and(fluxFoldChange > prctile(fluxFoldChange, 95), fluxFoldChange >= 1.5);
    enrichment_table_24h = FEA(a(:, 1), subSystems);
    enrichment_table_48h = FEA(a(:, 2), subSystems);
    enrichment_table_72h = FEA(a(:, 3), subSystems);
    a = and(fluxFoldChange < prctile(fluxFoldChange, 5), fluxFoldChange <= 0.8);
    enrichment_table_24h = [enrichment_table_24h; FEA(a(:, 1), subSystems)];
    enrichment_table_48h = [enrichment_table_48h; FEA(a(:, 2), subSystems)];
    enrichment_table_72h = [enrichment_table_72h; FEA(a(:, 3), subSystems)];
    enrichment_table_24h = [repmat({''}, size(enrichment_table_24h, 1), 1), enrichment_table_24h];
    enrichment_table_24h(strcmp(enrichment_table_24h(:, 2), 'p-value'), 1) = {'UP (FC > 95th percentile and >= 1.5)'; 'DOWN (FC < 5th percentile and <= 0.8)'};
    enrichment_table_48h = [repmat({''}, size(enrichment_table_48h, 1), 1), enrichment_table_48h];
    enrichment_table_48h(strcmp(enrichment_table_48h(:, 2), 'p-value'), 1) = {'UP (FC > 95th percentile and >= 1.5)'; 'DOWN (FC < 5th percentile and <= 0.8)'};
    enrichment_table_72h = [repmat({''}, size(enrichment_table_72h, 1), 1), enrichment_table_72h];
    enrichment_table_72h(strcmp(enrichment_table_72h(:, 2), 'p-value'), 1) = {'UP (FC > 95th percentile and >= 1.5)'; 'DOWN (FC < 5th percentile and <= 0.8)'};
    xlswrite(['covid19_FVA_results_' data_label], enrichment_table_24h, 'enrichment24h');
    xlswrite(['covid19_FVA_results_' data_label], enrichment_table_48h, 'enrichment48h');
    xlswrite(['covid19_FVA_results_' data_label], enrichment_table_72h, 'enrichment72h');

    % metabolite info to interpret the formulas
    model.mets = strrep(model.mets, '[]', '');
    t4 = [{'Metabolite ID', 'Metabolite name'}; ...
        model.mets, model.metNames];
    xlswrite(['covid19_FVA_results_' data_label], t4, 'metabolites');
    
    % aggregate pathway-level mean fold changes
    fluxFoldChange_pathways_all = [fluxFoldChange_pathways_all fluxFoldChange_pathways];
end


%% Plot pathway fold changes
L = 10; % number of datapoints
indexValue = 1; % value for which to set a particular color
topColor = [0.5 0 0]; % color for maximum data value (dark red = [0.5 0 0])
indexColor1 = [1 0 0]; % color for intermediate data value (red = [1 0 0])
indexColor = [1 1 1]; % color for null value (white = [1 1 1])
bottomcolor = [0 0 1]; % color for minimum data value (blue = [0 0 1])
% Calculate where proportionally indexValue lies between minimum and maximum values
largest = max(max(fluxFoldChange_pathways_all));
smallest = min(min(fluxFoldChange_pathways_all));
index = L*abs(indexValue-smallest)/(largest-smallest);
% Create color map ranging from bottom color to index color
% Multiplying number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
            linspace(bottomcolor(2),indexColor(2),100*index)',...
            linspace(bottomcolor(3),indexColor(3),100*index)'];
customCMap3 = [linspace(indexColor(1),indexColor1(1),100*(index))',...
            linspace(indexColor(2),indexColor1(2),100*(index))',...
            linspace(indexColor(3),indexColor1(3),100*(index))'];
% Create color map ranging from index color to top color
% Multiplying number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor1(1),topColor(1),100*(L-2*index))',...
            linspace(indexColor1(2),topColor(2),100*(L-2*index))',...
            linspace(indexColor1(3),topColor(3),100*(L-2*index))'];
customCMap = [customCMap1; customCMap3; customCMap2];  % Combine colormaps

% split pathways for the two heatmaps
idx = contains(uniqueSubSystems, 'translation and secretion', 'IgnoreCase', true);
uniqueSubSystems1 = uniqueSubSystems(~idx);
fluxFoldChange_pathways1 = fluxFoldChange_pathways_all(~idx, :);
uniqueSubSystems2 = uniqueSubSystems(idx);
fluxFoldChange_pathways2 = fluxFoldChange_pathways_all(idx, :);
% fix entries
uniqueSubSystems2(strcmp(uniqueSubSystems2, 'Protein translation and secretion')) = {'Common'};
for i = 1:length(uniqueSubSystems2)
    uniqueSubSystems2{i} = strrep(uniqueSubSystems2{i}, ' translation and secretion', '');
end
% re-sort entries
[uniqueSubSystems2, sort_idx] = sort(uniqueSubSystems2);
fluxFoldChange_pathways2 = fluxFoldChange_pathways2(sort_idx, :);

subplot('Position', [0.45 0.01 0.15 0.99])
h1 = heatmap({'24hpi','48hpi','72hpi','p24hpi','p48hpi','p72hpi'}, uniqueSubSystems1, fluxFoldChange_pathways1, 'Colormap', customCMap, ...
    'ColorbarVisible', 'off', 'ColorLimits', [smallest largest]);
set(gca, 'FontSize', 5, 'FontName', 'Arial')
subplot('Position', [0.75 0.35 0.15 0.65])
h2 = heatmap({'24hpi','48hpi','72hpi','p24hpi','p48hpi','p72hpi'}, uniqueSubSystems2, fluxFoldChange_pathways2, 'Colormap', customCMap, ...
    'ColorLimits', [smallest largest]);
set(gca, 'FontSize', 5, 'FontName', 'Arial')
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperOrientation = 'portrait';
fig.PaperSize = [10, 17];
print('-fillpage', 'time_heatmap', '-dpdf')


%%
%%
function resultCell = FEA(rxnSet, group)
% Significance analysis - Flux enrichment analysis using hypergeometric
% 1-sided test and FDR correction for multiple testing
%
% USAGE:
%
%    resultCell = FEA(model, rxnSet, 'subSystems')
%
% INPUTS:
%    model:           COBRA structure model
%    rxnSet:          reaction set to be enriched (vector of reaction indices e.g. 1:10)
%    group:           model.group structure e.g.
%                    'subSystems' : FEA looks for significantly enriched subsystems in rxnSet
%
% OUTPUT:
%    resultCell:    cell structure of enriched groups
%
% EXAMPLE:
%
%    load ecoli_core_model;
%    resultCell = FEA(modelEcore, 1:10, 'subSystems');
%
% .. Author: Marouen BEN GUEBILA 04/2016
% Modified by Guido Zampieri 11/2020

if nargin < 2
    error('The function FEA must be called with reaction set and group as arguments')
end
if ~isvector(rxnSet)
    error('Please provide the indices of the reactions e.g. 1:10')
end
if ~iscell(group)
    error('Please provide the group name as cell array of characters e.g. the subSystem field of any metabolic model ')
end

%Temporary Warning until FEA statistics are checked for multiple classes.
if iscell(group{1}) %Potentially multiple subSystems
    if any(cellfun(@numel, group) > 2)
        warning('Multiple subSystems detected for some reactions. FEA statistics might not be correct.\n Please consider using only one subSystem per reaction.')
    end
end

% compute frequency of enriched terms
%groups = eval(['model.' group]);
groups = group;
if iscell([groups{:}])
   [uniquehSubsystemsA] = unique([groups{:}]);
   presenceindicator = false(numel(uniquehSubsystemsA),numel(group));
   for i = 1:numel(groups)
       presenceindicator(:,i) = ismember(uniquehSubsystemsA,groups{i});
   end   
   [K,~] = find(presenceindicator);
else
    %This works only for fields which have a single entry.
    [uniquehSubsystemsA, ~, K] = unique(groups);
end
% fetch group
%enRxns = eval(['model.' group '(rxnSet)']);
enRxns = group(rxnSet);
m = length(uniquehSubsystemsA);
allSubsystems = zeros(1, m);

% look for unique occurences
if iscell([enRxns{:}])
   [uniquehSubsystems] = unique([enRxns{:}]);
   presenceindicator = false(numel(uniquehSubsystems),numel(group));
   for i = 1:numel(enRxns)       
        presenceindicator(:,i) = ismember(uniquehSubsystems,enRxns{i});
   end   
   [J,~] = find(presenceindicator);
else
    %This works only for fields which have a single entry.
    [uniquehSubsystems, ~, J] = unique(enRxns);
end

occ = histc(J, 1:numel(uniquehSubsystems));
[l, p] = intersect(uniquehSubsystemsA, uniquehSubsystems);
allSubsystems(p) = occ;

% compute total number of reactions per group
nRxns = histc(K, 1:numel(uniquehSubsystemsA));  % the number of reactions per susbsystem

% Compute p-values
% gopvalues = hygepdf(allSubsystems', max(nRxns), max(allSubsystems), nRxns);
gopvalues = hygecdf(allSubsystems'-1, repmat(sum(nRxns), length(nRxns), 1), nRxns, repmat(sum(allSubsystems), length(nRxns), 1), 'upper');

% take out the zeros for one-sided test
nonZerInd = find(allSubsystems);

% sort p-values
[m, rxnInd] = sort(gopvalues);

% intersect non zero sets with ordered pvalues
[~, nonZeroInd] = intersect(rxnInd, nonZerInd);
orderedPval = rxnInd(sort(nonZeroInd));

% Build result cell
% initilize variable
resultCell = cell(length(orderedPval) + 1, 5);
resultCell(1, :) = {'p-value', 'Adjusted p-value', 'Pathway', 'Enriched set size', 'Total set size'};

% P values
resultCell(2:end, 1) = num2cell(gopvalues(orderedPval));

% correct for multiple testing with FDR
resultCell(2:end, 2) = num2cell(mafdr(cell2mat(resultCell(2:end, 1)), 'BHFDR', true));

% Group name
resultCell(2:end, 3) = uniquehSubsystemsA(orderedPval);

% Test size
resultCell(2:end, 4) = num2cell(allSubsystems(orderedPval))';

% Total group size
resultCell(2:end, 5) = num2cell(nRxns(orderedPval));

end
