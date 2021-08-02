
%% Load data
label = 'RNAseq'; % 'RNAseq' or 'proteomics'
load('recon2.2s_sarscov2.mat');
load('iAB_AMO1410_SARS-CoV-2.mat');
if strcmp(label, 'RNAseq')
    load('RNA0vs24.mat');
    load('RNA0vs48.mat');
    load('RNA0vs72.mat');
    expr0vs24 = RNA0vs24;
    expr0vs48 = RNA0vs48;
    expr0vs72 = RNA0vs72;
elseif strcmp(label, 'proteomics')
    load('prot0vs24.mat')
    load('prot0vs48.mat')
    load('prot0vs72.mat')
    expr0vs24 = prot0vs24;
    expr0vs48 = prot0vs48;
    expr0vs72 = prot0vs72;
end


%% Add viral biomass pseudo-reaction
% defined in iAB_AMO1410_SARS-CoV-2, except for the lipids which are defined
% in Table 1 at https://doi.org/10.3390/genes12060796
vbof_idx = macrophage_SARS_CoV_2.S(:, strcmp(macrophage_SARS_CoV_2.rxns, 'VBOF')) ~= 0;
vbof_coeffs = macrophage_SARS_CoV_2.S(vbof_idx, strcmp(macrophage_SARS_CoV_2.rxns, 'VBOF'));
vbof_mets = macrophage_SARS_CoV_2.mets(vbof_idx);
vbof_mets = strrep(vbof_mets, '_DASH_', '_'); % homogenise metabolite ids
for i = 1:length(vbof_mets)
    str = vbof_mets{i}(1:end-2);
    vbof_mets{i} = [str '[c][]'];
end
vbof_coeffs = [vbof_coeffs; [-0.03152; -0.02110; -0.00374; -0.00102; -0.02093]]; % add lipids
vbof_mets = [vbof_mets; {'pchol_hs[c][]'; 'pe_hs[c][]'; 'pail_hs[c][]'; 'ps_hs[c][]'; 'chsterol[c][]'}];
model = addReaction(model, 'SARS-CoV-2_biomass_reaction', 'metaboliteList', vbof_mets, 'stoichCoeffList', vbof_coeffs, 'reversible', false, 'printLevel', 0);


%% Compute reaction expression
[reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length] = compute_reaction_expression(model);


%% Set growth conditions
model.lb(strcmp(model.rxns,'EX_o2(e)')) = -1000; % manually set unlimited oxygen uptake
model.lb(strcmp(model.rxns,'EX_o2s(e)')) = 0; % avoid consumption of superoxide. See Table 2 at https://doi.org/10.3389/fphys.2016.00327
model.ub(strcmp(model.rxns,'EX_o2(e)')) = 0; % avoid production of oxygen. See Table 2 at https://doi.org/10.3389/fphys.2016.00327


%% Obtain the flux distributions
changeCobraSolver('ibm_cplex', 'LP');
objectives = model.rxns; % the reaction list for FVA
genes = model.genes; % the list of HGNC ids from the model
minFlux = zeros(length(objectives), 4);
maxFlux = zeros(length(objectives), 4);
vbof_loads = [1/3; 2/3; 1; 0];
for t = 1:4
    disp(t)
    if t == 1
        genes_in_dataset = expr0vs24.hgnc_id; % gets the list of HGNC ids from the covid data
        expr_profile = expr0vs24.FC; % gene expression data
    elseif t == 2
        genes_in_dataset = expr0vs48.hgnc_id; % gets the list of HGNC ids from the covid data
        expr_profile = expr0vs48.FC; % gene expression data
    elseif t == 3
        genes_in_dataset = expr0vs72.hgnc_id; % gets the list of HGNC ids from the covid data
        expr_profile = expr0vs72.FC; % gene expression data
    elseif t == 4 % reference case  
        genes_in_dataset = genes; % gets the list of HGNC ids from the covid data
        expr_profile = ones(numel(genes),1); % gene expression data
    end
    
    x = ones(numel(genes),1);
    for i=1:numel(genes)
        position = find(strcmp(genes_in_dataset,genes{i}));
        if ~isempty(position)
            if length(position) > 1
                x(i) = mean(expr_profile(position));
            else
                x(i) = expr_profile(position);
            end
        end 
    end
    
    [f, F] = evaluate_objective_FVA(x, model, genes, reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length, objectives, vbof_loads(t));
    minFlux(:, t) = f;
    maxFlux(:, t) = F;
end

save(['simulation_results_' label])



%% Functions
%%
function [reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length] = compute_reaction_expression(model)

    genesets = model.grRules;
    genesets = regexprep(genesets,' AND ',' and '); 
    genesets = regexprep(genesets,' OR ',' or '); 

    reaction_expression = cell(length(genesets),1);
    reaction_expression(:) = {''};
    for i = 1:length(genesets)
        str_geneset = genesets{i};
        aux = associate_genes_reactions(str_geneset);
        reaction_expression{i} = aux; 
    end
    reaction_expression=strrep(reaction_expression,' ','');

    genes = model.genes;
    len = NaN(1, length(genes));
    for i = 1:length(genes)
        len(i)=length(genes{i});
    end
    [~, ixs_geni_sorted_by_length] = sort(len,'descend');
    reaction_expression_aux = reaction_expression;
    for i = 1:numel(ixs_geni_sorted_by_length)
        j = ixs_geni_sorted_by_length(i);
        matches = strfind(reaction_expression_aux,genes{j});
        pos_genes_in_react_expr{j} = find(~cellfun('isempty', matches));
        reaction_expression_aux(pos_genes_in_react_expr{j}) = strrep(reaction_expression_aux(pos_genes_in_react_expr{j}),genes{j},'');
    end

end

function str_output = associate_genes_reactions(str)
    % Sometimes grRules have no parentheses, which means that we need
    % to find a smart way to make sure that AND is solved before OR when we
    % substitute MIN and MAX respectively. This means that in the final expression, the MINs have to be
    % calculated before the MAXs.
    % To do so, we substitute the ORs first (which become MAXs), and then the ANDs inside
    % the MAXs. This is to ensure that we have an expression that first solves
    % the  ANDs (which are internal) and then solves the ORs (which are
    % external), thus respecting the common rule that AND is solved before OR
   
    while ( ~isempty(findstr(' or ', str)) ) % loops until all the AND and OR are not found because they have been substituted by MIN and MAX
    
        i = 1;
       while ( (strcmp(str(i:i+3),' or ')==0) )
            i = i+1; % while it does not find any 'and' and any 'or', keeps scrolling the array
        end
        
        str = substitute(str,i);
    end
 
    while ( ~isempty(findstr(' and ', str) ) ) % loops until all the AND and OR are not found because they have been substituted by MIN and MAX
    
        i = 1;
        
        while ( (strcmp(str(i:i+4),' and ')==0)    )
            i = i+1; % while it does not find any 'and' and any 'or', keeps scrolling the array
        end
        
        str = substitute(str,i);
    end
    
    if (isempty(findstr('max(', str)) && isempty(findstr('min(', str)) && ~strcmp(str,''))
        str_output = str;  
    else
        str_output = str; %if str is empty or is a combination of max/min of genes, we leave it as it is
    end
end

function str = substitute(str,i)
    i = i+1;
    % i is now positioned on the initial character of either 'and' or 'or'
    if (str(i)=='a')
        found_and = 1;
    else
        found_and = 0 ;
    end

    bracket_found = 0;
    j = i;
    while (  (strcmp(str(j),'(')==0) || (bracket_found~=-1) ) && (strcmp(str(j),',')==0 || (bracket_found~=0) ) && (j>1)
        j = j-1;
        if (str(j)==')')
            bracket_found = bracket_found+1; % signals further parentheses found along the path
        end
        if (str(j)=='(')
            bracket_found = bracket_found-1; % signals the closure of parentheses found along the path
        end
    end
    if (bracket_found == -1 || strcmp(str(j),',')~=0)
        j = j+1;
    end
    if (found_and == 1)
        k = i+3;
    else
        k = i+2;
    end

    bracket_found = 0;
    while ( (strcmp(str(k),')')==0) || ( bracket_found~=-1 )) && (strcmp(str(k),',')==0 || (bracket_found~=0) ) && (k<length(str))
        k = k+1;
        if (str(k)=='(')
            bracket_found = bracket_found+1; % signals further parentheses found along the path
        end
        if (str(k)==')')
            bracket_found = bracket_found-1; % signals the closure of parentheses found along the path
        end
    end
    if bracket_found == -1 ||  strcmp(str(k),',')~=0
        k = k-1;
    end

    if (found_and == 1)
        str_new = strrep( str, str(j:k), [' min(',str(j:i-1),', ',str(i+4:k), ') '] );
    else
        str_new = strrep( str, str(j:k), [' max(',str(j:i-1),', ',str(i+3:k), ') '] );
    end
    str = str_new;
end

function [minFlux, maxFlux] = evaluate_objective_FVA(x, model, genes, reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length, objectives, vbof_load)

yt=x';      % x' is the transpose of x, that is the gene expression array

eval_reaction_expression = reaction_expression;

for i=ixs_geni_sorted_by_length %loop over the array of the non-1 gene expressions, in order to replace the names of genes in geni_reazioni.mat with their values. All the gene set expressions with only 1s as gene values , will be left empty and at the end of this loop everything empty will be substituted with 1 anyway. This avoids looping over all the genes yt, which is very expensive
    posizioni_gene = pos_genes_in_react_expr{i};
    for j=1:length(posizioni_gene) %for each string of reaction_expression, we replace the substring 'bXXXX' with the number representing its gene expression
        eval_reaction_expression{posizioni_gene(j)} = strrep(eval_reaction_expression{posizioni_gene(j)}, genes{i}, num2str(yt(i),'%.15f'));  %Matlab strangely truncates decimal digits when using num2str. Addimg %.12f at least ensures that 12 decimal digits are included in the number converted into string
    end
end
eval_reaction_expression( cellfun(@isempty, eval_reaction_expression) ) = {'1.0'};  %replaces all the empty cells of gene expression (e.g. exchange reactions) with 1, i.e. gene expressed nomally

num_reaction_expression = zeros(1,length(eval_reaction_expression));
for i=1:length(num_reaction_expression)
    str = eval_reaction_expression{i};
    
    num_parenthesis = numel(strfind(str,')'));
    while (num_parenthesis > 32) %if there are more than 32 parentheses, matlab is unable to run EVAL. So we need to reduce these parentheses manually by starting to eval smaller pieces of the string
        to_replace = 'min.\d*+\.+\d*,\d*+\.+\d*.|max.\d*+\.+\d*,\d*+\.+\d*.|min..\d*+\.+\d*.,\d*+\.+\d*.|max..\d*+\.+\d*.,\d*+\.+\d*.|min..\d*+\.+\d*.,.\d*+\.+\d*..|max..\d*+\.+\d*.,.\d*+\.+\d*..|min.\d*+\.+\d*,.\d*+\.+\d*..|max.\d*+\.+\d*,.\d*+\.+\d*..';  %searches for all the strings of kind min(NUM.NUM,NUM.NUM) or max(NUM.NUM,NUM.NUM) or  min((NUM.NUM),NUM.NUM) or max((NUM.NUM),NUM.NUM) or  min((NUM.NUM),(NUM.NUM)) or max(NUM.NUM,(NUM.NUM)) or  min(NUM.NUM,(NUM.NUM)) or max((NUM.NUM),(NUM.NUM))
        substrings_to_replace = regexp(str, to_replace, 'match');
        if isempty(substrings_to_replace)
            num_parenthesis = 0; %if num_parenthesis > 32 and there is nothing that can be replaced with regexp, we force this, in order to avoid an endless loop. Later, eval will catch an exception as it cannot evaluate when num_parenthesis>32
        else
            for j = 1:numel(substrings_to_replace)
                ss_rep = substrings_to_replace{j};
                str = strrep(str,ss_rep,num2str(eval(ss_rep),'%.15f'));
            end
            num_parenthesis = numel(strfind(str,')'));
       end
    end
    
    str = regexprep(str,'/','');
    
    num_reaction_expression(i) = eval(str);   %evaluates the cells like they are numerical expressions (so as to compute min and max of gene expressions)
    
end

gamma = ones(1,length(reaction_expression));

for i=1:length(num_reaction_expression)   %loop over the array of the geneset expressions
        model.lb(i) = model.lb(i)*(num_reaction_expression(i)^gamma(i)); % model default x gene expression to the power of gamma
        model.ub(i) = model.ub(i)*(num_reaction_expression(i)^gamma(i));
end

% focus solution space
model = changeObjective(model, 'biomass_reaction');
FBAsolution = optimizeCbModel(model);
max_biomass = FBAsolution.f;
model = changeObjective(model, 'GNDc');
FBAsolution = optimizeCbModel(model);
max_GNDc = FBAsolution.f;
model = changeObjective(model, 'ACOAO7p');
FBAsolution = optimizeCbModel(model);
max_ACOAO7p = FBAsolution.f;
model = changeObjective(model, 'r0173');
FBAsolution = optimizeCbModel(model);
max_r0173 = FBAsolution.f;
model = changeRxnBounds(model, 'biomass_reaction', 0.5*max_biomass, 'l');
model = changeRxnBounds(model, 'GNDc', 0.8*max_GNDc, 'l');
model = changeRxnBounds(model, 'ACOAO7p', 0.8*max_ACOAO7p, 'l');
model = changeRxnBounds(model, 'r0173', 0.8*max_r0173, 'l');
model = changeRxnBounds(model, 'SARS-CoV-2_biomass_reaction', vbof_load*(0.5/3)*max_biomass, 'l');

% FVA
model.c = zeros(length(model.c), 1);
[minFlux, maxFlux] = fluxVariability(model, 100, 'max', objectives, 0);

end