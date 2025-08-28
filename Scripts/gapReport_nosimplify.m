function [noFluxRxns, noFluxRxnsRelaxed, subGraphs, notProducedMets, minToConnect,...
    neededForProductionMat, canProduceWithoutInput, canConsumeWithoutOutput, ...
    connectedFromTemplates, addedFromTemplates]=gapReport_nosimplify(model, templateModels)
% gapReport_nosimplify
%   Performs a gap analysis and summarizes the results 
%
%   model                       a model structure
%   templateModels              a cell array of template models to use for
%                               gap filling (optional)
%
%   noFluxRxns                  cell array with reactions that cannot carry
%                               flux
%   noFluxRxnsRelaxed           cell array with reactions that cannot carry
%                               flux even if the mass balance constraint is 
%                               relaxed so that it is allowed to have 
%                               net production of all metabolites
%   subGraphs                   structure with the metabolites in each of
%                               the isolated sub networks
%   notProducedMets             cell array with the metabolites that
%                               couldn't have net production
%   minToConnect                structure with the minimal number of
%                               metabolites that need to be connected in 
%                               order to be able to produce all other 
%                               metabolites and which metabolites each of
%                               them connects
%   neededForProductionMat      matrix where n x m is true if metabolite n
%                               allows for production of metabolite m
%   canProduceWithoutInput      cell array with metabolites that could be
%                               produced even when there is no input to the
%                               model
%   canConsumeWithoutOutput     cell array with metabolites that could be
%                               consumed even when there is no output from
%                               the model
%   connectedFromTemplates      cell array with the reactions that could be
%                               connected using the template models
%   addedFromTemplates          structure with the reactions that were
%                               added from the template models and which 
%                               model they were added from
%
% Usage: [noFluxRxns, noFluxRxnsRelaxed, subGraphs, notProducedMets, ...
%    minToConnect, neededForProductionMat, connectedFromTemplates, ...
%    addedFromTemplates] = gapReport_nosimplify(model, templateModels)

if nargin<2
    templateModels=[];         % if no template models are provided
    connectedFromTemplates=[]; % initialize empty
    addedFromTemplates=[];     % initialize empty
end

% Print header with model info
fprintf(['Gap analysis for ' model.id ' - ' model.name '\n\n']);

% If model has an "unconstrained" field, simplify it for proper flux checks
if isfield(model,'unconstrained')
    calculateINOUT=true;       % flag: check for producible/consumable mets w/o exchange
    closedModel=model;         % keep original for input/output checks
    model=simplifyModel(model);% remove unconstrained exchange reactions
else
    canConsumeWithoutOutput={}; % initialize empty if no unconstrained field
    canProduceWithoutInput={};
    calculateINOUT=false;
end

% Create a relaxed version of the model where metabolites can accumulate
model2=model;
model2.b=[model2.b inf(numel(model2.mets),1)];

% Identify reactions that can carry flux in strict and relaxed cases
I=haveFlux_nosimplify(model);       % strict case
noFluxRxns=model.rxns(~I);          % reactions with no flux under strict balance
J=haveFlux_nosimplify(model2);      % relaxed case
noFluxRxnsRelaxed=model2.rxns(~J);  % reactions with no flux under relaxed balance

% Remove blocked reactions to get reduced models
bModel=removeReactions(model,~I,true,true);   % strict case
cModel=removeReactions(model2,~J,true,true);  % relaxed case

% Print overview of reaction/metabolite coverage
fprintf('***Overview\n');
fprintf([num2str(numel(model.rxns)-sum(I)) ' out of ' ...
    num2str(numel(model.rxns)) ' reactions cannot carry flux (' ...
    num2str(numel(model.rxns)-sum(J)) ' if net production of all ' ...
    'metabolites is allowed)\n']);
fprintf([num2str(numel(model.mets)-numel(bModel.mets)) ' out of ' ...
    num2str(numel(model.mets)) ' metabolites are unreachable (' ...
    num2str(numel(model.mets)-numel(cModel.mets)) ' if net production of ' ...
    'all metabolites is allowed)\n']);

% Identify isolated subnetworks of metabolites
fprintf('\n***Isolated subnetworks\n');
subGraphs=getAllSubGraphs(model);
fprintf(['A total of ' num2str(size(subGraphs,2)) ...
    ' isolated sub-networks are present in the model\n']);
for i=1:size(subGraphs,2)
    fprintf(['\t' num2str(i) '. ' num2str(sum(subGraphs(:,i))) ...
        ' metabolites\n']);
end

% Check which metabolites can (or cannot) be produced and connectivity reqs
fprintf('\n***Metabolite connectivity\n');
[notProducedMets, ~, neededForProductionMat, minToConnect] = ...
    checkProduction_nosimplify(model,true,model.comps,false);
fprintf(['To enable net production of all metabolites, a total of ' ...
    num2str(numel(minToConnect)) ' metabolites must be connected\n']);
fprintf('Top 10 metabolites to connect:\n');
for i=1:min(10,numel(minToConnect))
    fprintf(['\t' num2str(i) '. ' minToConnect{i} '\n']);
end

% If the model had unconstrained reactions, check for leaks/sinks
if calculateINOUT==true
    fprintf('\n***Mass balancing\n');
    produced=canProduce(closedModel);         % mets produced without input
    canProduceWithoutInput=closedModel.mets(produced);
    consumed=canConsume(closedModel);         % mets consumed without output
    canConsumeWithoutOutput=closedModel.mets(consumed);
    fprintf([num2str(numel(canConsumeWithoutOutput)) ...
        ' metabolites could be consumed without any outputs\n' ...
        num2str(numel(canProduceWithoutInput)) ...
        ' metabolites could be produced without any inputs\n']);
end

% If template models are provided, attempt automated gap-filling
if ~isempty(templateModels)
    fprintf('\n***Automated gap-filling\n');
    [connectedFromTemplates, ~, addedFromTemplates] = ...
        fillGaps(model,templateModels);
    
    % List all template model IDs
    t=templateModels{1}.id;
    for i=2:numel(templateModels)
        t=[t ', ' templateModels{i}.id];
    end
    
    % Print how many reactions could be reconnected by gap-filling
    fprintf([num2str(numel(connectedFromTemplates)) ...
        ' unconnected reactions can be connected by including ' ...
        num2str(numel(addedFromTemplates)) ' reactions from\n' t '\n']);
end
end
