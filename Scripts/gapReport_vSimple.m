function gapReport_vSimple(model)
% gapReport_vSimple  Perform a basic gap analysis on a metabolic model.
%
%   gapReport_vSimple(model)
%
%   PURPOSE:
%       This function runs two common gap-detection checks on a genome-scale
%       metabolic model (GEM) without requiring model simplification.
%       It is designed for quick use on draft models, even if incomplete.
%
%   CHECKS PERFORMED:
%       1. Dead-end metabolites:
%            - Metabolites that are only produced or only consumed.
%       2. Blocked reactions:
%            - Reactions that cannot carry flux under default constraints
%              (both maximum and minimum achievable flux are ~0).
%
%   INPUT:
%       model   - struct, GEM in COBRA or RAVEN format
%                 Required fields: .mets, .rxns, .S
%                 Optional fields: .c (objective coefficients)
%
%   OUTPUT:
%       None (results are printed to the MATLAB command window)
%
%   NOTES:
%       - Uses RAVEN's detectDeadEnds and solveLP functions.
%       - Blocked reaction detection is based on linear programming flux tests.
%       - No modifications are made to the model.
%

    %-----------------------------
    % 1. Dead-end metabolites
    %-----------------------------
    % Use RAVEN's built-in dead-end detection to find metabolites that
    % cannot participate in steady-state flux (produced-only or consumed-only).
    deadEnds = detectDeadEnds(model);
    fprintf('Dead-end metabolites: %d\n', numel(deadEnds));
    if ~isempty(deadEnds)
        % Print list of metabolite IDs
        disp(model.mets(deadEnds));
    end

    %-----------------------------
    % 2. Blocked reactions
    %-----------------------------
    % A reaction is considered blocked if its maximum and minimum possible
    % flux (given default bounds) are both ~0, or if LP solving fails.
    tol = 1e-7; % tolerance for "zero flux"
    blockedRxns = false(numel(model.rxns), 1);

    % Reset objective function if present
    if isfield(model, 'c')
        model.c(:) = 0;
    end

    % Test each reaction individually
    for i = 1:numel(model.rxns)
        % Create a temporary model with the i-th reaction as objective
        tempModel = model;
        tempModel.c(:) = 0;
        tempModel.c(i) = 1;  % maximize flux through reaction i
        solMax = solveLP(tempModel, 1); % optimization direction: maximize

        % Now minimize the same reaction by maximizing negative flux
        tempModel.c(i) = -1;
        solMin = solveLP(tempModel, 1); % maximize of -flux = minimize flux

        % Mark as blocked if solver failed OR both directions ~0
        if solMax.stat ~= 1 || solMin.stat ~= 1 || ...
           (abs(solMax.f) < tol && abs(solMin.f) < tol)
            blockedRxns(i) = true;
        end
    end

    % Report blocked reactions
    fprintf('Blocked reactions: %d\n', sum(blockedRxns));
    if any(blockedRxns)
        % Print list of blocked reaction IDs
        disp(model.rxns(blockedRxns));
    end
end
