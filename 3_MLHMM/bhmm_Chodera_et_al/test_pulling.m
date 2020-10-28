% Test Bayesian HMM with real pulling data.

clear;

dataDirectory = '../datasets/elms-high-resolution/';
stateIndex = 47;

% Read bead positions.
AXPos = textread(sprintf('%s/%d.AXPos', dataDirectory, stateIndex));
BXPos = textread(sprintf('%s/%d.BXPos', dataDirectory, stateIndex));

% Read trap separation.
trapSep = textread(sprintf('%s/%d.trapSep', dataDirectory, stateIndex));

% Compute bead-to-bead extension.
o_t = trapSep + BXPos - AXPos;
T = length(o_t);

% Pack data into cellarray.
data{1} = o_t';

% Set number of states.
nstates = 4;

% Determine maximum-likelihood HMM.
%disp('Determining maximum-likelihood model...');
%mlmodel = mlhmm(data, model.nstates);

options = bhmm_default_options();
options.maximumIterations = 10; % maximum number of allowed iterations
options.convergenceTolerance = 1.0e-4; % relative convergence tolerance in likelihood

% Options for transition matrix estimation.
options.reversible = true; % infer reversible transition matrices
options.diagonally_dominant = false; % don't enforce diagonally-dominant transition matrices
options.verbosity = 2; % set verbosity level
options.equilibrium = false; % trajectory data is not initially drawn from equilibrium

options.tau = 4e-4; % seconds between observations

% Generate an initial (poor) model from data.
%initial_model = generate_initial_model(data, nstates, options);

% Determine maximum-likelihood HMM.
disp('Determining maximum-likelihood model...');
mlmodel = mlhmm(data, nstates, options);

% Initialize Bayesian HMM.
disp('Sampling from Bayesian posterior...');
nsamples = 1000;
models = bhmm(data, mlmodel, nsamples, options);

% Plot
nburnin = 50; % number of samples to discard to equilibration
plot_state_assignments(data, models(nburnin+1:end), options);
options.sameaxis = false;
analyze_states(models(nburnin+1:end, options));
