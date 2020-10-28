function [ mlmodels, state_trajectory, arxv ] = mlhmmINTsegments( XY, iTrace, iThreshs, init, varargin )
%mlhmmINTsegments: Maximum-likelihood Hidden Markov Model Analyis of a
%single particle XY-time-trace, divided into similar-intensity segments. 

%% parse input
p = inputParser;
addRequired(p,'XY');
addRequired(p, 'iTrace');
addRequired(p, 'iThreshs');
addRequired(p, 'init');
%addOptional(p, 'xyLims', [Inf Inf]);
addParameter(p, 'options', []);
addParameter(p, 'sigmas',[]);
addParameter(p, 'initial_models', []);
parse(p, XY, iTrace, iThreshs, init, varargin{:})

XY = p.Results.XY;
iTrace = p.Results.iTrace;
iThreshs = p.Results.iThreshs;
[segments, segmInds] = iSegments(iTrace, iThreshs);

if ~isempty(p.Results.sigmas)
    sigmas = p.Results.sigmas;
    modelmode = 'sigma';
elseif ~isempty(p.Results.initial_models)
    models = p.Results.initial_models(segmInds);
    modelmode = 'model';
end

if ~isempty(p.Results.options)
    mlhmmOptions = p.Results.options;
else
    mlhmmOptions.verbosity = 2;
    mlhmmOptions.convergenceTolerance = 1e-4;
    mlhmmOptions.reversible = 1;
    mlhmmOptions.maximumIterations = 100;
    mlhmmOptions.equilibrium = 1;
    mlhmmOptions.use_java = 1;
    mlhmmOptions.tau = 0.1;
    mlhmmOptions.assign_states = 1;
end

if ~isempty(segments)
    XY = XY(:,1:segments(end,2));
    switch modelmode
        case 'sigma'
        [mlmodels, state_trajectory] = mlhmmXYsegments(XY, segments, 2, mlhmmOptions, sigmas);
        case 'model'
        [mlmodels, state_trajectory] = mlhmmXYsegments(XY, segments, 2, mlhmmOptions, [], models);
    end
    arxv.XY = XY;
else
    mlmodels = [];
    state_trajectory = [];
    arxv.XY = [];
end

arxv.segments = segments + init - 1;
arxv.models = cell(size(mlmodels));
for i = 1:length(mlmodels)
    arxv.models{i}.Tij = mlmodels{i}.Tij;
    arxv.models{i}.Pi = mlmodels{i}.Pi;
    arxv.models{i}.states = cell(2,1);
    for j = 1:2
        arxv.models{i}.states{j}.mu = mlmodels{i}.states{j}.mu;
        arxv.models{i}.states{j}.sigma = mlmodels{i}.states{j}.sigma;
    end
end

