function [ models ] = make_consensus_models( state_trajectories, XY, iEdges, medI )

    % options: not sure if all are necessary
    options.verbosity = 2;
    options.convergenceTolerance = 1e-4;
    options.reversible = 1;
    options.maximumIterations = 100;
    options.equilibrium = 1;
    options.use_java = 1;
    options.tau = 0.1;
    options.assign_states = 1;
    
    
    ssqdSeg = cell(2,length(iEdges));
    for i = 1:numel(ssqdSeg)
        ssqdSeg{i} = [0 0];
    end
    
    % determine segment edges and indices -> store in cell arrays
    segments = cell(1,length(medI));
    segmInds = cell(1,length(medI));
    for i = 1:length(medI)
        [segments{i}, segmInds{i}] = iSegments(medI{i}, iEdges);
    end
    
    % Compute expected transition quantities and standard deviations.  
    Nij = zeros(2);
    nSeg = zeros(size(ssqdSeg));
  
    for i = 1:numel(state_trajectories)
        p_ti = [double((state_trajectories{i}==1))' double((state_trajectories{i}==2))'];
        for t = 1:size(p_ti,1)-1
            Nij = Nij + p_ti(t,:)' * p_ti(t+1,:);
        end
        if ~isempty(segments{i})
            for j = 1:length(segmInds{i})
                for k = 1:2
                    tmpIND = intersect(find(state_trajectories{i}==k),segments{i}(j,1):segments{i}(j,2));
                    for coord = 1:2
                        ssqdSeg{k,segmInds{i}(j)}(coord) = ssqdSeg{k,segmInds{i}(j)}(coord) + sum((XY{i}(coord,tmpIND) - mean(XY{i}(coord,tmpIND))).^2);
                    end
                    nSeg(k,segmInds{i}(j)) = nSeg(k,segmInds{i}(j)) + numel(tmpIND);
                end
            end 
        end
    end

    % generate output: models
    models = cell(length(iEdges),1);
    
    models{1}.Tij = transition_matrix_mle(Nij, options);
    models{1}.logTij = log(models{1}.Tij);
    models{1}.Pi = stationary_probability(models{1}.Tij);
    models{1}.logPi = log(models{1}.Pi);
    for i = 2:numel(models)
        models{i} = models{1};
    end
    for i = 1:numel(models)
        models{i}.nstates = 2;
        for k = 2:-1:1
            models{i}.states{k,1}.mu = [0 0];
            models{i}.states{k,1}.sigma = sqrt(ssqdSeg{k,i}./(nSeg(k,i)-1));
        end
    end
end

