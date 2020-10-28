function show_models(models)
% Print the parameters of the specified model.
%
% show_model(model)
%
% ARGUMENTS
%  model (struct) - model to be shown

% Display transition matrix.
disp('Tij = ');
disp(models{1}.Tij);

% Display states.
for segmentIndex = 1:length(models)
    for stateIndex = 1:models{segmentIndex}.nstates
      state = models{segmentIndex}.states{stateIndex};
      fprintf('Segment %4d, state %4d : Pi = %8.6f  mu = %12.7f %12.7f  sigma = %12.7f %12.7f\n', ...
          segmentIndex, stateIndex, models{segmentIndex}.Pi(stateIndex), state.mu(1), state.mu(2), state.sigma(1), state.sigma(2));
    end
end

return
