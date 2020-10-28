% Bayesian hidden Markov model (BHMM) analysis of timeseries data
%
% The data is assumed to come from a fixed number of states each emitting real-valued signals corresponding to 
% one-dimensional i.i.d. normal random variates.  Each state is characterized by a mean and variance, inferred 
% from the data.  Transitions between the states are governed by a transition matrix, and dynamics is assumed
% to be Markovian.
%
%
% EXAMPLE USE
%
% % Find maximum-likelihood model for data.
% maximum_likelihood_model = mlhmm(data, nstates);
%
% % Generate samples from Bayesian HMM.
% bayesian_model_samples = bhmm(data, maximum_likelihod_model);
%
%
% PROVIDED FUNCTIONS
%
% * Test harnesses
%
% test_bhmm                     - test Bayesian HMM on synthetic data
% test_pulling                  - test Bayesian HMM on real pulling data
%
% * Hidden Markov model construction and sampling:
%
% bhmm_default_options          - generate an options structure for BHMM containing default parameters
% mlhmm                         - determine maximum-likelihood HMM using EM algorithm
% bhmm                          - generate samples from Bayesian HMM
%
% * Analysis of BHMM output:
%
% analyze_rates                 - analyze transition rate distributions
% analyze_states                - analyze state emission probabilities (NOT IMPLEMENTED YET)
%
% * Auxiliary functions:
%
% em_gaussian_mixture           - fit a mixture of Gaussians using EM algorithms (used for initialization)
%
% plot_state_assignments        - plot colored state assignments from a particular model
%
% update_emission_probabilities - update state emission parameters by sampling E ~ P(E | S,O)
% update_transition_matrix      - update transition matrix by sampling T ~ P(T | S)
% update_state_trajectories     - update state tranjectories by sampling S ~ P(S | T,E,O)
% rate_matrix_sample            - generate an uncorrelated sample from the rate matrix
% transition_matrix_mle         - generate maximum-likelihood estimate of transition matrix
% transition_matrix_sample      - generate a new (potentially correlated) sample from the transition matrix
%
% viterbi                       - implementation of the Viterbi algorithm for computing most likely state sequences
% forward_backward_algorithm    - implementation of the forward-backward algorithm for computing probabilities of state sequences
% baum_welch                    - implementation of the Baum-Welch algorithm (UNTESTED)
% sample_state_trajectory       - sample state sequence from Bayesian posterior using Nina's algorithm
% 
% pack_vector, unpack_vector    - methods used by rate matrix sampling
% argmax                        - index of maximum argument in a vector
% count_transitions             - count number of observed transitions between states given a state sequence
% draw                          - draw one of a set of possible outcomes
% exportfig                     - render nice figure (from Mathworks)
% logsum                        - compute logarithm of sum of exponentials
% 
% generate_test_model           - generate a test model and synthetic data
% show_model                    - print the parameters of the current model
% show_trajectory               - print a state sequence trajectory
%
% compute_log_likelihood        - compute log likelihood of an HMM given the data



