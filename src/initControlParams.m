function con = initControlParams(varargin)
%initControlParams Initialize control parameters of GSH with optional overrides
%
% Usage:
%   con = initControlParams(); % uses all default values
%   con = initControlParams('epsilon', 0.5, 'csq', 0.2); % overrides defaults
%                                                           for those inputs
%
% Output:
%   con - struct with control parameters (epsilon, sig, csq, c, gm)
%
% Default values are:
%   epsilon = 0.7
%   sig = 1;
%   csq = 0.1
%   gm = 50
%
% - Aditya

% Default parameters
defaults = struct(...
    'epsilon', 0.7, ...
    'sig',     1, ...
    'csq',     0.1, ...
    'gm',      50 ...
);

% Override defaults with user-specified values
user = struct(varargin{:});
fields = fieldnames(defaults);
for k = 1:numel(fields)
    if isfield(user, fields{k})
        defaults.(fields{k}) = user.(fields{k});
    end
end

% Compute derived quantity
defaults.c = sqrt(defaults.csq);

% Output control parameter struct
con = struct('epsilon', defaults.epsilon, ...
             'sig',     defaults.sig, ...
             'csq',     defaults.csq, ...
             'c',       defaults.c, ...
             'gm',      defaults.gm);
end