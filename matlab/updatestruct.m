function S=updatestruct(S, varargin)

for i=1:length(varargin)
  name =varargin{i};
  S = setfield(S, name, evalin('caller', name));
end
