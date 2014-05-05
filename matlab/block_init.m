function [block,opt]=block_init(data, opt)


% Number of samples for each location
nSample = cellfun(@length, getfieldarray(data, 'time'))';

% Find events
if isempty(opt.tev)
  tev = findeventsbydiffandgap(data,1/24);
else
  tev= opt.tev;
end

% Number of radioactive nuclides
N = opt.N;

% Number of events
E = length(tev);

% Number of observations
M = sum(nSample);

% Number of locations
L = length(data);

% Vectorize the observations
yy = zeros(M,1);
tt = zeros(M,1);
ll = zeros(M,1);
I  = cell(1,L);
ix0 = 0; 
for ii=1:L
  I{ii}=ix0+(1:nSample(ii)); ix0=I{ii}(end);
  yy(I{ii})=data(ii).data;
  tt(I{ii})=data(ii).time;
  ll(I{ii})=ii;
end

% Initialize variables
A = randcg(zeros(N,E), eye(N), eye(N), zeros(N,E),[],[], ones(N,E), 1);

%B = randcg(zeros(L,E), std(yy).^2*eye(L), eye(L), zeros(L,E), [], [], ones(L,E), 1);
B = randcg(zeros(L,E), max(1e-3,std(yy).^2)*eye(L), eye(L), zeros(L,E), [], ...
	   [], ones(L,E), 1);

bias = zeros(L,1);

sigma0 = max(1e-3,0.1*sum(yy.^2));
sigma = sigma0/M*ones(L,1);

% sigma_alpha = opt.sigma_alpha*ones(N,1);
sigma_beta  = opt.sigma_beta*ones(E,1);

sigma_bias  = opt.sigma_bias;

Tmask = tt(:,ones(1,E),ones(1,N))>tev(ones(M,1),:,ones(1,N));

nlogP=zeros(1,opt.Tmax);
nlpLoss=zeros(L,opt.Tmax);

res_rmse=zeros(1,opt.Tmax);

kk=0;

block = archive('kk','nSample','tev','E','M','L','yy','tt','ll','I','A','B','bias','sigma','sigma_beta','sigma_bias','Tmask','nlogP','nlpLoss','res_rmse');


if ~isempty(opt.ttime)
  Mtest = length(opt.ttime);
  block.Mtest = Mtest;
  block.tt_test=kron(ones(L,1), opt.ttime);
  block.ll_test=kron((1:L)', ones(Mtest,1));
  
  block.Tmask_test = block.tt_test(:,ones(1,E),ones(1,N))>tev(ones(Mtest*L,1),:,ones(1,N));

  block.test1 = zeros(Mtest,L);
  block.test2 = zeros(Mtest,L);
end  


