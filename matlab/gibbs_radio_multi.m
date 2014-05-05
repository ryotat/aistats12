function [block,lambda_nucl,misc]=gibbs_radio_multi(data, lambda_nucl, varargin)

gitstring=gitlog;
stream = RandStream.getDefaultStream();
randstate=stream.State;

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'N', 10,...
                      'prior_nucl', [0.5 20],...
                      'prior_sigma', [0 1e+9],...
		      'prior_sigma_alpha', [0 1e+9],...
                      'prior_sigma_beta', [0 1e+9],...
                      'prior_sigma_bias', [0 1e+9],...
                      'sigma_alpha', 1,...
                      'sigma_beta', 1,...
                      'sigma_bias', 1,...
                      'Tmax', 100,...
                      'ttime', [],...
                      'tdata', [],...
                      'Ite', [],...
		      'xdate', [],...
                      'ktest',inf,...
                      'tev',[],...
                      'prop_variance',[],...
		      'init', []);

% Number of data blocks
nBlocks = length(data);

% Number of radioactive nuclides
N = opt.N;

if isempty(opt.init)
  for ii=1:nBlocks
    block{ii}=block_init(data{ii}, opt); 
  end

  % Initialize lambda_nucl
  if ~exist('lambda_nucl','var') || isempty(lambda_nucl)
    lambda_nucl=gamrnd(opt.prior_nucl(1)*ones(N,1),opt.prior_nucl(2)*ones(N,1));
  end

  % Initialize sigma_alpha
  sigma_alpha = opt.sigma_alpha*ones(N,1);


else
  lambda_nucl=opt.init.lambda_nucl;
  block = opt.init.block;
  sigma_alpha = opt.init.sigma_alpha;

  for ii=1:nBlocks
    [Mtest,L]=size(block{ii}.test1);
    block{ii}.kk    = 0;
    block{ii}.test1 = zeros(Mtest,L);
    block{ii}.test2 = zeros(Mtest,L);
  end
end

% Total number of events
nEvents = sum(cellfun(@(x)x.E, block));

nlogP=zeros(1,opt.Tmax);
nlogPmin = inf;
nlpB=zeros(nBlocks,opt.Tmax);

err_test=zeros(1,opt.Tmax);

Lstat = zeros(opt.Tmax, N);
Sstat = zeros(opt.Tmax, N);

% $$$ figure(1)
% $$$ clf;
% $$$ for ii=1:block{1}.L
% $$$   %  ax(ii)=subplot(block{1}.L,1,ii);
% $$$   ax(ii)=subplotxl(block{1}.L,1,[ii,1], 0.01, 0);
% $$$ end

% alpha_nucl1 = zeros(N,1);
% alpha_nucl2 = zeros(N,1);

for kk=1:opt.Tmax
  % Sample lambdas by random walk    
  fprintf('Sampling lambdas...\n');
  [alpha,block,prop_variance]=randomWalkMulti(log(lambda_nucl),block,...
           @(alpha, c, block)energylambda_c(alpha, c, block, sigma_alpha,opt));
  prop_variance
  lambda_nucl=exp(alpha);
  
  [mm,ix]=sort(-sigma_alpha);
  disp(lambda_nucl(ix))

 
%  if kk>opt.ktest
%    alpha_nucl1=alpha_nucl1+log(lambda_nucl);
%    alpha_nucl2=alpha_nucl2+log(lambda_nucl).^2;
%  end
  
  % nlogP value for lambda_nucl (nlogP(kk) initialize here)
  nlogP(kk)=sum(-(opt.prior_nucl(1)-1)*log(lambda_nucl)+lambda_nucl/opt.prior_nucl(2));
  
  % Sample each block
  res2=0;
  for ii=1:nBlocks
    fprintf('============= Sampling over block %d ========================\n',ii);
    block{ii}=block_sample(block{ii}, lambda_nucl, sigma_alpha, 1, opt);
    nlpB(ii,kk)=block{ii}.nlogP(block{ii}.kk);
    nlogP(kk) = nlogP(kk)+nlpB(ii,kk);
    res2 = res2+block{ii}.res2;
  end

  % Sample sigma_alpha
  fprintf('Sampling sigma_alpha...\n');
  if kk>0
  for jj=1:N
    A2=0;
    for ii=1:nBlocks
      A2=A2+sum(block{ii}.A(jj,:).^2);
    end
    sigma_alpha(jj)=1/gamrnd(opt.prior_sigma_alpha(1)+nEvents/2, opt.prior_sigma_alpha(2)/(1+opt.prior_sigma_alpha(2)*A2/2));
  end
  sigma_alpha
  end

  % nlogP value for A and sigma_alpha
  for jj=1:N
    A2=0;
    for ii=1:nBlocks
      A2=A2+sum(block{ii}.A(jj,:).^2);
    end
    nlogP(kk) = nlogP(kk)+(opt.prior_sigma_alpha(1)+nEvents/2+1)*log(sigma_alpha(jj))...
        +(1/opt.prior_sigma_alpha(2)+A2/2)/sigma_alpha(jj);
  end

  if nlogP(kk)<nlogPmin
    nlogPmin = nlogP(kk);
    mblock = cellfun(@(S)rmfield(S,{'Ecurv','Curv','zz','Tmask','Tmask_test'}),block,'UniformOutput',0);
    map=struct('block',{mblock},'lambda_nucl',lambda_nucl, ...
               'sigma_alpha',sigma_alpha,'nlogP',nlogP(kk),'kk',kk);
  end
  
  fprintf('[[%d]] Negative logP=%g; Sum of residuals=%g\n', kk, nlogP(kk), res2);

% $$$   figure(1);
% $$$   for ii=1:nBlocks
% $$$     % ax=subplot(nBlocks,1,ii);
% $$$     err_test(kk)=block_plot(ax, block{ii}, opt);
% $$$   end
  
% $$$   figure(2)
% $$$   plot([nlogP(1:kk); sum(block{1}.nlpLoss(:,1:kk),1)]','-x','linewidth',2);
% $$$   grid on;

%  figure(3)
%  Astat(kk,:)=sum(block{1}.A.^2);
%  Bstat(kk,:)=sum(block{1}.B.^2);
   Lstat(kk,:)=lambda_nucl';
   Sstat(kk,:)=sigma_alpha;
%  subplot(4,1,1);
%  semilogy(1:kk, Astat, 'linewidth',2);
%  subplot(4,1,2);
%  semilogy(1:kk, Bstat, 'linewidth',2);
%  subplot(2,1,1);
%  semilogy(1:kk, lambda_log(:,1:kk)', 'linewidth', 2);
%  subplot(2,1,2);
%  semilogy(1:kk, Sstat(1:kk,:), 'linewidth', 2);
  
end

% alpha_nucl_mean = alpha_nucl1/(kk-opt.ktest);
% alpha_nucl_std = sqrt(alpha_nucl2/(kk-opt.ktest)-alpha_nucl_mean.^2);

misc=struct('nlogP',nlogP,...
            'sigma_alpha',sigma_alpha,...
            'Lstat', Lstat,...
            'Sstat', Sstat,...
            'err_test', err_test,...
            'prop_variance', prop_variance,...
            'map', map,...
            'opt',opt,...
            'gitstring',gitstring,...
            'randstate',randstate);



function [energy, hessian, block]= energylambda_c(alpha,c,block,sigma_alpha,opt)

energy=sum(-(opt.prior_nucl(1)-1)*alpha+exp(alpha)/opt.prior_nucl(2));
hessian=exp(alpha)/opt.prior_nucl(2);

for ii=1:length(block)
  [eng, hess, block{ii}]=block_energylambda_c(alpha,c,block{ii},sigma_alpha,opt);
  energy  = energy + eng;
  hessian = hessian + hess;
end
