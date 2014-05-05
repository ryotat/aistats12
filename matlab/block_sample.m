function block=block_sample(block, lambda_nucl, sigma_alpha, Tmax, opt)

extract(block);
sigma=block.sigma;

N=length(lambda_nucl);

% $$$ muAsave = zeros(E+1,N);
% $$$ vvAave = zeros(E+1,N);
% $$$ muBsave = zeros(E+1,L);
% $$$ vvBave = zeros(E+1,L);

while kk<block.kk+Tmax
  kk=kk+1;

  in = mat2cell([sigma, nSample],ones(L,1),2);
  sigmaall = cell2mat(cellfun(@(x)x(1)*ones(x(2),1), in, 'UniformOutput',0));

  zz = sum(zz,3);
  zall = sum(zz,2)+bias(ll);


  
  for rep=1:5
    % Sample A
    fprintf('Sampling A...\n');
    for jj=1:E
      res = yy-zall+zz(:,jj);
      T = bsxfun(@times, B(ll,jj), squeeze(Ecurv(:,jj,:)));
      Tdsigma=bsxfun(@rdivide,T,sigmaall);
      C = (Tdsigma'*T)+diag(1./sigma_alpha);
      mu = C\(Tdsigma'*res);
      
      Cinv = inv(C);
 %     muAsave(jj,:)=mu';
%      vvAsave(jj,:)=diag(Cinv)';
      
      A(:,jj)=randcg(mu, Cinv, eye(N), zeros(N,1), [], [], A(:,jj),1);
      
      
      % Current prediction
      zz(:,jj) = (permute(Ecurv(:,jj,:),[1,3,2])*A(:,jj)).*B(ll,jj);
      zall = sum(zz,2)+bias(ll);
    end
    
    % Sample B
    fprintf('Sampling B (and the bias term)...\n');
    for jj=1:L
      % Sum over nuclides
      T=[sum(permute(A(:,:,ones(1,nSample(jj))),[3,2,1]).*Ecurv(I{jj},:,:),3),ones(nSample(jj),1)];
      C=T'*T/sigma(jj)+diag(1./[sigma_beta; sigma_bias]);
      mu = C\(T'*yy(I{jj})/sigma(jj));

      Cinv = inv(C);
%      muBsave(:,jj)=mu;
%      vvBsave(:,jj)=diag(Cinv);
     
      vv=randcg(mu, Cinv, eye(E+1), zeros(E+1,1), [], [], [B(jj,:),bias(jj)]',1);
      B(jj,:) = vv(1:end-1)';
      bias(jj) = vv(end);
    end
    
    % Update current prediction
    zz = sum(permute(A(:,:,ones(1,M)),[3,2,1]).*B(ll,:,ones(1,N)).*Ecurv,3);
    zall = sum(zz,2)+bias(ll);

      % Sample sigma
  res = yy-zall;
  fprintf('Sampling sigma\n');  
  for jj=1:L
    res2 = sum(res(I{jj}).^2);
    sigma(jj)=1/gamrnd(opt.prior_sigma(1)+nSample(jj)/2, 1/(1/opt.prior_sigma(2)+res2/2));
  end

  end


  % Sample sigma_beta
  fprintf('Sampling sigma_beta\n');
  sigma_beta=1/gamrnd(opt.prior_sigma_beta(1)+L*E/2, 1/(1/opt.prior_sigma_beta(2)+sum(B(:).^2)/2));
  sigma_beta=sigma_beta(ones(E,1));
  
  % Sample sigma_bias
  fprintf('Sampling sigma_bias\n');
  sigma_bias=1/gamrnd(opt.prior_sigma_bias(1)+L/2, 1/(1/opt.prior_sigma_bias(2)+sum(bias.^2)/2));
  
  
  % Compute prediction
  if exist('Mtest','var') && kk > opt.ktest
  Curv_test = -repmat(shiftdim(lambda_nucl,-2), [Mtest*L, E, 1])...
      .*max(0,tt_test(:,ones(1,E),ones(1,N))-tev(ones(Mtest*L,1),:,ones(1,N)));
  Ecurv_test = Tmask_test.*exp(Curv_test);    
  ztest=zeros(Mtest,L);
  ztest(:)=sum(sum(permute(A(:,:,ones(1,Mtest*L)),[3,2,1]).*B(ll_test,:,ones(1,N)).*Ecurv_test,3),2)+bias(ll_test);
  
  test1 = test1+ztest;
  test2 = test2+ztest.^2;
  end

  % nlogP value for the likelihood and sigma (noise variance)
  nlogP(kk)=0;
  res = yy-zall;
  for jj=1:L
    nlogP(kk)=nlogP(kk)+(opt.prior_sigma(1)+nSample(jj)/2+1)*log(sigma(jj))...
              +(1/opt.prior_sigma(2)+sum(res(I{jj}).^2)/2)/sigma(jj);
    nlpLoss(jj,kk)=nSample(jj)*log(sigma(jj))/2;
  end


  % nlogP value for B and sigma_beta
  for jj=1:E
    nlogP(kk) = nlogP(kk) + (opt.prior_sigma_beta(1)+L/2+1)*log(sigma_beta(jj))...
        +(1/opt.prior_sigma_beta(2)+0.5*sum(B(:,jj).^2))/sigma_beta(jj);
  end

  % nlogP value for bias
  nlogP(kk) = nlogP(kk) + (opt.prior_sigma_bias(1)+L/2+1)*log(sigma_bias)...
      +(1/opt.prior_sigma_bias(2)+sum(bias.^2)/2)/sigma_bias;
  
  fprintf('[%d] Negative logP=%g; Sum of residuals=%g\n', kk, nlogP(kk), sum((yy-zall).^2));

  % fprintf('Sum of residuals=%g\n', sum((yy-zall).^2));
  
  
  res_rmse(kk) = sqrt(mean(res.^2));
end


res2 = sum(res.^2);

block = updatestruct(block, 'kk','A','B','bias','sigma','sigma_bias','sigma_beta','nlogP','nlpLoss','res2','zall','res_rmse');

if exist('Mtest','var')
  block.test1     = test1;
  block.test2     = test2;
end

% $$$ block.muAsave = muAsave;
% $$$ block.vvAsave = vvAsave;
% $$$ block.muBsave = muBsave;
% $$$ block.vvBsave = vvBsave;
