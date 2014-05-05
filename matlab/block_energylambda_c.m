function [energy, hess, block]= block_energylambda_c(alpha,c,block,sigma_alpha,opt)

extract(block);
sigma=block.sigma;
in = mat2cell([sigma, nSample],ones(L,1),2);
sigmaall = cell2mat(cellfun(@(x)x(1)*ones(x(2),1), in, 'UniformOutput',0));


N=length(alpha);
lambda_nucl = exp(alpha);

if ~isempty(c)
    %% Compute the basis function exp(-lambda_nucl(t-te))
    Curv(:,:,c) = -lambda_nucl(c)*max(0,tt(:,ones(1,E))-tev(ones(M,1),:));
    Ecurv(:,:,c) = Tmask(:,:,c).*exp(Curv(:,:,c));

    %% Current prediction (per event and per nuclide)
    zz(:,:,c) = (permute(A(c,:,ones(1,M)),[3,2,1]).*B(ll,:)).*Ecurv(:,:,c);
else         
    % Compute the basis function exp(-lambda_nucl(t-te))
    Curv = -repmat(shiftdim(lambda_nucl,-2), [M, E, 1])...
           .*max(0,tt(:,ones(1,E),ones(1,N))-tev(ones(M,1),:,ones(1,N)));
    Ecurv = Tmask.*exp(Curv);    
    zz = permute(A(:,:,ones(1,M)),[3,2,1]).*B(ll,:,ones(1,N)).*Ecurv;
end

% Current prediction (all events)
zall = sum(sum(zz,2),3)+bias(ll);

energy=0.5*sum((yy-zall).^2./sigmaall);

% nlogP value for A
N=size(A,1);
for jj=1:N
  eng=sum(A(jj,:).^2)/(2*sigma_alpha(jj));
  
  energy = energy+eng;
end


Curvzz = Curv.*zz;
hess=-permute(sum(Curvzz.*(1+Curv),2),[3,1,2])*((yy-zall)./sigmaall)...
     +(permute(sum(Curvzz,2),[3,1,2]).^2)*(1./sigmaall);



block.Curv  = Curv;
block.Ecurv = Ecurv;
block.zz    = zz;
