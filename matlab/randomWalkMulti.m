function [x,block,prop_variance]=randomWalkMulti(x,block,fun,varargin)

N=length(x);
opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'iterations', 5,...
                      'prop_variance', ones(N,1),...
                      'varfact',1);

fprintf('Computing initial energy...\n');
[nlogP_old, hess, block] = fun(x,[],block);

prop_variance=opt.varfact./abs(hess);

accept=0;
accept_l=zeros(N,1);
for l = 1:opt.iterations
    disp(['Operating on loop ' num2str(l) '/' num2str(opt.iterations)]);    
    for c=1:N        
        xnew=x;
        xnew(c)=xnew(c)+sqrt(prop_variance(c))*randn;
        block_new = block;
        %for jj=1:length(block)
        %block_new{jj}.A(c,:)=exp(xnew(c)-x(c))*block{jj}.A(c,:);
        %end
        
%        fprintf('Computing energy for iter=%d nucl=%d...\n',l,c);
        [nlogP_new, hess_new, block_new] = fun(xnew,c,block_new);
        prop_variance_new=opt.varfact./abs(hess_new);
       
        d2=(xnew-x).^2;
%        fprintf('nlogP_old=%g nlogP_new=%g\n', nlogP_old, nlogP_new);
        if rand<exp(-(nlogP_new...
                      +0.5*sum(d2./prop_variance_new+log(prop_variance_new))...
                      -nlogP_old...
                      -0.5*sum(d2./prop_variance+log(prop_variance))))
            
          x=xnew;
            block=block_new;
            nlogP_old=nlogP_new;
            prop_variance=prop_variance_new;
            accept=accept+1;
            accept_l(c)=accept_l(c)+1;
            
        end        
    end   

end

disp(['Accepted ' num2str(accept) '/' num2str(opt.iterations*N)]);

  