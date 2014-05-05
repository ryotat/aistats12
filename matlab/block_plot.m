function err=block_plot(ax, block, opt)

opt=set_defaults(opt, 'plot_tev',0);

err=nan;

extract(block);

fprintf('In block_plot kk=%d\n', kk);


if kk>opt.ktest
  test_mean = test1/(kk-opt.ktest);
  test_std = sqrt(test2/(kk-opt.ktest)-test_mean.^2);

  if ~isempty(opt.Ite)
    err=sum(compute_error(test_mean, opt.tdata, opt.Ite));
  end
end

for jj=1:L
  if ~isnan(ax(jj))
  axes(ax(jj));
  if ~isempty(opt.ttime) && kk > opt.ktest
    plot(tt(I{jj}), yy(I{jj}), 'x',...
         opt.ttime, test_mean(:,jj),'-','linewidth',2);
    hold on;
    if ~isempty(opt.tdata)
      plot(opt.ttime, opt.tdata(:,jj),'m-.','linewidth',2);
    end
    
    plot(opt.ttime, test_mean(:,jj)-test_std(:,jj), '--',...
         opt.ttime, test_mean(:,jj)+test_std(:,jj), '--',...
         'color',[1 0 0]);

    if opt.plot_tev
      plot([1;1]*tev, (ylim)'*ones(1,length(tev)), 'm--', 'linewidth', 2);
    end
    
    hold off;
    
    if ~isempty(opt.xdate)
      xlim([min(opt.xdate),max(opt.xdate)]);
      grid on;
      set(gca,'xtick', opt.xdate,...
	      'xticklabel',cellfun(@(x)datestr(x,'dd/m'),num2cell(opt.xdate),...
                                   'UniformOutput',0));
  
    end
    
  else
    plot(tt(I{jj}), yy(I{jj}),...
         tt(I{jj}), zall(I{jj}), 'linewidth', 2);
  end
  end
end
