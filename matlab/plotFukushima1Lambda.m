function plotFukushima1Lambda(misc, labels)

for ii=1:length(labels)
  if labels(ii).name(end)=='*'
    labels(ii).name(end)=[];
  end
end



Tmax=misc.opt.Tmax;

% [val,ind]=sort(mean(misc.Sstat(1:120,:)),'descend');
%misc.Sstat=misc.Sstat(:,ind);
%misc.Lstat=misc.Lstat(:,ind);
col=[0.2 0.8 0.8;
    0.8 0.2 0.8;
    0.8 0.8 0.2;    
    1 0 0;
    0 1 0;
    0 0 1;
    0.5 0.5 0.5; 
    0.2 0.2 0.8;
    0.2 0.8 0.2;
    0.8 0.2 0.2   
    ];
q=0.01;
hold on;
minS=min(min((misc.Sstat)));
for c=1:size(misc.Lstat,2)
        fill([1:Tmax Tmax:-1:1]', [log10(log(2)./misc.Lstat(:,c))+sqrt((misc.Sstat(:,c))-minS)*q; flipud(log10(log(2)./misc.Lstat(:,c))-sqrt((misc.Sstat(:,c))-minS)*q)],col(c,:),'EdgeAlpha',0,'FaceAlpha',0.5);        
end
hold on;
for c=1:size(misc.Lstat,2)
  plot(1:Tmax, log10(log(2)./misc.Lstat(:,c)),'-','color',col(c,:),'LineWidth',0.5);
end
  ylim([-3 2])
grid on;

hold on;
for ii=1:length(labels)
  plot([Tmax-10 Tmax], [1;1]*log10(labels(ii).halflife), 'k-');
  text(Tmax+mod(ii,2)*50, log10(labels(ii).halflife), labels(ii).name);
end

xticks=get(gca,'xtick');
set(gca,'ytick',[-3 -2 -1 0 1 2])
[hx,hy]=format_ticks(gca,cellfun(@num2str, num2cell(xticks),'UniformOutput',0),{'10^{-3}','10^{-2}','10^{-1}','10^{0}','10^1','10^{2}'})
set(hx,'fontsize',16);
set(hy,'fontsize',16);

text(125,-3.4,'Iterations','FontWeight','bold','FontSize',20);
text(-20,-0.5,'Half lives (days)','FontWeight','bold','Rotation',90,'FontSize',20)
