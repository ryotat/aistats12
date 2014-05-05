RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

load fukushima1.mat

Tmax=1000;

tinterval = datenum(2011,3,11)+[0,30];
% Convert into struct
for ii=1:length(I)
  time = Datenum(J==ii);
  numb = Numbers(J==ii);
  
  ix=find(time>=tinterval(1) & time<=tinterval(2));
  
  data(ii)=struct('location',uq{ii},...
                  'locationeng',uqeng{ii},...
                  'time', time(ix),...
                  'data', numb(ix));
end

% Remove location '-'
data(1)=[];

% Take only the locations that have at least ten measurements
nSample = cellfun(@length, getfieldarray(data, 'time'));
ix = find(nSample<10);
data(ix)=[];
nSample(ix)=[];

% 
nte = 1000;
ttime = linspace(datenum(2011,3,11), datenum(2011,4,10), nte)';


%
% $$$ tev = [datenum(2011,3,12,7,0,0),...   % releasing pressure
% $$$        datenum(2011,3,12,15,36,0),... % after-shake and white smoke #1
% $$$        datenum(2011,3,13,9,20,0),...  % manual vent
% $$$        datenum(2011,3,14,11,1,0),...  % #3 explosion
% $$$        datenum(2011,3,15,6,10,0),... % #2 (and also #4?) explosion (according to report#24)
% $$$        datenum(2011,3,15,9,38,0),... % #4 fire
% $$$        datenum(2011,3,16,5,45,0),... % #4 fire
% $$$        datenum(2011,3,16,8,30,0),... % #3 smoke (report #26)
% $$$        ];

tev = [];

% Gibbs Sample 
[block,lambda_nucl_est,misc]=gibbs_radio_multi({data}, [], ...
                                               'tev', tev,...
                                    'N', 10,...
                                    'Tmax', Tmax,...
                                    'ttime', ttime,...
                                    'ktest', 500);
len=length(data)

misc.context = getcontext;

save 'result_fukushima1_1hour_1000.mat' block lambda_nucl_est misc

% Visualize Prediction
figure
set(gcf,'position',[-422, 243, 1358, 330]);

xtick = datenum(2011,3,11)+(0:30);
ticklabel = datestr(xtick, 'mm/dd'); ticklabel(:,1)=[];

for ii=1:len
  clf;
  ax=nan*ones(1,len);
  ax(ii)=gca;
  block_plot(ax, block{1}, misc.opt);
  xlim([min(xtick), max(xtick)]);
  set(gca,'xtick',xtick,...
        'xticklabel',ticklabel);
  text(0,0.8,data(ii).locationeng,'unit','normalized');
  grid on;
  set(gcf,'paperpositionmode','auto');
  print('-dpng',sprintf('fukushima1-%d.png', ii));
end

% Plot observation RMSE
figure, plot(block{1}.res_rmse);
ylim([250 300]);
grid on;
title('Observation RMSE');
xlabel('Iterations');
ylabel('RMSE');

% Plot lambdas
load halftimes
for ii=1:length(name)
  labels(ii)=struct('name',name{ii},'halflife',half(ii));
end
figure, plotFukushima1Lambda(misc, labels);

