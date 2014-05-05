function tev = findeventsbydiffandgap(data, delta)

if ~exist('delta','var')
  delta=1/24;
end

tev=[];

for ii=1:length(data)
  M=max(data(ii).data)-min(data(ii).data);
  time = data(ii).time;
  ixdiff = find(diff(time)<1 & diff(data(ii).data)>M*0.02);
  ixgap  = find(diff(time)>1)+1;

  tev = [tev, time([ixdiff; ixgap])'];
end
tev=unique(tev);

% Remove overlaps
dtev = diff(tev);
while min(dtev)<delta
  ix = find(dtev < delta);

  tev(ix) = (tev(ix)+tev(ix+1))/2;
  tev(ix+1)=[];
  
  dtev = diff(tev);
end


