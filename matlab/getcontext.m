function context=getcontext()

context=dbstack(1,'-completenames');

if length(context)==0
  context=struct('file',pwd,'name','base','line',0);
end



