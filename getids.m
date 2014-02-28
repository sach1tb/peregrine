function ids=getids(k, Xh, Xi)
if ~isempty(Xh)
    ids=find(Xh(1:Xi.nX:end,end)~=0);
else
    ids=[];
end
if isempty(ids), ids=[]; end