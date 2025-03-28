% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function index_pairs = fcnmatchSURF(f1,f2,rangeRatio)
%f1 = Nx128;
%f2 = Mx128;
dotProdRatio = sqrt(rangeRatio); %approximation to go from true distance metric to dot product metric
N = size(f1,1);

%NORM L2 WAY --------------------------------------------------------------
match = zeros(1,N);
dotprods = f2*f1'; %=cos(theta)
[m1, i] = max(dotprods);  for j=1:N; dotprods(i(j),j) = 0; end;  m2= max(dotprods);
ratio =  acos(m1)./acos(m2);
passflag = ratio < dotProdRatio;


if sum(passflag)<100 %not enough matches
    j=0;
    while sum(passflag)<100 && j<10
        j=j+1;
        dotProdRatio = dotProdRatio*1.1;
        passflag = ratio < dotProdRatio;
    end
end


match(passflag) = i(passflag);
v1 = find(match>0);
v2 = match(v1);
index_pairs = [v1' v2'];


if numel(v2)~=numel(unique(v2)) %we have duplicates in v2
    s = sortrows([index_pairs, ratio(passflag)'], 3);
    [~,i]=unique(s(:,2),'rows','first');
    index_pairs = s(i,1:2);
end


% c=find(match>0);
% dotProdRatio = dotProdRatio(c(:))';
% mp=[loc1(c(:),1) loc1(c(:),2) loc2(match(c(:)),1) loc2(match(c(:)),2) dotProdRatio];
% mp_sorted=sortrows(mp,5); %sort on valsRatio (ascending)
% [mp_unique,index]=unique(mp_sorted(:,1:4),'rows','first');
% matchpoints = sortrows([mp_unique mp_sorted(index,5)],5); % unique positions, and sorted by distRatio
end

