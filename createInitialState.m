function s=createInitialState(phi,size)

s=double(rand(size)<phi)+1;


sites=size(1)*size(2);

psites=round(phi*sites);
r=(s==2);
q=~r;
p1=nnz(r);
while p1~=psites
    
if p1>psites
    r1=find(r);
    s(r1(ceil(length(r1)*rand(p1-psites,1))))=1;
elseif p1<psites
    r1=find(q);
    s(r1(ceil(length(r1)*rand(psites-p1,1))))=2;
end
r=(s==2);
q=~r;
p1=nnz(r);
end


% s=ones(size);
% tic
% while nnz(s==2)<psites
%     s(ceil(sites*rand))=2;
% end
% toc
% nnz(s==2)/sites