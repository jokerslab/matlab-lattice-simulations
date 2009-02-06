%Plot Lattice
%Dan Jung
function a=plotBonds(sites,varargin)
D=ndims(sites);
if length(varargin)<1
    vals=[-1 1];
    colors={'black' 'none'};
elseif length(varargin)<2
    vals=varargin{1};
    %colors={'black' 'red' 'green' 'cyan' 'blue' 'magenta' 'yellow' [.5 .5 .5]};
    
    cmap=colormap('cool');

    for i=1:length(vals)
        colors{i}=cmap(floor((i-1)*64/(length(vals)-0.999)+1),:);
    end
    %colors={[.5 0 0] [1 0 0] [1 .5 .5] [.5 .5 .5] [.5 1 .5] [0 1 0] [0 .5 0]};
else
    vals=varargin{1};
    colors=varargin{2};
end

if D==2
    [r,c]=size(sites);
    d=max(r,c);
    msize=250/d;
    for i=1:length(vals)
        [y,x]=find(sites==vals(i));
            y=-y;
        p=mod(y,2);  %1 on odd row--downward bond
        q=1-p;  %1 on even row --right bond
        d=p+2*q;
        y=(y-p)/2;
        x=[x';x'+q'];
        y=[y';y'-p'];
        z=zeros(2,length(x));
        
        line(x,y,z,'color',colors{i},'linewidth',4)
    end
    
end

