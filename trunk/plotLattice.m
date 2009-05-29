%Plot Lattice
function a=plotLattice(sites,varargin)
cla
D=ndims(sites);
if length(varargin)<1
    vals=[-1 1];
    colors={'black' 'none'};
elseif length(varargin)<2
    vals=varargin{1};
    colors={'none' 'black' 'green' 'cyan' 'blue' 'magenta' 'yellow' [.5 .5 .5]};
    if length(vals)>8
        
        cmap=colormap('lines');
        for i=1:length(vals)
            im=mod(i,64)+1;
            colors{i}=cmap(im,:);
        end
    end
else
    vals=varargin{1};
    colors=varargin{2};
end

if D==2
    [r,c]=size(sites);
    d=max(r,c);
    if d>60
       colormap([1 1 1;0 0 0;0 1 0;0 0 1])
       image(sites); 
    else
    msize=max(250/d,4);
    for i=1:length(vals)
        [y,x]=find(sites==vals(i));
            y=-y;

        line('xdata',x,'ydata',y,'markersize',msize,'marker','o', ...
    'linestyle','none','markerfacecolor',colors{i},'markeredgecolor','black')
    end
    axis([.5 c+.5 -r-.5 -.5])
    axis equal
    end
end

