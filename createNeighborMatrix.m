function neighbors=createNeighborMatrix(rows,cols)
neighbors=zeros(8,rows,cols);
for i=1:rows
    for j=1:cols
        neighbors(:,i,j)=neighbors2d(i,j,rows,cols);
    end
end