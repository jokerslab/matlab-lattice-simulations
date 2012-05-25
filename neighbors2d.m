function n=neighbors2d(i,j,rows,columns)
%neighbors2d creates a matrix whose entries return a vector of nearest
%neighbor indices
        if nargin<4
            columns=rows;
            rows=j;
            index=i;
            
            [i j]=ind2sub([rows columns],index);
        else
        index=(j-1)*rows+i;
        end
        n=zeros(1,8);
        if (i~=1) && (i~=rows)
            n(1:2)=[index-1 index+1];
            if (j~=1) && (j~=columns)
                n(5:8)=[index-rows-1 index+rows-1 index-rows+1 index+rows+1];
            elseif j==1
                n(5:8)=[index+(columns-1)*rows-1 index+rows-1 index+(columns-1)*rows+1 index+rows+1];
            else
                n(5:8)=[index-rows-1 index-(columns-1)*rows-1 index-rows+1 index-(columns-1)*rows+1 ];
            end 
            
        elseif i==1
            n(1:2)=[index+rows-1 index+1];
            if (j~=1) && (j~=columns)
                n(5:8)=[index-1 index+rows+rows-1 index-rows+1 index+rows+1];
            elseif j==1
                n(5:8)=[index+(columns-1)*rows+rows-1 index+rows+rows-1 index+(columns-1)*rows+1 index+rows+1];
            else
                n(5:8)=[index-1 index-(columns-1)*rows-1+rows index-rows+1 index-(columns-1)*rows+1 ];
            end 
        else
            n(1:2)=[index-1 index-rows+1];
            if (j~=1) && (j~=columns)
                n(5:8)=[index-rows-1 index+rows-1 index-2*rows+1 index+1];
            elseif j==1
                n(5:8)=[index+(columns-1)*rows-1 index+rows-1 index+(columns-1)*rows+1-rows index+1];
            else
                n(5:8)=[index-rows-1 index-(columns-1)*rows-1 index-2*rows+1 index-(columns-1)*rows+1-rows ];
            end 
        end
        if (j~=1) && (j~=columns)
            n(3:4)=[index-rows index+rows];
        elseif j==1
            n(3:4)=[index+(columns-1)*rows index+rows];
        else
            n(3:4)=[index-rows index-(columns-1)*rows];
        end
        
        
    end