function [ b ] = relaxation(a)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

b=a;
flag=0;
original=[100 100 100 100 100 100;
    100  95  95  85  65  50;
    100  80  80  75  60  50;
    100  70  70  60  55  50;
    100  40  40  35  30  50;
    100   0   0   0   0  50]

residual=zeros(6);

while(flag == 0)
    for i=2:5
        for j=2:5
            residual(i,j)=original(i-1,j)+original(i+1,j)+original(i,j-1)+original(i,j+1)-4*original(i,j);
        end
    end
    
    residual
    
    for i=2:5
        for j=2:5
            original(i,j)=original(i,j)+0.25*residual(i,j);
            original(i,j)=round(original(i,j));
        end
    end
    original
    flag=test(residual);
end
residual

    function [result] = test(A)
        result=1;
        
        for i=2:5
            for j=2:5
                if(abs(residual(i,j))>2)
                    result=0;
                end
            end
        end
        
    end




end

