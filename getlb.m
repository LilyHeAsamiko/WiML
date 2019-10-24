function Y = getlb(I11, lt, a, b, winL1, winL2)   
    Y = zeros(a, b); 
    for i = 1:winL1:(a-1)*winL1+1
        for j = 1:winL2:(b-1)*winL2+1
            batch = I11(i:i+winL1-1, j:j+winL2-1);
            if sum(sum(batch >= lt))/(winL1*winL2)> 0.5
                Y(uint32((i+winL1)/winL1), uint32((j+winL1)/winL2)) = 1;
            else
                Y(uint32((i+winL2)/winL1), uint32((j+winL2)/winL2)) = 0;                
            end
        end
    end
end