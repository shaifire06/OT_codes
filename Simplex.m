%simplex
clc;
close all;
% Input Parameters
C = [2 5];%objective function
a = [1 4; 3 1; 1 1];%information
b = [24; 21; 9];

[m,n] = size(a); % n is the number of variables 
s = eye(m);% s=[1 0 0;0 1 0;0 0 1]
A = [a s b];%the format in which table is written

cost = zeros(1,n+m+1);%cost=[0 0 0 0 0 0]
cost(1:n) = C; %in the cost we only need three zeros for s1,s2,s3

bv = n+1:1:n+m; %basic variables 

zjcj = cost(bv)*A - cost;
zcj = [zjcj; A];

simplexTable=array2table(zcj);
simplexTable.Properties.VariableNames(1:size(zcj,2))={'x1','x2','s1','s2','s3','Sol'}

flag=true;

while flag
    if any(zjcj<0)
        fprintf("The Current BFS is not Optimal \n");
        zc = zjcj(1:end-1);
         [Enter_val, pvt_col] = min(zc);
    if all (A(:,pvt_col)<=0)
        error("LPP is Unbounded all entries are <=0 in column %d",pvt_col);
    else
        sol = A(:,end)
        column = A(:,pvt_col)
        for i = 1:m
            if column(i) > 0
                ratio(i) = sol(i)./column(i);
            else
                ratio(i) = inf;
            end
        end
        [leaving_val, pvt_row] = min(ratio);
    end

    bv(pvt_row) = pvt_col;%new basic variable 
    pvt_key = A(pvt_row, pvt_col);

    A(pvt_row,:) = A(pvt_row,:)./pvt_key;

    for i = 1:m
        if (i~=pvt_row)
            A(i,:) = A(i,:) - A(i, pvt_col).*A(pvt_row,:);
        end
    end

    zjcj = zjcj - zjcj(pvt_col).*A(pvt_row,:);
    zcj = [zjcj;A];

    table = array2table(zcj);

    table.Properties.VariableNames(1:n+m+1) = {'x1','x2','s1','s2','s3','Sol'}

    else
        flag = false;
        fprintf("The Current BFS is Optimal \n")
    end
end
