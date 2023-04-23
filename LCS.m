%least cost method of Transportation Problem
cost=[2 7 4;3 3 1;5 5 4;1 6 2];
% disp (cost);
a=[5 8 7 14];
% disp (a);
b=[7 9 18];
% disp (b);
% step 1 to check if balanced or not 
if sum(a)==sum(b)
    fprintf('Balanced\n');
else
    fprintf('Unbalanced\n');
    %in case of unbalanced dummy row or column is added 
    if sum(a)<sum(b) %add dummy row 
        cost(end+1,:)=zeros(1,length(b));
        a(end+1)=sum(b)-sum(a);
    elseif sum(b)<sum(a) %add dummy column 
        cost(:,end+1)=zeros(length(a),1);
        b(end+1)=sum(a)-sum(b);
    end
end
icost=cost;
x=zeros(size(cost));
[m,n]=size(cost);
bfs=m+n-1; %basic feasible solutions 
for i=1:size(cost,1) %for rows
    for j=1:size(cost,2) %for columns 
        hh=min(cost(:));
        [row_index,col_index]=find(hh==cost);
        x11=min(a(row_index),b(col_index));
        [val,ind]=max(x11); 
        ii=row_index(ind);
        jj=col_index(ind); 
%         y11=min(a(ii),b(jj)); 
        x(ii,jj)=val; 
        a(ii)=a(ii)-val; 
        b(jj)=b(jj)-val;
        cost(ii,jj)=inf;
   end
end

fprintf('Initial BFS= \n');
ib=array2table(x);
disp(ib);
totalbfs=length(nonzeros(x));
if totalbfs==bfs
    fprintf('Non Degenerate BFS \n');
else
    fprintf('Degenerate BFS \n');
end
initialcost=sum(sum(icost.*x));
fprintf('Cost= %d\n',initialcost);
