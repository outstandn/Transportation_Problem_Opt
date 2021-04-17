%Written by: Tobechukwu Onyedi all rights reserved
%Transportation Problem: North West Corner and UV Optimization (Special
%Case)
%No part of this code is allowed to be used.
%if any part of the code must be used, reference should be stated.
%%Defining Problem and initializing solution and back up

clc; clear var; clear all; close all;

%transportationMatrix = [1 4 10 7 9 14;9 5 17 16 12 16;8 6 14 10 17 11;11 8 20 10 15 16;3 15 18 2 19 0];
% transportationMatrix = [4	 6	 1	 16	 10	 18; 9	 8	 6	 12	 13	 12;3	 11	 5	 9	 9	 9; 10	 13	 4	 13	 10	 22;14	 11	 16	 12	 8  0];
% transportationMatrix = [1 3 4 5 20; 5 2 10 3 30;3 2 1 4 50; 6 4 2 6 20;0 0 0 0 5;30 20 60 15 0];
%transportationMatrix = [13 22 19 21 16 20 1; 18 17 24 18 22 27 1; 20 22 23 24 17 31 1;...
%                        14 19 13 30 23 22 1; 21 14 17 25 14 23 1; 17 23 18 20 16 24 1; 1 1 1 1 1 1 0];
%transportationMatrix = [32 60 200 20; 40 68 80 30; 120 104 60 45; 30 35 30 0];
%transportationMatrix = [90 100 130 20; 100 140 100 15; 100 80 80 10; 5 20 20 0];
transportationMatrix = [10 0 20 11 20; 12 7 9 20 25; 0 14 16 18 15; 10 15 15 20 0];

                    
backUpTransportationMatrix = transportationMatrix; %get size of intial Transportation problem
[row, column]=size(transportationMatrix); %get the size of the matrix
disp(transportationMatrix);


CostMatrix = transportationMatrix(1:row-1,1:column-1);
rowSupplyMatrix = transportationMatrix(1:row-1,column);
colDemandMatrix = transportationMatrix(row,1:column-1);

backUpRowSupplyMatrix = rowSupplyMatrix; %backup supply matrix
backUpColDemandMatrix = colDemandMatrix; %backup demand matrix

backUpCostMatrix = CostMatrix; %backup initial cost matrix
leastCostMatrix = zeros(row-1,column-1); %create leastcost matrix with same size as cost matrix
leastCostMatrix(1:row-1,1:column-1) = -1; %replace matrix with -1

U = zeros(row-1,1); U(2:row-1)= intmax; %create U zero matrix and set to intmax except U(0)
V = zeros(1,column-1); V(:)= intmax; %create V zero matrix and set to intmax


% Check if balance or not then calculate initia feasible solution
supplySum=sum(backUpRowSupplyMatrix);
demandSum=sum(backUpColDemandMatrix);
disp('-------------------------------------------------------------------');
%% Calculate initial feasible solution
if (supplySum == demandSum)
    disp('Balanced Transportation Problem.');
    disp('Solving Using Least Cost Method.');
    while (supplySum ~= 0) || (demandSum ~= 0)
        
        %get minimum cost and the row, column index
        minCostColumn = min(CostMatrix(CostMatrix > 0),[],1);
        tempminCost = min(minCostColumn);
        [Min,Indexs] = min(CostMatrix(:));
        [minRow, minCol] = ind2sub(size(CostMatrix),Indexs);
        
        %use index to get minimum of supply and demand
        tempGoods = min(rowSupplyMatrix(minRow), colDemandMatrix(minCol));
        
        %update the solution matrix
        leastCostMatrix(minRow,minCol) = tempGoods;
        
        %set the cost matrix to max value so that it is not read again.
        CostMatrix(minRow, minCol) = intmax;
        
        %update the remainder of supply and check if "zero"
        rowSupplyMatrix(minRow)= rowSupplyMatrix(minRow)- tempGoods;
        if rowSupplyMatrix(minRow) == 0 %delete row if "zero"
            CostMatrix(minRow,1:column-1)=intmax;
        end
        
        %update the remainder of demand and check if "zero"
        colDemandMatrix(minCol)= colDemandMatrix(minCol)- tempGoods;
        if colDemandMatrix(minCol) == 0 %delete column if "zero"
            CostMatrix(1:row-1,minCol)= intmax;
        end
        
        %recalculate sum of supply and demand (condition to stop loop).
        supplySum=0; demandSum=0;
        for i=1:row-1
            supplySum=supplySum + rowSupplyMatrix(i);
        end
        for j=1:column-1
            demandSum=demandSum + colDemandMatrix(j);
        end
        disp(leastCostMatrix)
    end
    LleastCostMatrix(1:row-1,1:column-1) = leastCostMatrix.*(leastCostMatrix >-1);
else
    disp('Unbalanced Transportation Problem');
end
transportationMatrix = backUpTransportationMatrix;
%% Compute Transportation Cost
disp(LleastCostMatrix);
transportationCost=0;
for i=1:row-1
    for j=1:column-1
        transportationCost = transportationCost+(transportationMatrix(i,j).*LleastCostMatrix(i,j));
    end
end

disp(['The transportation cost is: ',num2str(transportationCost)]);
backUpLeastCostMatrix = LleastCostMatrix; %Back up initial solution
%% Optimization
disp('-------------------------------------------------------------------');
disp('---------------------Starting Optimization!------------------------');
disp('-------------------------------------------------------------------');
iteration = 0; ifailed = 0; dummyRow = 0; dummyColumn = 0;
CostMatrix = backUpCostMatrix;

while (true)
    %Reset Values of U,V and P
    U = zeros(row-1,1);
    U(2:row-1)= intmax; %create U zero matrix and set to intmax except U(1)

    V = zeros(1,column-1);
    V(1:column-1)= intmax; %create V zero matrix and set to intmax
    %Check for accuracy of solution
    AB= sum(sum (leastCostMatrix > -1));
    ABC = sum(backUpColDemandMatrix > -1);
    ABD = sum(backUpRowSupplyMatrix > -1);
    maxRC= max(row,column) -1;

    if (AB == (ABC + ABD-1))
        disp('Northwest solution fulfilled the condition Row + Column - 1');
        disp('Calculating "U" & "V"');
        
        %Calculate U & V
        U(1) = 0;
        while (true)
            for i=1:row-1
                for j=1:column-1
                    if ((leastCostMatrix(i,j) > -1) && (U(i) == intmax) && (V(j) == intmax))
                        continue
                    elseif ((leastCostMatrix(i,j) > -1) && (U(i) == intmax))
                        U(i) = CostMatrix(i,j) - V(j);
                    elseif ((leastCostMatrix(i,j) > -1) && (V(j) == intmax))
                        V(j) = CostMatrix(i,j) - U(i);
                    end                    
                end
            end
            %check if intmas is disappeared from UV arrays the break
            if isempty(U(U == intmax)) && isempty(V(V == intmax))
                disp('U is: '); disp(U);
                disp('V is: '); disp(V);
                break;
            end
        end 
        %Calculate Penaltiy
        disp('Calculating "Penalty"');
        P = zeros(row-1, column-1); %Penalty matrix
        for i=1:row-1
            for j=1:column-1
                if leastCostMatrix(i,j) == -1				
					%Used to calculate Penaltiy
                    P(i,j) = U(i) + V(j) - transportationMatrix(i,j); 
                end
            end
        end
        disp('Penalty is: '); disp(P);

		%Find Maximum P and its Indexes
		[pMax,Ind] = max(P(:));
		[pMaxRow, pMaxCol] = ind2sub(size(P),Ind);
        fprintf('Maximum Penalty is: %d at index [%d, %d]\n', pMax,pMaxRow, pMaxCol);
        

		%If PMax is greater than Zero then some optimization can be done
        if pMax > 0
            %get a logical (1,0) array of element that is > -1
            loop = leastCostMatrix > -1;
            loop(pMaxRow, pMaxCol) = 1; %input loop starting point for optimization
            while(true)
                %Sum this 1 and 0 row wise and column wise
                loopRow = sum(loop,2);            
                %Eliminate rows that has a sum of 1, it signifies there no path continuity
                for i=1:row-1
                    if loopRow(i) == 1
                        loop(i,1:column-1)= 0;
                        loopRow(i) = 0;
                    end
                end
                
                %Sum this 1 and 0 column wise
                loopColumn = sum(loop,1);
                %Eliminate columns that has a sum of 1,it signifies there no path continuty
                for j=1:column-1
                    if loopColumn(j) == 1
                        loop(1:row-1,j)= 0;
                        loopColumn(j) =0;
                    end
                end
                
          % Check again if there are still failed paths using modulus 2
                loopRow = sum(loop,2);
                loopColumn = sum(loop,1);
                modRow = mod(loopRow, 2)==1;
                modColumn = mod(loopColumn, 2)==1;
                
                if isempty(modRow(modRow == 1))&& isempty(modColumn(modColumn == 1))
                    break;
                end
            end
            
            disp('path found is as shown below');
            disp(loop);

            LleastCostMatrix = zeros(row,column);
            LleastCostMatrix(1:row-1,1:column-1) = leastCostMatrix(1:row-1,1:column-1);
            
            temp = LleastCostMatrix;
            LleastCostMatrix = leastCostMatrix; %save leastCostMatrix
            leastCostMatrix = temp;

            %Get Minimum element and index from created loop
            newCost= loop(1:row-1,1:column-1).* leastCostMatrix(1:row-1,1:column-1);
            tempMinNewCost = min(newCost(newCost > 0));
            leastCostMatrix = LleastCostMatrix; %restore leastCostMatrix

            index = find(loop); %gets index of element > 0 in a colimn vector
            indexTempMinNewCost = find(newCost == tempMinNewCost);
            while (true) %rotate index until index is at top
                if index(1) == Ind
                    break;
                else
                    index = circshift(index,1); %rotate matrix
                end
            end

            s = [row-1, column-1]; %s= size of matrix
            [I,J] = ind2sub(s,index); %generates path from 1's left in an array 'index'
            circle = [I,J]; 

            %get minimum negative element
            key = zeros(length(circle)/2,1);
            columnSelector = 0;
            if pMaxRow < pMaxCol
                columnSelector = 1;
            elseif pMaxRow > pMaxCol || pMaxRow == pMaxCol
                columnSelector = 2;
            end
            switch (columnSelector)
                case 1
                    b = 0; %used to swap between column 1 and 2
                    for i = 1:length(circle)-1 % number of rows-1
                        for j = i+1:length(circle) %next row with reference to the first row
                            if circle(i,b+1) == circle(j,b+1)
                                %do a swap
                                confirmedCircle = circle(i+1,:);
                                circle(i+1,:) = circle(j,:);
                                circle(j,:) = confirmedCircle;
                                b = mod(b+1,2); %swaps odd or even to form the loop
                            end
                        end
                    end
                case 2
                b = 1; %used to swap between column 1 and 2
                    for i = 1:length(circle)-1 % number of rows-1
                        for j = i+1:length(circle) %next row with reference to the first row
                            if circle(i,b+1) == circle(j,b+1)
                                %do a swap
                                confirmedCircle = circle(i+1,:);
                                circle(i+1,:) = circle(j,:);
                                circle(j,:) = confirmedCircle;
                                b = mod(b+1,2); %swaps odd or even to form the loop
                            end
                        end
                    end
            end
            for i= 1:length(circle)
                if mod(i,2)== 0 && leastCostMatrix(circle(i,1), circle(i,2)) > -1%if even, sub
                    key(i/2)= leastCostMatrix(circle(i,1), circle(i,2)); %store the negative element
                end
            end
            minKey = min(key(key > 0)); %minimum negative element

            %input update the costmatrix by adding and substration as need
            for i= 1:length(circle)
                if mod(i,2)== 0 %if even , sub - 1%if odd, add
                    leastCostMatrix(circle(i,1), circle(i,2)) = leastCostMatrix(circle(i,1), circle(i,2))- minKey;
                    if leastCostMatrix(circle(i,1), circle(i,2)) == 0
                        leastCostMatrix(circle(i,1), circle(i,2)) = leastCostMatrix(circle(i,1), circle(i,2))- 1;
                    end
                elseif mod(i,2)== 1 %if even , add
                    leastCostMatrix(circle(i,1), circle(i,2)) = leastCostMatrix(circle(i,1), circle(i,2))+ minKey;
                    if leastCostMatrix(circle(i,1), circle(i,2)) < minKey
                        leastCostMatrix(circle(i,1), circle(i,2)) = leastCostMatrix(circle(i,1), circle(i,2))+ 1;
                    end
                end
            end
            LleastCostMatrix = leastCostMatrix.*(leastCostMatrix>-1);
            iteration = iteration + 1;
            disp('The Temp optimized solution looks like this: ');
            disp(LleastCostMatrix);
            
            %Calculating Temp optimized Cost
            transportationcost = 0;
            for i=1:row-1
                for j=1:column-1
                    transportationcost =transportationcost +(transportationMatrix(i,j)*LleastCostMatrix(i,j));
                end
            end
            disp(['The Temp optimized transportation cost is: ',num2str(transportationcost)]);
        else
            fprintf('Optimal solution reached after %d iteration\n', iteration);
            break
        end
    else
        disp('Northwest solution FAILED! the condition M + N - 1');
        disp('Adding Dummy input to recalculate');
        if ifailed == 0
            findDummy = leastCostMatrix > -1; 
            %Sum this 1 and 0 row wise and column wise
            findDummyRow = sum(findDummy,2);
            findDummyColumn = sum(findDummy,1);
            
            % 2 was used because it is almost always 1 for first and last row/column
            for i = 2:row-2 
                for j=2:column-2
                    if findDummyRow(i) == 1 && findDummyColumn(j)==1
                        leastCostMatrix(i,j) = 0;
                        dummyRow = i;
                        dummyColumn = j;
                    end
                end
            end

            ifailed = ifailed +1;

            LleastCostMatrix = zeros(row-1,column-1);
            for i=1:row-1
                for j=1:column-1
                    if leastCostMatrix(i,j)>-1
                        LleastCostMatrix(i,j) = leastCostMatrix(i,j);
                    end
                end
            end
        else
            leastCostMatrix(dummyRow,dummyColumn) = 0;
        end
        %display updated cost matrix
        if dummyRow ~= 0 && dummyColumn ~= 0
            LleastCostMatrix(dummyRow,dummyColumn) = 0;
        end
        disp('Updated Cost Matrix');
        disp(LleastCostMatrix);
        disp('---------------------------------------------------------------');
    end
    
end
disp('-------------------------------------------------------------------');
%% Calculate Optimized Cost
transportationcost = 0;
for i=1:row-1
    for j=1:column-1
        transportationcost =transportationcost +(transportationMatrix(i,j)*LleastCostMatrix(i,j));
    end
end
disp('Final Cost Matrix');
disp(LleastCostMatrix);
disp(['The transportation cost is: ',num2str(transportationcost)]);
