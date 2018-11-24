function [S,SD]=MaxPerXtra(M,N,d)

% Finds maximum assignment of M using Hungarian method

% M(i,j,1)=0; % finiteness
% M(i,j,2)=0; % degree
% M(i,j,3)=0; % constant

% for 'large negitive lambda' so that entries of M are sorted
% lexigraphically with the second component in reverse order.

% OUTLINE

% Subtract biggest element from each row

% Subtract biggest element from each col

% While no assignment

% Check for assignment

% Mark rows/cols

% Subtract Biggest element from allo unmarked components and add to all double
% marked components

% loop

% OUTPUT = optimal assignment S and degree of assignment ValA which is
% equal to the multipilicity of -inf as an eigenvalue of the matrix
% polynomial

for i=1:1:N
    for j=1:1:N
        Zeros(i,j)=0; % initialization of zeors array
    end
    ColZeros(i)=0; % col and row counts of zeros
    RowZeros(i)=1;
end
TotalZeros=0;
for i=1:1:N
    [jj,a,b,c]=BiggestRowEntry(M(i,:,:),N,d); % find biggest entry in the row
    for j=1:1:N % now subtract from row
        if (M(i,j,1)==1) % only do anything if vector element is finite
            M(i,j,2)=M(i,j,2)-b;
            M(i,j,3)=M(i,j,3)-c;
        end
    end
    Zeros(i,jj)=1; % location or row i's zero
    ColZeros(jj)=ColZeros(jj)+1;
    TotalZeros=TotalZeros+1;
end

for i=1:1:N
    if (ColZeros(i)==0)
        [jj,a,b,c]=BiggestColEntry(M(:,i,:),N,d);
        for j=1:1:N
            if (M(j,i,1)==1)
                M(j,i,2)=M(j,i,2)-b;
                M(j,i,3)=M(j,i,3)-c;
            end
        end
        Zeros(jj,i)=1;
        RowZeros(jj)=RowZeros(jj)+1;
        ColZeros(i)=ColZeros(i)+1;
        TotalZeros=TotalZeros+1;
    end
end

% NOW EACH ROW AND COL CONTAINS AT LEAST ONE ZERO WE ATTEMPT A ZERO MARKING
AssignedZeros=0;
while(AssignedZeros<N)
    % try assignment
    [AssignedZeros,Zeros]=Assign(Zeros,RowZeros,ColZeros,TotalZeros,N,d);
    if (AssignedZeros<N)
        % determine marking
        % subtract biggest unmarked entry from unmarked entry and add it to the
        % double marked entries
        [M,Zeros,ColZeros,RowZeros,TotalZeros]=MarkandSwitch(M,Zeros,ColZeros,RowZeros,TotalZeros,N,d);
    end
end

% Finally we have the optiaml assignment

for i=1:1:N
    for j=1:1:N
        if (Zeros(i,j)==3)
            S(i)=j; 
            SD(i)=M(i,j,2); % stores which edge is in the optimal assignment
        end
    end
end

end


function [M,Zeros,ColZeros,RowZeros,TotalZeros]=MarkandSwitch(M,Zeros,ColZeros,RowZeros,TotalZeros,N,d)



% First mark rows with no assigned zeros
% mark cols with zeros in marked rows
% mark row with assigend zeros in marked cols
% ...

% index rows cols with 0, 1-marked, 2-marked and checked

unchecked=N;
for i=1:1:N
    Rflag(i)=1;
    Cflag(i)=0;
    for j=1:1:N
        if (Zeros(i,j)==3)
            Rflag(i)=0;
            unchecked=unchecked-1;
        end
    end
end

while(unchecked>0)
    
    % set to zero if we find a row or col with index 1
    for i=1:1:N
        if (Rflag(i)==1)
            for j=1:1:N
                if (Zeros(i,j)>0.5)
                    if (Cflag(j)==0)
                        Cflag(j)=1;
                        unchecked=unchecked+1;
                    end
                end
            end
            unchecked=unchecked-1;
            Rflag(i)=2;
        end
        if (Cflag(i)==1)
            for j=1:1:N
                if (Zeros(j,i)==3)
                    if (Rflag(j)==0)
                        Rflag(j)=1;
                        unchecked=unchecked+1;
                    end
                end
            end
            unchecked=unchecked-1;
            Cflag(i)=2;
        end
    end
    
end

[ii,jj,b,c]=BiggestUnMarked(M,Rflag,Cflag,N,d); % Find biggest unmarked entry in M
for i=1:1:N
    for j=1:1:N
        if (Zeros(i,j)>0.5) % reset zeros
            Zeros(i,j)=1;
        end
        if (Rflag(i)==2 && Cflag(j)==0) % and subtract it from all unmarked
            if (M(i,j,1)==1) % only bother if finite
                M(i,j,2)=M(i,j,2)-b;
                M(i,j,3)=M(i,j,3)-c;
            end
        end
        if (Rflag(i)==0 && Cflag(j)==2) % and add it to all double marked
            if (M(i,j,1)==1) % only bother if finite
                M(i,j,2)=M(i,j,2)+b;
                M(i,j,3)=M(i,j,3)+c;
            end
            if (Zeros(i,j)==1)
            Zeros(i,j)=0; % loose any double marked zeros
            RowZeros(i)=RowZeros(i)-1;
            ColZeros(j)=ColZeros(j)-1;
            TotalZeros=TotalZeros-1;
            end
        end
    end
end
Zeros(ii,jj)=1; % new zero from unmarked maximium
RowZeros(ii)=RowZeros(ii)+1;
ColZeros(jj)=ColZeros(jj)+1;
TotalZeros=TotalZeros+1;
 
end








function [ii,jj,b,c]=BiggestUnMarked(M,Rflag,Cflag,N,d)


% returns value and position of the biggest unmarked entry in M

seenany=0;
for i=1:1:N
    for j=1:1:N
        if (Rflag(i)==2 && Cflag(j)==0)
            if (seenany==0)
                seenany=1;
                a=M(i,j,1); % have to choose the first on we see incase all infinite
                b=M(i,j,2);
                c=M(i,j,3);
                ii=i;
                jj=j;
            else
                if (M(i,j,1)>a) % if j is fininte but current best is infinite then take it
                    a=M(i,j,1);
                    b=M(i,j,2);
                    c=M(i,j,3);
                    ii=i;
                    jj=j;
                elseif (M(i,j,1)==a) % if same finiteness then compare lambda coefficiant
                    if (M(i,j,2)<b) % reverse order of lambda coefficiant since assumption is that lambda is very large and negitive
                        a=M(i,j,1);
                        b=M(i,j,2);
                        c=M(i,j,3);
                        ii=i;
                        jj=j;
                    elseif (M(i,j,2)==b && M(i,j,3)>c) % if same lambda coefficiant then compare constant term
                        a=M(i,j,1);
                        b=M(i,j,2);
                        c=M(i,j,3);
                        ii=i;
                        jj=j;
                    end
                end
            end
        end
    end
end
            

end









function [AssignedZeros,Zeros]=Assign(Zeros,RowZeros,ColZeros,TotalZeros,N,d)



% try to assign N zeors, zeros=2 if marked but unassigned, zeros=3 if
% assigned

AssignedZeros=0;
MarkedZeros=0;

while(MarkedZeros<TotalZeros)
    
   ii=0; % Search for most isolated zero
   jj=0;
   br=N;
   bc=N;
   for i=1:1:N
       if (RowZeros(i)>0 && RowZeros(i)<br)
           br=RowZeros(i);
           ii=i;
       end
       if (ColZeros(i)>0 && ColZeros(i)<bc)
           bc=ColZeros(i);
           jj=i;
       end
   end
   if (br<bc)
       for i=1:1:N
           if (Zeros(ii,i)==1)
               jj=i;
           end
       end
   else
       for i=1:1:N
           if (Zeros(i,jj)==1)
               ii=i;
           end
       end
   end
   % Mark zeros appropriately
   for i=1:1:N 
       if (Zeros(ii,i)==1) % remove row
           Zeros(ii,i)=2;
           RowZeros(ii)=RowZeros(ii)-1;
           ColZeros(i)=ColZeros(i)-1;
           MarkedZeros=MarkedZeros+1;
       end
       if (Zeros(i,jj)==1) % remove col
           Zeros(i,jj)=2;
           RowZeros(i)=RowZeros(i)-1;
           ColZeros(jj)=ColZeros(jj)-1;
           MarkedZeros=MarkedZeros+1;
           
       end
   end
   
   Zeros(ii,jj)=3; % assign ii,jj
   AssignedZeros=AssignedZeros+1;

end

end










function [jj,a,b,c]=BiggestRowEntry(M,N,d)



% retruns the position and value of the biggest element in the vector M(1,:,123)

a=M(1,1,1);
b=M(1,1,2);
c=M(1,1,3);
jj=1;
for j=1:1:N
    if (M(1,j,1)>a) % if j is fininte but current best is infinite then take it
        a=M(1,j,1);
        b=M(1,j,2);
        c=M(1,j,3);
        jj=j;
    elseif (M(1,j,1)==a) % if same finiteness then compare lambda coefficiant
        if (M(1,j,2)<b) % reverse order of lambda coefficiant since assumption is that lambda is very large and negitive
            a=M(1,j,1);
            b=M(1,j,2);
            c=M(1,j,3);
            jj=j;
        elseif (M(1,j,2)==b && M(1,j,3)>c) % if same lambda coefficiant then compare constant term
            a=M(1,j,1);
            b=M(1,j,2);
            c=M(1,j,3);
            jj=j;
        end
    end
end
end










function [jj,a,b,c]=BiggestColEntry(M,N,d)



% retruns the position and value of the biggest element in the vector M(:,1,123)

a=M(1,1,1);
b=M(1,1,2);
c=M(1,1,3);
jj=1;
for j=1:1:N
    if (M(j,1,1)>a) % if j is fininte but current best is infinite then take it
        a=M(j,1,1);
        b=M(j,1,2);
        c=M(j,1,3);
        jj=j;
    elseif (M(j,1,1)==a) % if same finiteness then compare lambda coefficiant
        if (M(j,1,2)<b) % reverse order of lambda coefficiant since assumption is that lambda is very large and negitive
            a=M(j,1,1);
            b=M(j,1,2);
            c=M(j,1,3);
            jj=j;
        elseif (M(j,1,2)==b && M(j,1,3)>c) % if same lambda coefficiant then compare constant term
            a=M(j,1,1);
            b=M(j,1,2);
            c=M(j,1,3);
            jj=j;
        end
    end
end
end

