function Mf = fillmiss(M,val)

% FILLMISS Interpolation values in a matrix.
%	MF = FILLMISS(M) Interpolates missing vales
%	of matrix M (marked by NaN) from a set of nearest
%	available elements.
%	MF = FILLMISS(M) Returns matrix MF equal to input M
%	except for the elements of M marked by NaN, where
%	MF is interpolated from available neighbouring elements.
%
%  Missing values are calculated as weighted sum of linear
%  interpolations from nearest available points.
%  Altogether 5 estimates from column-wise and 5 for row-wise
%  1-d linear interpolation are calculated.
%  Weights are such that for the best case (isolated missing
%  points away from the boundary) the interpolation is equvalent
%  to average of 4-point Lagrangian polynomial interpolations
%  from nearest points in a row and a column.

%  Kirill K. Pankratov,  kirill@plume.mit.edu
%  12/19/94

% Handle input .............................................

if nargin==0, help fillmiss, return, end

if nargin==1
    % If no missing values, output=input, quick exit.
    % If all values are missing, give up  (output=input)
    % and prevent endless recursion ....................
    if all(all(isnan(M)))|~any(any(isnan(M)))
        Mf = M;
        return
    end
    
    % Interpolation coefficients .......................
    % [left right left-left left-right right-right]
    coef =  [eps eps 1/2 4/3 1/2];
    
    % Sizes ...................
    n_miss = length(find(isnan(M)));
    szM = size(M);
    szMf2 = szM(2)+4;
    is_vec = szM==1;
    
    % Auxillary ...............
    v4 = -1:2;
    o4 = ones(1,4);
    o5 = ones(5,1);
    o2 = ones(2,1);
    omiss = ones(n_miss,1);
    wsum = zeros(n_miss,2);
    a = wsum;
    
    % Interpolate from both rows and columns, if possible
    for jj = find(~is_vec)
        
        bM = zeros(2,szM(3-jj));   % Make margins
        Mf = M; if jj==2, Mf = M'; end
        Mf = [bM; isnan(Mf); bM];
        
        Mf = Mf(:);
        miss = find(Mf==1);       % Missing ##
        exis = find(~Mf);         % Available ##
        Mf = cumsum(~Mf);
        i_m = Mf(miss);
        
        Mf = M; if jj==2, Mf = M'; end
        Mf = [bM*nan; Mf; bM*nan];
        
        % Indices ................
        I = i_m(:,o4)+v4(omiss,:);
        % Quartets of neib. pts for each missing pts.
        I = reshape(exis(I),n_miss,4); 
        
        % Make 5 estimates .........................
        W = miss(:,o4)-I;
        A = zeros(n_miss,5);
        A(:,1:2) = reshape(Mf(I(:,2:3)),n_miss,2);
        A(:,3:5) = (W(:,1:3)-W(:,2:4));
        A(:,3:5) = reshape(Mf(I(:,2:4)),n_miss,3).*W(:,1:3);
        A(:,3:5) = A(:,3:5)-reshape(Mf(I(:,1:3)),n_miss,3).*W(:,2:4);
        A(:,3:5) = A(:,3:5)./(W(:,1:3)-W(:,2:4));
        
        % Calculate weights ......
        W = [abs(W(:,2:3)) abs(W(:,1:3))+abs(W(:,2:4))];
        W = (~isnan(A))./W;
        
        W = W.*coef(omiss,:);
        wsum(:,jj) = sum(W')';
        wsum(:,jj) = wsum(:,jj)-(wsum(:,jj)==0);
        W = W./wsum(:,jj*ones(1,5));
        i_m = find(isnan(A));
        A(i_m) = zeros(size(i_m));
        A = A.*W;
        a(:,jj) = A*o5;
    end
    
    % Correspondence between row and column numbering .....
    exis = ceil(miss/szMf2);
    exis = exis+(miss-(exis-1)*szMf2)*szMf2;
    [exis,i_m] = sort(exis);
    wsum(:,1) = wsum(i_m,1);
    a(i_m,1) = a(:,1);
    
    % Combine estimates from rows and columns .............
    wsum = wsum+(wsum==-1);
    exis = wsum*o2;
    i_m = exis==0;
    exis(i_m) = exis(i_m)+nan;
    wsum = wsum./exis(:,o2);
    exis = (a.*wsum)*o2;
    
    % Insert interpolated values into Mf ..................
    Mf(miss) = exis;
    
    % Remove NaNs at the margins
    Mf = Mf(3:szM(jj)+2,:);
    if jj==2, Mf = Mf'; end
    
    % If there are still missing pts, repeat the procedure
    if any(any(isnan(Mf))), Mf = fillmiss(Mf); end
    
elseif nargin==2 & strcmp(val,'lin')
    % here we do a linear case
    Mf=M;
    in=isnan(Mf);
    [sy,sx]=size(Mf);
    if sy>2 & sx>2 % We have a matrix
        [jj,ii]=find(in==1);
        while ~isempty(jj)
            Mf(1,1)=mnanmean([Mf(1,2) Mf(2,2) Mf(2,1)]);
            Mf(1,end)=mnanmean([Mf(1,end-1) Mf(2,end-1) Mf(2,end)]);
            Mf(end,1)=mnanmean([Mf(end,2) Mf(end-1,2) Mf(end-1,1)]);
            Mf(end,end)=mnanmean([Mf(end,end-1) Mf(end-1,end-1) Mf(end-1,end)]);
            
            for i=2:length(ii)-1
                Mf(1,ii(i))=mnanmean([ Mf(1,ii(i)-1) ...
                        Mf(1,ii(i)+1) Mf(1,ii(i)-1) ...
                        Mf(2,ii(i)) Mf(2,ii(i)+1)]);
                Mf(end,ii(i))=mnanmean([ Mf(end,ii(i)-1) ...
                        Mf(end,ii(i)+1) Mf(end,ii(i)-1) ...
                        Mf(end-1,ii(i)) Mf(end-1,ii(i)+1)]);
            end
            for j=2:length(jj)-1
                Mf(jj(j),1)=mnanmean([Mf(jj(j)-1,1) Mf(jj(j)-1,2) ...
                        Mf(jj(j),2) Mf(jj(j)+1,1) Mf(jj(j)+1,2)]);	
                Mf(jj(j),end)=mnanmean([Mf(jj(j)-1,end) Mf(jj(j)-1,end-1) ...
                        Mf(jj(j),end-1) Mf(jj(j)+1,end) Mf(jj(j)+1,end-1)]);
            end
            
            in(:,1)=0;in(:,end)=0;in(1,:)=0;in(end,:)=0;
            Mf(in)=0;
            
            in=isnan(Mf);
            [jj,ii]=find(in==1);
        end
        
    else % we have a vector
        ii=find(in==1);
        for i=1:length(ii)
            %fac1=
            %fac2=ii(i)+2;
            
            if ii(i)==1
                fac1=ii(i); fac2=ii(i)+2;
            elseif ii(i)==2
                fac1=ii(i)-1; fac2=ii(i)+2;
            elseif (ii(i)+2)==length(Mf)+1
                fac1=ii(i)-2; fac2=ii(i)+1;
            elseif (ii(i)+2)==length(Mf)+2
                fac1=ii(i)-2; fac2=ii(i);
            else
                fac1=ii(i)-2; fac2=ii(i)+2;
            end
            
            Mf(ii(i))=mnanmean(Mf([fac1:fac2])); 
            %elseif ii(i)==length(Mf)
            %	Mf(ii(i))=nanmean([Mf(ii(i)-1) Mf(ii(i)-2)]); 
            %      else
            %	Mf(ii(i))=nanmean([Mf(ii(i)-1) Mf(ii(i)+1)]); 
            %end
        end
    end
    
    

end

  