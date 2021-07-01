function [pdf, Coords] = compute_pdf(Data,N)
% based on compute_pdf_fixedbins.m code
% Compute pdf 1D, 2D, or 3D

%N = number of bins 


nTup = size(Data,1);
dim = size(Data,2);

if dim==1
    pdf = zeros(1,N);
elseif dim==2
    pdf=zeros(N,N);
elseif dim==3
    pdf=zeros(N,N,N);
end

Coords = zeros(dim,N);
Edges = zeros(dim,N+1);

%%%%%%%%%%%%%%%%%%%%% Determine bin coordinates %%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:dim
    %[PDF(i,:),Edges(i,:)] = histcounts(Data(:,i),N,'Normalization','probability'); %% we need just Edges to be match with our
    % existing thresholds in quantization.
        
    Edges(i,:)= linspace(min(Data(:,i)-10^-8), max(Data(:,i)),N+1);    
    Coords(i,:)=(Edges(i,1:end-1)+Edges(i,2:end))./2;
    %Coords(i,1)=0;
        
end

         
    for i = 1:dim        
        dat = Data(:,i);
        bindata = ones(size(dat));
        edges = Edges(i,:);
        
        for e = 1:N
            bindata(dat>=edges(e) & dat<edges(e+1))= e;
            if e==N
            bindata(dat>=edges(e+1))=e;  
            end
        end

        BinData(:,i)=bindata; %(i,j,k) bin numbers for each data point
    end
    
 
    if dim==1
        C=zeros(1,N);
        for n = 1:nTup
            dat = BinData(n);
            C(dat)=C(dat)+1;
        end
    elseif dim==2
        C=zeros(N,N);
        for n = 1:nTup
            dat = BinData(n,:);
            C(dat(1),dat(2))=C(dat(1),dat(2))+1;
        end
    elseif dim==3
        C=zeros(N,N,N);
        for n = 1:nTup
            dat = BinData(n,:);
            C(dat(1),dat(2),dat(3))=C(dat(1),dat(2),dat(3))+1;
        end
    end
 
pdf = C./nTup;  
    
end



