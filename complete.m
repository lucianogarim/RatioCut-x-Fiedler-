Tree_Generator_main                        %Program to generate trees
G=graph(adj);
plot(G);
arestas=G.Edges;                            %Import edges from an object
arestas = splitvars(arestas, 'EndNodes');
arestas = arestas{:,:};
arestas=transpose(arestas);
ss=arestas(1,:);
rr=arestas(2,:);
n=size(transpose(ss),1)+1;
size_sets=zeros(n,2);
for j=1:n
    try
        ss(j)=[];
        rr(j)=[];
        GG=graph(ss,rr);
        bins = conncomp(GG);        %take the connected components
        s=sum(bins>1.5);
        scomp=sum(bins < 1.5);
        size_sets(j,1:2)=[s,scomp]
        for i=1:n
            if bins(i)==1
                x(i)= sqrt((scomp)/(s*n));
            else
                x(i)=-sqrt((s)/(scomp*n));
            end
        end
        x=transpose(x);
        rcut(j)= (x' * laplacian * x);
        x=x';
        ss=arestas(1,:);
        rr=arestas(2,:);
        continue
    end
    % Nothing to do
end

[M,I] = min(rcut);              %Take the min RatioCut, and it's position
   ss(I)=[];
    rr(I)=[];
    GG=graph(ss,rr);
    bins = conncomp(GG);        
    s=sum(bins>1.5);
    scomp=sum(bins < 1.5);
    for i=1:n
        try
            if bins(i)==1
                x(i)= sqrt((scomp)/(s*n));
            else
                x(i)=-sqrt((s)/(scomp*n));
            end
            continue
        end
    end
    x=transpose(x);
    rcut(I)= (x' * laplacian * x);
    
    tiledlayout(2,2)
    nexttile
    p = plot(G);                        %Plot the result
    b=zeros(n,1);
    for i=1:n
        if x(i)>=0
            b(i)=i;
        end
    end
    b = nonzeros(b');
    highlight(p,b,'NodeColor','r')
    layout(p,'force3')
view(3)
title('RatioCut Algorithm');


nexttile


diff=abs(size_sets(1:n-2,1)-size_sets(1:n-2,2));
var=[diff transpose(rcut)];

tt=stackedplot(var);
title('Column 1: Difference between the size of sets. Column2: Values of each Ratio Cut')


[autv,autva]=eigs(laplacian,size(laplacian,1));
U2=autv(:,size(laplacian,1)-1);

size_set=zeros(size(laplacian,1),2);
%Fiedler's algorithm (RatioCut formulation)
for j=1:size(laplacian,1)
    t(j)=U2(j,1);                       %set the threshold according to the vector components
    x = (U2 > t(j));                    %taking the elements higher than the threshold
    x = (2 * x) - 1;                    %adjusting the vector x to values ??of 1 and -1 (indicators)
    s=sum(x>-1);                        %defining the sizes of the sets
    scomp=sum(x < 1);
    size_set(j,1:2)=[s,scomp]
    n=s+scomp;
    for i=1:n
        if x(i)>0
            x(i)= sqrt((scomp)/(s*n));  %building the characteristic vectors
        else
            x(i)=-sqrt((s)/(scomp*n));
        end
    end
    rrcut(j) = (x' * laplacian * x);    %calculating the ratiocut %%% The generated partition is in vector x
end
rrcut=transpose(rrcut);

[M,I] = min(rrcut(rrcut>0));
t(I)=U2(I,1);
x = (U2 > t(I));
x = (2 * x) - 1;
s=sum(x>-1);
scomp=sum(x < 1);
n=s+scomp;
for i=1:n
    if x(i)>0
        x(i)= sqrt((scomp)/(s*n));
    else
        x(i)=-sqrt((s)/(scomp*n));
    end
end
rrcut(I) = (x' * laplacian * x);


nexttile
G=graph(adj);
p = plot(G);                        %Plot the result
b=zeros(n,1);
for i=1:n
    if x(i)>=0
        b(i)=i;
    end
end
b = nonzeros(b');
highlight(p,b,'NodeColor','r')
layout(p,'force3')
view(3)
title('Fiedler Algorithm');

nexttile


diff=abs(size_set(:,1)-size_set(:,2));
var=[diff rrcut];

tt=stackedplot(var);
title('Column 1: Difference between the size of sets. Column2: Values of each Fiedler Cut')


