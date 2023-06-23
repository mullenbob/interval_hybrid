%program to generate truss data for N by M cross braced truss

%N1 and N2 are the number of nodes in x and y direction, respectively 
clear
outfile=fopen('xbrace.txt','w');
n1 =2   % number of nodes in x direction    
n2= 2 % number of nodes in y direction    
nnode=n1*n2   % number of nodes
nele=(n1-1)*(n2)+(n2-1)*n1+(n1-1)*(n2-1)*2  % number of elements total
betax=zeros(nnode,1);
betay=zeros(nnode,1);
alpha1=.03;
alpha2=.03;
alpha3=.03;
alpha4=.03;

ND=2;
group=zeros(1,nele);
alpha=ones(1,nele)*0.005;
force=zeros(ND,nnode);
IFIX=zeros(ND,nnode,1,'int8');
ix=zeros(2,nele,1,'int32');
xleng=1.
yleng=.75;
E=2.00E+08	
area=.001
fload=10.
xc=zeros(nnode,1);
yc=zeros(nnode,1);
delx=xleng/(n1-1);
dely=yleng/(n2-1);

% generate nodal coordiantes
i3=0;
yy=0.;

for i2=1:n2
xx=0.;
for i1=1:n1
    i3=i3+1;
    xc(i3)=xx;
    yc(i3)=yy;
    xx=xx+delx;
    if (i2 == 1)
        %fix bottom nodes in y direction
        IFIX(1,i3)=1;
        IFIX(2,i3)=0;
    else
        IFIX(1,i3)=1;
        IFIX(2,i3)=1;

    end
    if (i1 == 1)
        %put loads on right hand side in x direction of 10. kip
        if (i2 ~= 1)
            force(1,i3)=fload;
            betax(i3)=.05;
        else
            force(1,i3) = 0;
    

    end
        
    end
end
yy=yy+dely;
end

J=0;
%generate element nodal data
i3=0;
%start with horizontal elements
node1=1;
node2=node1+1;

for i2=1:(n2)
    for i1=1:(n1-1)
    i3=i3+1;
    ix(1,i3)=node1;
    ix(2,i3)=node2;
    group(i3)=1;
    alpha(i3)=alpha1;
    node1=node1+1;
    node2=node2+1;
    end
    node1=(i2)*n1+1;
    node2=node1+1;
end

%now do vertical  elements
node1=1;
node2=node1+n1;

for i1=1:(n1)
    for i2=1:(n2-1)
    i3=i3+1;
    ix(1,i3)=node1;
    ix(2,i3)=node2;
    group(i3)=2;
    alpha(i3)=alpha2;
    node1=node2;
    node2=node1+n1;
    end
    node1=(i1)+1;
    node2=node1+n1;
end   
% now do bottom left to upper right  elements
node1=1;
node2=node1+n1+1;
for i1=1:(n1-1)
    for i2=1:(n2-1)
        i3=i3+1;
   ix(1,i3)=node1;
    ix(2,i3)=node2;
    group(i3)=3;
    alpha(i3)=alpha3;
    node1=node1+n1;
    node2=node2+n1;
    end
    node1=i1+1;
    node2=node1+n1;
end
% now do bottom left to upper right  elements
node1=2;
node2=n1+1;
for i1=1:(n1-1)
    for i2=1:(n2-1)
        i3=i3+1;
   ix(1,i3)=node1;
    ix(2,i3)=node2;
    group(i3)=4;
    alpha(i3)=alpha4;
    node1=node1+n1;
    node2=node2+n1;
    end
    node1=i1+1;
    node2=node1+n1;
end
%correct first node to be pinned
       IFIX(1,1)=0;
fprintf(outfile,'%d   %d \n',nnode,nele);
for i1=1:nnode
    fprintf(outfile,'%d  %g %g %d %d %g %g %g %g \n',i1,xc(i1),yc(i1),IFIX(1,i1),IFIX(2,i1),force(1,i1),force(2,i1), betax(i1),betay(i1));
end
for i1=1:nele
    fprintf(outfile,'%d %d %d %g %g %g %g \n',i1,ix(1,i1),ix(2,i1),area,E,alpha(i1),group(i1));
end
fprintf(outfile,'4 \n 1 .01 \n 2 .01 \n 3 .01 \n 4 .01 \n ');

      