clear
clearvars -global
global inp;
global out;
global tx;
global ty;
global force;
global fxA;
global fyA;
global iii;
global nnd
global nel
global ndof
global telement
global tE;
global tA;
global resxA;
global resyA;
global delta;
global iii;
iii=0;
inp= fopen('Hierarchical100.inp','r');
out =fopen('Hierarchical100optga.out','w'); 
% inp= fopen('popova5.inp','r');
% out =fopen('popova5optga.out','w'); 
readtrussx(inp,out);
x0=tE/1.0e+6;  %set initial values for modulus values
lb=x0*(1-.005);  % set bounds using .005
ub=x0*(1+.005);
tic
A=[];
Aeq=[];
b=[];
beq=[];
disp=trusssolve(x0)
tic
options = optimoptions('fmincon','Algorithm','interior-point','Display','iter','MaxFunctionEvaluations',30000000,'ConstraintTolerance',1e-6,'FunctionTolerance',1.e-10);
%options = optimoptions('fmincon','PlotFcn','optimplotconstrviolation','Algorithm','sqp','Display','iter','MaxFunctionEvaluations',30000000);
%[x,fval,exitflag,output] = fmincon(@funx,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);

[x,fval,exitflag,output] = fmincon(@funx,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);
fprintf(out,' fmincon solution min %15.10e max %15.10e  fvalue min %15.10e  time %s function evals %d\n', min(x),max(x),fval,toc,iii);
iii=0;
tic;
[x,fval,exitflag,output] = fmincon(@funy,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);

fprintf(out,' fmincon solution min %15.10e max %15.10e  fvalue max %15.10e  time %s function evals %d\n', min(x),max(x),fval,toc,iii);
fprintf(out,'\nNodal information - optimization model\n');
fprintf(out,'Node   X      Y     Restraints       Fx          Fy           U-x    U-y\n');
ii=1;

for i=1:nnd
    
     fprintf(out,'%2d   %4.1f   %4.1f    %d    %d     %9.1f    %9.1f           %10.8g      %10.8g\n',i,tx(i),ty(i),resxA(i),resyA(i),fxA(i),fyA(i),delta(ii),delta(ii+1)); 
ii=ii+2;
end


    function [c,ceq] = nonlcon(x);
c = [];
ceq = [ ];
end
function disp=funx(optin);
global iii;
iii=iii+1;
% if (mod(iii,1000)== 0)
%    
% min(optin)
% max(optin)
% end
disp=trusssolve(optin);
return
end
function disp=funy(optin);
global iii;
iii=iii+1;

disp=-trusssolve(optin);
return
end
function readtrussx(inp,out)
   
% set up global variables to exchange inifor between functions
global inp;
global out;
global tx;
global ty;
global force;
global nnd
global nel
global ndof
global telement
global tE;
global tA;
global fxA;
global fyA;
global resxA;
global resyA;



% read in data in same forame as 
                                                                            
[mat1] = fscanf(inp,'%d',2);
nnd= mat1(1) %number of nodes - assembled model
nel  = mat1(2) %number of elements
ndof=2*nnd

%x and y coodinates of nodes for assembled model
tx = zeros(nnd,1);
ty = zeros(nnd,1);
resxA=zeros(nnd,1);
resyA=zeros(nnd,1);
fxA=zeros(nnd,1);
fyA=zeros(nnd,1);

%reading nodal information for the assembled model
[mat1] = fscanf(inp,'%d %f %f %d %d %e %e  %f %f',[9,nnd]);
mat1 = mat1';
tx = mat1(:,2);
ty = mat1(:,3);
%nodal information for assembled model
resxA = mat1(:,4); %x restraint for assembled model
resyA = mat1(:,5); %x restraint for assembled model
%x and y components of forces for assembled model
fxA = mat1(:,6);
fyA = mat1(:,7);

%reading element information for the assembled model
tA       = zeros(nel,1);
tE       = zeros(nel,1);
telement = zeros(nel,2);
[mat2] = fscanf(inp,'%d %d %d %f %e %f %d',[7,nel]);
mat2 = mat2';
telement = mat2(:,2:3); %nodal connectivity matrix for assembled model
tA       = mat2(:,4);
tE       = mat2(:,5);
%alfaA    = mat2(:,6); %uncertainty of E in assembled mode
return
end
function temp=trusssolve(tEnew)
global ndof;
global delta;
global nel;
global tx;
global ty;
global force;
global nnd
global nel
global ndof
global telement
global tE;
global tA;
global fxA;
global fyA;
global resxA;
global resyA;
K=zeros(ndof);
  %assenble stiffness matrix
  for e=1:nel
     i=telement(e,1);
     j=telement(e,2);
     txx=tx(j)-tx(i);
     tyy=ty(j)-ty(i);
     l=sqrt(txx^2+tyy^2);
     c=txx/l;
     s=tyy/l;
     k11=(tEnew(e)*1.0e6*tA(e)/l)*[(c)^2 s*c;s*c s^2];
    K((2*i-1):2*i,(2*i-1):2*i)=K((2*i-1):2*i,(2*i-1):2*i)+k11;
    K((2*j-1):2*j,(2*j-1):2*j)=K((2*j-1):2*j,(2*j-1):2*j)+k11;
    K((2*i-1):2*i,(2*j-1):2*j)=K((2*i-1):2*i,(2*j-1):2*j)-k11;
    K((2*j-1):2*j,(2*i-1):2*i)=K((2*j-1):2*j,(2*i-1):2*i)-k11;
  end
%
% force fixed dof
for e=1:nnd
    
if (resxA(e) == 0)
        ii=e*2-1;
        for i=1:nnd
            K(ii,i)=0.;
            K(i,ii)=0.;
        end
        K(ii,ii)=1.0;
end
if (resyA(e) == 0)
        ii=e*2;
        for i=1:nnd
            K(ii,i)=0.;
            K(i,ii)=0.;
        end
        K(ii,ii)=1.0;
end
end
% build force vector
ii=0;
for i=1:nnd
ii=ii+1;
force(ii)=fxA(i);
ii=ii+1;
force(ii)=fyA(i);
end

delta=K\force';
temp=delta(2*nnd-1)*10000;
return  
end
