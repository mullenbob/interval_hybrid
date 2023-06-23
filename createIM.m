function [aa,Do,Cm]=createIM (ne,nnA,ndof,ndof1,strainflag,ndofA,xA,yA,elementA,A,E,resxA,resyA,alfaA,gammao)
% this function is used in the hybrid method to create the interval fixed point matricies that do not change
a = zeros (ne,ndof);
aA = zeros(ne,ndofA);
Le=zeros(ne,1);
cosine=zeros(ne,1);
sine=zeros(ne,1);
Doo=zeros(ne,1); % changeto vector from diaganol 
%allocated space for D arrays
Do=Doo*infsup(1.,1.);
D=Do;
DoD=Do;
Dom=Doo;
 for e=1:ne
  connA=elementA(e,:);
  x1 = xA(connA(1));
  y1 = yA(connA(1));
  x2 = xA(connA(2));
  y2 = yA(connA(2));
  le=sqrt((x2-x1)^2+(y2-y1)^2);
  Le(e) = le;
  c=(x2-x1)/le;
  s=(y2-y1)/le;
  cosine(e) = c;
  sine(e) = s;
 
  dof = [2*e-1 2*e];
 
  dofA=[2*connA(1)-1 2*connA(1) 2*connA(2)-1 2*connA(2)];
  
  a(e,dof)=[1 -1];
 
  aA(e,dofA)= [c s -c -s]' ; %R*[1 0 -1  0]'; % again   why not just vector [c s -c -c]'
  %____________diagonal member matrix_______________________

 %Doo(e)=A(e)*E(e)/le; 
 Do(e)=A(e)*E(e)/le*infsup(1.0,1.0);
 %D(e)=(A(e)*E(e)/le*Alfa(e));
 %DoD(e)=Do(e)-D(e) ; 
 Dom(e) = (E(e)*A(e)/le)*(1.+alfaA(e)*gammao(e));
 %Dom2(e,e) = (E(e)*A(e)/le)*(1.+alfaA(e)*gammao(e));
end
aa=zeros(ne,4*ne+ndof1+strainflag*ndof);
aa(:,1:ndof)=a;
%----------------formulation of the new approach-------------------------


% one constraint term for each EBE node
i2=1;
cons1=zeros(ndof+ndof1,ndof);
for i1=1:ne
    cons1(i2,i2)=1;
    i3=elementA(i1,1)+ne;
    i4=2*i3-1;
    cons1(i4,i2)=-cosine(i1);
    i4=i4+1;
    cons1(i4,i2)=-sine(i1);
    i3=elementA(i1,2)+ne;
    i2=i2+1;
    cons1(i2,i2)=1.;
    i4=2*i3-1;
    cons1(i4,i2)=-cosine(i1);
    i4=i4+1;
    cons1(i4,i2)=-sine(i1);
    i2=i2+1;
end
if (strainflag== 1)
    i2=0;
    cons2=zeros(ndof,ne);
    for i1=1:ne
    i2=i2+1;
    cons2(i2,i1)=-1/Le(i1);
    i2=i2+1;
    cons2(i2,i1)=1./Le(i1);
    end
istrain=ndof+ndof1+ndof+ndof;
end
inostrain=ndof+ndof1+ndof;
%two centered matriese K no group and Km with 1+gamma*delta term
if (strainflag == 1)
    Km=zeros(istrain,istrain);
    Km(1:ndof,inostrain+1:inostrain+ne)=cons2;
    Km(inostrain+1:inostrain+ne,1:ndof)=cons2';
    Km(inostrain+1:inostrain+ne,inostrain+ne+1:inostrain+ndof)=-eye(ne,ne);
    Km(inostrain+ne+1:inostrain+ndof,inostrain+1:inostrain+ne)=-eye(ne,ne);
    % K=Km;
else
    Km=zeros(inostrain,inostrain);
    % K=Km;
end 

    Km(1:ndof,1:ndof)=a'*diag(Dom)*a;
       
    Km(1:ndof+ndof1,ndof+ndof1+1:inostrain)=cons1;
    Km(ndof+ndof1+1:inostrain,1:ndof+ndof1)=cons1';
    % K(1:ndof,1:ndof)=a'*diag(Doo)*a;
    % K(1:ndof+ndof1,ndof+ndof1+1:inostrain)=cons1;
    % K(ndof+ndof1+1:inostrain,1:ndof+ndof1)=cons1';
bn2 = 0;
for i=1:nnA
     if(resxA(i)==0)
           bn2 = bn2+1;
           ifix3(bn2) = ndof+ 2*i-1;
     end
     if(resyA(i)==0)
          bn2 = bn2+1;
          ifix3(bn2) = ndof+ 2*i;
    end
end
Km(ifix3,:)=zeros(bn2,4*ne+ndof1+strainflag*ndof);
Km(:,ifix3)=zeros(4*ne+ndof1+strainflag*ndof,bn2);
Km(ifix3,ifix3)=eye(length(ifix3));
% K(ifix3,:)=zeros(bn2,4*ne+ndof1+strainflag*ndof);
% K(:,ifix3)=zeros(4*ne+ndof1+strainflag*ndof,bn2);
% K(ifix3,ifix3)=eye(length(ifix3));
%check for centered solution
Cm=inv(Km);
return
end