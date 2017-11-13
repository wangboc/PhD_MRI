function ImgRecon=SENSEArbitrary_new(KspaceDataWeighted,WeightingFunctions,ReduceFactor,InitImg,trajectory,area,mask)
% running time 1.4740e+003 for 3 iterations

%transform the matrix KspaceDataWeighted which has the size nc*nk1 to vector m 
%nk1 is less than nk which is original number of K-space data because KspaceDataWeighted is down-sampled data 
[n1,n2,nc]=size(KspaceDataWeighted);

%calculate the diagnoal matrix I which has the size N^2*N^2
[n1,n2,Coilnum]=size(WeightingFunctions);

weight=sqrt(sum(abs(WeightingFunctions).^2,3));
I=(abs(weight)>eps).*(ones(n1,n2)./(weight+eps))+(abs(weight)>eps).*ones(n1,n2);

%calculate the diagnoal matrix D which has the size (nc*nk1)*(nc*nk1)

%KspaceDataWeighted=FullKspaceData(:,(off+1):ReduceFactor:(Dk-off));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CG algorithm
%RequiredAcc=1/100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%a=I*E(H)*D*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:nc 
    %xspaceimage(:,:,k)=INVNUFFT1(Initialimg,KspaceDataWeighted(k,:),trajectory,D1);
    xspaceimage(:,:,k)=EH(trajectory,KspaceDataWeighted(:,:,k),area);
    xspaceimageWeighted(:,:,k)=xspaceimage(:,:,k).*conj(WeightingFunctions(:,:,k));
end
%load my_data3 %just load xspaceimageWeighted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=I.*sum(xspaceimageWeighted,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RequiredAcc=1e-10;
b=InitImg;
p=a;
r=a;
count=0;
delta=r(:)'*r(:)/(a(:)'*a(:))

while (count < 10 & delta > RequiredAcc)
    
    count=count+1;
        
        q=MatrixAfunction3(p,WeightingFunctions,trajectory,I,area); 
        %call function MatrixAfunction to calculate  I*E(H)*D*E*I
        b=b+(r(:)'*r(:))/(p(:)'*q(:))*p;
        
        tempr=r(:);
        r=r-(r(:)'*r(:))/(p(:)'*q(:))*q;
        p=r+(r(:)'*r(:))/(tempr'*tempr)*p;
        delta=r(:)'*r(:)/(a(:)'*a(:))
        %figure, imshow(abs(rot90(I.*b,1)),[])
        b=b.*mask;
end
        
ImgRecon=I.*b;







