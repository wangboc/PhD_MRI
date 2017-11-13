function KY=ReducedSample(k_space_full,R,Km,ACSL,DL,DH,DM)

% ************************************************************
%  Matrix
% ************************************************************
AL=ACSL;
%k_space_red=k_space_full(1:R:end,:,:);
k_space_red(:,:,:)=k_space_full(1:R:Km,:,:);

[dk1,dk2,CoilNum] = size(k_space_red);
if AL==0
   k_space=k_space_red;
else
   k_space(1:floor((DL-2)/R)+1,:,:) = k_space_red(1:floor((DL-2)/R)+1,:,:);
   k_space(floor((DL-2)/R)+2:floor((DL-2)/R)+AL+1,:,:) = k_space_full(DL:DH,:,:);
   k_space(floor((DL-2)/R)+AL+2:dk1+ACSL-DM,:,:) = k_space_red(floor((DL-2)/R)+2+DM:dk1,:,:);
end
   %size(k_space)
%k_space=k_space_red;
for s = 1 : CoilNum
    for kv = 1 : dk2
        dv = (s-1)*(dk1+ACSL-DM)*dk2+(kv-1)*(dk1+ACSL-DM);
        KY(dv+1:dv+dk1+ACSL-DM,1) = k_space(:,kv,s);
    end
end    


