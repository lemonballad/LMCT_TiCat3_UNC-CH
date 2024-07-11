function [NAtoms,r0,N,e,AtNum,q,r,rho]=read_cube(current_cube)

NAtoms=dlmread(current_cube,'',[2 0 2 0]);

r0=dlmread(current_cube,'',[2 1 2 3]);
N=dlmread(current_cube,'',[3 0 5 0])';
e=dlmread(current_cube,'',[3 1 5 3]);

row_1=6;row_2=row_1+NAtoms-1;

AtNum=dlmread(current_cube,'',[6 0 6+NAtoms-1 0]);
q=dlmread(current_cube,'',[6 1 6+NAtoms-1 1]);
r=dlmread(current_cube,'',[6 2 6+NAtoms-1 4]);

row_1=row_2+1;col_1=0;

temp=dlmread(current_cube,'',row_1,col_1)';

[dim_1,dim_2]=size(temp);
area_temp=dim_1*dim_2;

temp=reshape(temp,1,area_temp);

next6=N(3)-mod(N(3),6)+6;
good_indices=mod(1:area_temp,next6)<=N(3)&mod(1:area_temp,next6)~=0;

good_temp=temp(good_indices);
rho=reshape(good_temp,N(3),N(2),N(1));
end