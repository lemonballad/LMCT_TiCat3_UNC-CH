function write_cube(header,filename,rho)
[N3,N2,N1]=size(rho);
fN3=floor(N3/6);
f_out=fopen(filename,'w+');
fprintf(f_out,'%s',header);
for n1=1:N1
    for n2=1:N2
        fprintf(f_out,'%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n',rho(1:fN3,n2,n1));
        fprintf(f_out,'%13.5E%13.5E%13.5E%13.5E%13.5E\n',rho(fN3+1:N3,n2,n1));
    end
end
fclose('all');
end