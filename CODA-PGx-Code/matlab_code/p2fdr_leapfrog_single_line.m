function fdr=p2fdr_leapfrog_single_line(pp)

param=[0,1.0000,1.6958,-0.9612];


fdr=param(1)+(param(2)-param(1))./(1+10.^((param(3)-pp)*param(4)));
