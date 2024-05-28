function fdr=p2fdr_leapfrog_mult_lines(pp,nn)

%nn can be from 1 to 4

if nn==1
    param=[0,1.0000,1.5789,-1.4931];
elseif nn==2
    param=[0,1.0000,2.3137,-0.7396];
elseif nn==3
    param=[0,1.0000,3.4675,-1.7506];
elseif nn==4
    param=[0,1.0000,4.4613,-1.7577];
end



fdr=param(1)+(param(2)-param(1))./(1+10.^((param(3)-pp)*param(4)));
