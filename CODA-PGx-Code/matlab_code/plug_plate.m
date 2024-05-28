function plug_plate(x,drugs,fig_title)

y=load('in.txt');

if size(y)~=[7,12]
    disp('this dont like right');
    return
end

xx=log10(x);
jj=1;
for n=1:3:10
    
    subplot(2,2,jj);
    yy=y(:,n:n+2);
    fixed_params=[min(mean(yy,2)),max(mean(yy,2)),NaN,nan];
    param=sigm_fit(xx,mean(yy,2),fixed_params);
    
    fprintf(['IC5' '\t' num2str(get_ic(5+min(mean(yy,2)),param)) '\n']);
    fprintf(['IC10' '\t' num2str(get_ic(10+min(mean(yy,2)),param)) '\n']);
    fprintf(['IC15' '\t' num2str(get_ic(15+min(mean(yy,2)),param)) '\n']);
    fprintf(['IC20' '\t' num2str(get_ic(20+min(mean(yy,2)),param)) '\n']);
    fprintf(['IC50' '\t' num2str(get_ic(50+min(mean(yy,2)),param)) '\n']);
    
    fprintf('\n');
    title(drugs(jj))
    jj=jj+1;
end

eval(['print -dpdf ' fig_title]);


function my_ic=get_ic(ic,param)

yval=(ic/100)*(param(2)-param(1))+param(1);
hh=((param(2)-param(1))/(yval-param(1))-1);
xval=-(log10(hh)/param(4))+param(3);

my_ic=10^xval;



    
    
    
    
