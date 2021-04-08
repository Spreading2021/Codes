clear
clc
datx = csvread("ITS.csv"); %Data with ~10000 generated paths
%datx = csvread("DataO.csv"); %Data original
n=80; %pack
iter=100000;

m=1;


for l=1:iter
disp(l)
for i=2:n+1
    un_val = 2;
    %un_val = 1;
    up_val = 9269; %degradation paths
    %up_val = 14;  %degradation paths
    r=randi([un_val up_val]);
    for j=2:11
        pack(j,i-1)=datx(j,r);
    end
    a{l}=datx(12,r);
    b{l}=datx(13,r);
    c{l}=datx(14,r);
end
promp = min(pack(11,:));
for i=1:n
    if pack(11,i)==promp
        col_ind=i;
    end
end 
      
promp_a(:,l)=pack(:,col_ind);


tim=[0 17 42 62 102 122 202 442 722 782 1102];
cf{l}=promp_a(:,l)';

packsTF_10{l}=spline(cf{l},tim,10);
packsTF_20{l}=spline(cf{l},tim,20);
packsTF_30{l}=spline(cf{l},tim,30);
packsTF_40{l}=spline(cf{l},tim,40);
packsTF_50{l}=spline(cf{l},tim,50);
packsTF_60{l}=spline(cf{l},tim,60);

TF10{l}=((10/(a{l}*b{l}))+(1/(2*a{l})))^(1/c{l});
TF20{l}=((20/(a{l}*b{l}))+(1/(2*a{l})))^(1/c{l});
TF30{l}=((30/(a{l}*b{l}))+(1/(2*a{l})))^(1/c{l});
TF40{l}=((40/(a{l}*b{l}))+(1/(2*a{l})))^(1/c{l});
TF50{l}=((50/(a{l}*b{l}))+(1/(2*a{l})))^(1/c{l});
TF60{l}=((60/(a{l}*b{l}))+(1/(2*a{l})))^(1/c{l});

err10{l}=(abs(packsTF_10{l}-TF10{l})/packsTF_10{l})*100;
err20{l}=(abs(packsTF_20{l}-TF20{l})/packsTF_20{l})*100;
err30{l}=(abs(packsTF_30{l}-TF30{l})/packsTF_30{l})*100;
err40{l}=(abs(packsTF_40{l}-TF40{l})/packsTF_40{l})*100;
err50{l}=(abs(packsTF_50{l}-TF50{l})/packsTF_50{l})*100;
err60{l}=(abs(packsTF_60{l}-TF60{l})/packsTF_60{l})*100;


end

err10f=mean(cell2mat(err10));

err20f=mean(cell2mat(err20));
err30f=mean(cell2mat(err30));
err40f=mean(cell2mat(err40));
err50f=mean(cell2mat(err50));
err60f=mean(cell2mat(err60));

errorf=[err10f err20f err30f err40f err50f err60f];

