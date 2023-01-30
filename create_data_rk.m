% RK from Man Institute
% 19 countries
mD = readmatrix("data_rk.csv");
mD(1:5,:)
% mD1 data with Straigts Times upto 2015/9/18
% mD2 data without ST: 2017/12/04 (2017/12/05 omitted)
mD2 = [mD(4:end-1,1:end-2) mD(4:end-1,end)];
mD2 = fillmissing(mD2,'previous');
mD1 = mD(4:4109,:);
mD1 = fillmissing(mD1,'previous');
for t=1:length(mD2)
    if (mD2(t,6)==0)
        mD2(t,6)=mD2(t-1,6);
    end
end
for t=1:length(mD1)
    if (mD1(t,6)==0)
        mD1(t,6)=mD1(t-1,6);
    end
end

iF=100;
sum(sum(isnan(mD1(end-2000-iF+1:end,:))))
sum(sum(isnan(mD2(end-2000-iF+1:end,:))))
[mD1(end-2000-iF+1,1) mD1(end,1)]
[mD2(end-2000-iF+1,1) mD2(end,1)]
save("data_rk.mat","mD1","mD2");
clear
clc