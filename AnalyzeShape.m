%% -----------------------------------------
% Make Figures for Nucleus Volume Paper
%
% Descirbed in: Periodicity of nuclear morphology in human fibroblasts
% Last Update: 9/26/14
% Author: Laura Seaman
% Contact: laseaman@umich.edu
% functions used: findAB.m, SKL.m
%% -----------------------------------------

%% Load Data
'loading data and calculating averages'
load('nucleus_data')
% includes: 'data','data_full','volume','eccentricity','times'
% data/data_full: 3x Num_samples x32x20 eigenvalues, sample, time point, nucleus
% volume/ratio/eccentricity/trce: num_samples x32x20 sample, time, nucleus
% times: 1x32 time (hr) of each time pt

num_samples = 6;
times3 = [times(1:5), times(7:32)];
f = 1; %figure counter
s = 50:50:1000; %locations for figure positioning

%% calculate lots of averages

%average for all samples
avgEV  = squeeze(nanmean(nanmean(data,2),4)); %3x32
EVc = bsxfun( @minus, data,nanmean(nanmean(data,3),4)) + 20; %c-mean centered

%average seperated by sample
avgEV_s  = nanmean(data,4); %3xNSx32
avgEcc_s = nanmean(eccentricity,3); %NSx32
avgVol_s = nanmean(volume,3); %NSx32
% f-full (3A5 filled in), c - mean centered
Vol_f = squeeze((4/3).*pi.*data_full(1,:,:,:).*data_full(2,:,:,:).*data_full(3,:,:,:));
Ecc_f = squeeze(sqrt((data_full(3,:,:,:).^2-data_full(1,:,:,:).^2)./data_full(3,:,:,:).^2));
avgVol_sfc = bsxfun( @rdivide, nanmean(Vol_f,3),nanmean(nanmean(Vol_f,3),2))+mean(nanmean(nanmean(Vol_f,3),2)); %NSx32
avgEcc_sfc = bsxfun( @minus, nanmean(Ecc_f,3),nanmean(nanmean(Ecc_f,3),2))+mean(nanmean(nanmean(Ecc_f,3),2)); %NSx32

stdEV  = nanstd(reshape(permute(data,[1,3,2,4]),3,32,[]),0,3); %3x6x32x20 -> 3x32
stdEV_s  = nanstd(data,0,4); %3x6x32x20 -> 3x6x32
stdEV_sc  = nanstd(EVc,0,4); %3x6x32x20 -> 3x6x32
stdEVc  = nanstd(reshape(permute(EVc,[1,3,2,4]),3,32,[]),0,3); %3x6x32x20 -> 3x32

%% Summary Fig

% make tubes by sample and for average - for average std pool all pre std calc
Nth=100;  % num points around each ellipse
theta=linspace(0,2*pi,Nth);
z_cent = linspace(150,00,6); t_sep = 105;
z_sub = 10; %the amount subtracted from the second set

f_tubes = figure('position',[s(f) 100 500 1000]); f = f+1;

%x-y average of all samples tube
Y1 = (avgEV(3,:)+stdEV(3,:)); Z1 = (avgEV(2,:)+stdEV(2,:));
Y2 = (avgEV(3,:)-stdEV(3,:)); Z2 = (avgEV(2,:)-stdEV(2,:));
y1=Y1'*cos(theta); z1=Z1'*sin(theta)+180;
y2=Y2'*cos(theta); z2=Z2'*sin(theta)+180;
x=repmat(times',1,Nth);
surf(x,y1,z1,'FaceColor',[.1 .5 1],'EdgeAlpha',0,'MeshStyle','row','FaceAlpha',.5);
hold on
surf(x,y2,z2,'FaceColor',[0 .2 1],'EdgeAlpha',.3,'MeshStyle','row','FaceAlpha',1);

%x-z average of all samples tube
Y1 = (avgEV(3,:)+stdEV(3,:)); Z1 = (avgEV(1,:)+stdEV(1,:));
Y2 = (avgEV(3,:)-stdEV(3,:)); Z2 = (avgEV(1,:)-stdEV(1,:));
y1=Y1'*cos(theta); z1=Z1'*sin(theta)+180-z_sub;
y2=Y2'*cos(theta); z2=Z2'*sin(theta)+180-z_sub;
x=repmat(times',1,Nth)+t_sep;
surf(x,y1,z1,'FaceColor',[.1 .5 1],'EdgeAlpha',0,'MeshStyle','row','FaceAlpha',.5);
surf(x,y2,z2,'FaceColor',[0 .2 1],'EdgeAlpha',.3,'MeshStyle','row','FaceAlpha',1);

%axis & scale bars
set(gca,'DataAspectRatio',[2 1 1])
axis off 
axis([-50 180 -15 15 -25 190])
text(-50,7,186,'All samples','FontSize',12) 
plot3(17:32,repmat(-5,1,16) ,repmat(-15,1,16),'Color','k','LineWidth',3)
text(18,-7,-17,'15 hours','FontSize',11) 
plot3(repmat(-17,1,11),zeros(11,1) ,-15:-5,'Color','k','LineWidth',3)
text(-35,6,-10,'10 \mum','FontSize',11)

%% tubes for each sample

for sample = [1,2,4,5,6];
    
    %x-y single sample tube
    Y1 = squeeze((avgEV_s(3,sample,:)+stdEV_s(3,sample,:)))';
    Z1 = squeeze((avgEV_s(2,sample,:)+stdEV_s(2,sample,:)))';
    Y2 = squeeze((avgEV_s(3,sample,:)-stdEV_s(3,sample,:)))';
    Z2 = squeeze((avgEV_s(2,sample,:)-stdEV_s(2,sample,:)))';
    y1=Y1'*cos(theta); z1=Z1'*sin(theta)+z_cent(sample); 
    y2=Y2'*cos(theta); z2=Z2'*sin(theta)+z_cent(sample);
    x=repmat(times',1,Nth);
    surf(x,y1,z1,'FaceColor',[.1 .5 1],'EdgeAlpha',0,'MeshStyle','row','FaceAlpha',.5);
    surf(x,y2,z2,'FaceColor',[0 .2 1],'EdgeAlpha',.3,'MeshStyle','row','FaceAlpha',1);
    text(-50,5,z_cent(sample)+7,{'Human';['Sample ',num2str(sample)]},'FontSize',12)
    
    %x-z single sample tube
    Y1 = squeeze((avgEV_s(3,sample,:)+stdEV_s(3,sample,:)))';
    Z1 = squeeze((avgEV_s(1,sample,:)+stdEV_s(1,sample,:)))';
    Y2 = squeeze((avgEV_s(3,sample,:)-stdEV_s(3,sample,:)))';
    Z2 = squeeze((avgEV_s(1,sample,:)-stdEV_s(1,sample,:)))';
    y1=Y1'*cos(theta); z1=Z1'*sin(theta)+z_cent(sample)-z_sub; 
    y2=Y2'*cos(theta); z2=Z2'*sin(theta)+z_cent(sample)-z_sub;
    x=repmat(times',1,Nth)+t_sep;
    surf(x,y1,z1,'FaceColor',[.1 .5 1],'EdgeAlpha',0,'MeshStyle','row','FaceAlpha',.5);
    surf(x,y2,z2,'FaceColor',[0 .2 1],'EdgeAlpha',.3,'MeshStyle','row','FaceAlpha',1);
end

%sample 3 - seperated because of missing data
sample = 3;
std_s3 = squeeze(cat(3,stdEV_s(:,3,1:5), stdEV_s(:,3,7:32)));
avg_EV_s3 = squeeze(cat(3,avgEV_s(:,3,1:5), avgEV_s(:,3,7:32)));

%x-y sample 3 tube
Y1 = squeeze(avg_EV_s3(3,:)+std_s3(3,:))';
Z1 = squeeze(avg_EV_s3(2,:)+std_s3(2,:))';
Y2 = squeeze(avg_EV_s3(3,:)-std_s3(3,:))';
Z2 = squeeze(avg_EV_s3(2,:)-std_s3(2,:))';
y1=Y1*cos(theta); z1=Z1*sin(theta)+z_cent(sample); 
y2=Y2*cos(theta); z2=Z2*sin(theta)+z_cent(sample);
x=repmat(times3',1,Nth);
surf(x,y1,z1,'FaceColor',[.1 .5 1],'EdgeAlpha',0,'MeshStyle','row','FaceAlpha',.5);
surf(x,y2,z2,'FaceColor',[0 .2 1],'EdgeAlpha',.3,'MeshStyle','row','FaceAlpha',1);
text(-50,5,z_cent(sample)+7,{'Human';['Sample ',num2str(sample)]},'FontSize',12)

%x-z sample 3 tube
Y1 = squeeze(avg_EV_s3(3,:)+std_s3(3,:))';
Z1 = squeeze(avg_EV_s3(1,:)+std_s3(1,:))';
Y2 = squeeze(avg_EV_s3(3,:)-std_s3(3,:))';
Z2 = squeeze(avg_EV_s3(1,:)-std_s3(1,:))';
y1=Y1*cos(theta); z1=Z1*sin(theta)+z_cent(sample)-z_sub; 
y2=Y2*cos(theta); z2=Z2*sin(theta)+z_cent(sample)-z_sub;
x=repmat(times3',1,Nth)+t_sep;
surf(x,y1,z1,'FaceColor',[.1 .5 1],'EdgeAlpha',0,'MeshStyle','row','FaceAlpha',.5);
surf(x,y2,z2,'FaceColor',[0 .2 1],'EdgeAlpha',.3,'MeshStyle','row','FaceAlpha',1);

%% Volume and Eccentricity plots

f_VE = figure('position',[s(f) 550 350 600]); f = f+1;
subplot(211)
plot(times,avgVol_s,'LineWidth',1.5)
xlabel('Time (hr)','FontSize',13), ylabel('Volume (\mum)','FontSize',13)
set(gca,'fontsize',11)
subplot(212)
plot(times, avgEcc_s,'LineWidth',1.5)
xlabel('Time (hr)','FontSize',13), ylabel('Eccentricity','FontSize',13)
set(gca,'fontsize',11)
legend('Sample 1','Sample 2','Sample 3','Sample 4','Sample 5','Sample 6','Location','SouthWest')

%% Calculate best fit

freq_range = 1/35:.002:1/3;

F_ecc = NaN(num_samples,length(freq_range));
alpha_ecc = NaN(num_samples,length(freq_range));
beta_ecc = NaN(num_samples,length(freq_range));
Fs_ecc = NaN(num_samples,length(freq_range));
F_vol = NaN(num_samples,length(freq_range));
alpha_vol = NaN(num_samples,length(freq_range));
beta_vol = NaN(num_samples,length(freq_range));
Fs_vol = NaN(num_samples,length(freq_range));
F_EV = NaN(3,num_samples,length(freq_range));
alpha_EV = NaN(3,num_samples,length(freq_range));
beta_EV = NaN(3,num_samples,length(freq_range));
Fs_EV = NaN(3,num_samples,length(freq_range));

sum_ts = NaN(6,5);

for sample = 1:num_samples;
    fi_ecc = avgEcc_s(sample,:)-nanmean(avgEcc_s(sample,:));
    fi_vol = avgVol_s(sample,:)-nanmean(avgVol_s(sample,:));
    fi_EV1 = reshape(avgEV_s(1,sample,:)-nanmean(avgEV_s(1,sample,:)),1,[]);
    fi_EV2 = reshape(avgEV_s(2,sample,:)-nanmean(avgEV_s(2,sample,:)),1,[]);
    fi_EV3 = reshape(avgEV_s(3,sample,:)-nanmean(avgEV_s(3,sample,:)),1,[]);
    c = 1;
    for w = freq_range
        [alpha_ecc(sample,c),beta_ecc(sample,c),F_ecc(sample,c),Fs_ecc(sample,c)] ...
            = findAB(fi_ecc,w,times);
        [alpha_vol(sample,c),beta_vol(sample,c),F_vol(sample,c),Fs_vol(sample,c)] ...
            = findAB(fi_vol,w,times);
        [alpha_EV(1,sample,c),beta_EV(1,sample,c),F_EV(1,sample,c),Fs_EV(1,sample,c)] ...
            = findAB(fi_EV1,w,times);
        [alpha_EV(2,sample,c),beta_EV(2,sample,c),F_EV(2,sample,c),Fs_EV(2,sample,c)] ...
            = findAB(fi_EV2,w,times);
        [alpha_EV(3,sample,c),beta_EV(3,sample,c),F_EV(3,sample,c),Fs_EV(3,sample,c)] ...
            = findAB(fi_EV3,w,times);
        c =c+1;
    end
    sum_ts(sample,1) = nansum(fi_vol.^2);
    sum_ts(sample,2) = nansum(fi_ecc.^2);
    sum_ts(sample,3) = nansum(fi_EV1.^2);
    sum_ts(sample,4) = nansum(fi_EV2.^2);
    sum_ts(sample,5) = nansum(fi_EV3.^2);
    
end
sum_ts

%% plot spectrums

% Avg spectrum
f_AvgSp = figure('position',[s(f) 100 900 1000]); f = f+1;
subplot(521)
FsMean = nanmean(Fs_vol);
plot(freq_range,FsMean)
title('Volume','FontSize', 13)
ylabel('Power','FontSize', 12)
fInd = find(FsMean== nanmax(FsMean));
hold on, plot(freq_range(fInd),FsMean(fInd),'r.', 'MarkerSize',15)
text(freq_range(fInd)+5*10^-3,FsMean(fInd),['P= ' num2str(1./freq_range(fInd)) ' hr'],'FontSize', 12)
fInd = find(FsMean== nanmax(FsMean(1:10)));
plot(freq_range(fInd),FsMean(fInd),'c.', 'MarkerSize',15)
text(freq_range(fInd)+1.3*10^-2,FsMean(fInd),['P= ' num2str(1./freq_range(fInd)) ' hr'],'FontSize', 12)

subplot(523)
FsMean = nanmean(Fs_ecc);
plot(freq_range,FsMean)
title('Eccentricity','FontSize', 13) 
ylabel('Power','FontSize', 12)
fInd = find(FsMean== nanmax(FsMean));
hold on, plot(freq_range(fInd),FsMean(fInd),'r.', 'MarkerSize',15)
text(freq_range(fInd)+5*10^-3,FsMean(fInd)+2*10^-5,['P= ' num2str(1./freq_range(fInd)) ' hr'],'FontSize', 12)
fInd = find(FsMean== nanmax(FsMean(3:10)));
plot(freq_range(fInd),FsMean(fInd),'c.', 'MarkerSize',15)
text(freq_range(fInd)+1*10^-2,FsMean(fInd)-2*10^-5,['P= ' num2str(1./freq_range(fInd)) ' hr'],'FontSize', 12)

subplot(525)
FsMean = nanmean(squeeze(Fs_EV(3,:,:)),1);
plot(freq_range,FsMean)
title('Long Axis','FontSize', 13)
xlabel('Frequency (1/hr)','FontSize', 12)
ylabel('Power','FontSize', 12)
fInd = find(FsMean== nanmax(FsMean));
hold on, plot(freq_range(fInd),FsMean(fInd),'r.', 'MarkerSize',15)
text(freq_range(fInd)+10^-2,FsMean(fInd)-.005,['P= ' num2str(1./freq_range(fInd)) ' hr'],'FontSize', 12)
fInd = find(FsMean== nanmax(FsMean(1:10)));
hold on, plot(freq_range(fInd),FsMean(fInd),'c.', 'MarkerSize',15)
text(freq_range(fInd)+1.5*10^-2,FsMean(fInd),['P= ' num2str(1./freq_range(fInd)) ' hr'],'FontSize', 12)

subplot(527)
FsMean = nanmean(squeeze(Fs_EV(2,:,:)),1);
plot(freq_range,FsMean)
title('Middle Axis','FontSize', 13)
ylabel('Power','FontSize', 12)
fInd = find(FsMean== nanmax(FsMean));
hold on, plot(freq_range(fInd),FsMean(fInd),'r.', 'MarkerSize',15)
text(freq_range(fInd)+10^-2,FsMean(fInd)-.1,['P= ' num2str(1./freq_range(fInd)) ' hr'],'FontSize', 12)
fInd = find(FsMean== nanmax(FsMean(1:10)));
hold on, plot(freq_range(fInd),FsMean(fInd),'c.', 'MarkerSize',15)
text(freq_range(fInd)+1.2*10^-2,FsMean(fInd),['P= ' num2str(1./freq_range(fInd)) ' hr'],'FontSize', 12)

subplot(529)
FsMean = nanmean(squeeze(Fs_EV(1,:,:)),1);
plot(freq_range,FsMean)
title('Short Axis','FontSize', 13)
ylabel('Power','FontSize', 12), xlabel('Frequency (1/hr)','FontSize', 12)
fInd = find(FsMean== nanmax(FsMean));
hold on, plot(freq_range(fInd),FsMean(fInd),'c.', 'MarkerSize',15)
text(freq_range(fInd)+2*10^-3,FsMean(fInd),['P= ' num2str(1./freq_range(fInd)) ' hr'],'FontSize', 12)
fInd = find(FsMean== nanmax(FsMean(10:50)));
hold on, plot(freq_range(fInd),FsMean(fInd),'r.', 'MarkerSize',15)
text(freq_range(fInd)+2*10^-3,FsMean(fInd),['P= ' num2str(1./freq_range(fInd)) ' hr'],'FontSize', 12)

subplot(522)
plot(freq_range,Fs_vol)
title('Volume','FontSize', 13)
subplot(524)
plot(freq_range,Fs_ecc)
title('Eccentricity','FontSize', 13)
%legend('Sample 1','Sample 2','Sample 3','Sample 4','Sample 5','Sample 6')
subplot(526)
plot(freq_range,squeeze(Fs_EV(3,:,:)))
title('Long Axis','FontSize', 13)
subplot(528)
plot(freq_range,squeeze(Fs_EV(2,:,:)))
title('Middle Axis 2','FontSize', 13)
subplot(5,2,10)
plot(freq_range,squeeze(Fs_EV(1,:,:)))
title('Short Axis','FontSize', 13)
xlabel('Frequency (1/hr)','FontSize', 12)

%% Plot Best Fits

avgVol_sc = 1/6.*nansum(bsxfun(@minus,avgVol_s,nanmean(avgVol_s,2)));
avgEcc_sc = 1/6.*nansum(bsxfun(@minus,avgEcc_s,nanmean(avgEcc_s,2)));
avgEV1_sc  = 1/6.*nansum(bsxfun(@minus,squeeze(avgEV_s(1,:,:)),nanmean(squeeze(avgEV_s(1,:,:)),2)));
avgEV2_sc  = 1/6.*nansum(bsxfun(@minus,squeeze(avgEV_s(2,:,:)),nanmean(squeeze(avgEV_s(2,:,:)),2)));
avgEV3_sc  = 1/6.*nansum(bsxfun(@minus,squeeze(avgEV_s(3,:,:)),nanmean(squeeze(avgEV_s(3,:,:)),2)));

fInd_vol = find(nanmean(Fs_vol)== nanmax(nanmean(Fs_vol)));
freq_vol = freq_range(fInd_vol); 
as_bc_vol = zeros(6,32);
si_vol = sin((2*pi*freq_vol).*times); ci_vol = cos((2*pi*freq_vol).*times);
fInd_ecc = find(nanmean(Fs_ecc)== nanmax(nanmean(Fs_ecc)));
freq_ecc = freq_range(fInd_ecc); 
as_bc_ecc = zeros(6,32);
si_ecc = sin((2*pi*freq_ecc).*times); ci_ecc = cos((2*pi*freq_ecc).*times);
fInd_EV1 = find(nanmean(squeeze(Fs_EV(1,:,:)))== nanmax(nanmean( squeeze(Fs_EV(1,:,:)) )));
freq_EV1 = freq_range(fInd_EV1); 
as_bc_EV1 = zeros(6,32);
si_EV1 = sin((2*pi*freq_EV1).*times); ci_EV1 = cos((2*pi*freq_EV1).*times);
fInd_EV2 = find(nanmean(squeeze(Fs_EV(2,:,:)))== nanmax(nanmean( squeeze(Fs_EV(2,:,:)) )));
freq_EV2 = freq_range(fInd_EV2); 
as_bc_EV2 = zeros(6,32);
si_EV2 = sin((2*pi*freq_EV2).*times); ci_EV2 = cos((2*pi*freq_EV2).*times);
fInd_EV3 = find(nanmean(squeeze(Fs_EV(3,:,:)))== nanmax(nanmean( squeeze(Fs_EV(3,:,:)) )));
freq_EV3 = freq_range(fInd_EV3); 
as_bc_EV3 = zeros(6,32);
si_EV3 = sin((2*pi*freq_EV3).*times); ci_EV3 = cos((2*pi*freq_EV3).*times);

for sample = 1:num_samples
    as_bc_vol(sample,:) = alpha_vol(sample,fInd_vol)*si_vol + beta_vol(sample,fInd_vol)*ci_vol;
    as_bc_ecc(sample,:) = alpha_ecc(sample,fInd_ecc)*si_ecc + beta_ecc(sample,fInd_ecc)*ci_ecc;
    as_bc_EV1(sample,:) = alpha_EV(1,sample,fInd_EV1)*si_EV1 + beta_EV(1,sample,fInd_EV1)*ci_EV1;
    as_bc_EV2(sample,:) = alpha_EV(2,sample,fInd_EV2)*si_EV2 + beta_EV(2,sample,fInd_EV2)*ci_EV2;
    as_bc_EV3(sample,:) = alpha_EV(3,sample,fInd_EV3)*si_EV3 + beta_EV(3,sample,fInd_EV3)*ci_EV3;
end
avgfit_vol = 1/6.*nansum(as_bc_vol);
avgfit_ecc = 1/6.*nansum(as_bc_ecc);
avgfit_EV1 = 1/6.*nansum(as_bc_EV1);
avgfit_EV2 = 1/6.*nansum(as_bc_EV2);
avgfit_EV3 = 1/6.*nansum(as_bc_EV3);

f_avgFit = figure('position',[s(f) 300 350 800]); f = f+1;
subplot(511)
plot(times,avgVol_sc,times,avgfit_vol,'LineWidth',2)
ylabel('\Delta Volume (\mum^3)','FontSize',12)
title('Volume','FontSize',14)
subplot(512)
plot(times,avgEcc_sc,times,avgfit_ecc,'LineWidth',2)
ylabel('\Delta Eccentricity','FontSize',12)
title('Eccentricity','FontSize',14)
subplot(513)
plot(times,avgEV3_sc,times,avgfit_EV3,'LineWidth',2)
ylabel('\Delta length (\mum)','FontSize',12)
title('Long Axis','FontSize',14)
subplot(514)
plot(times,avgEV2_sc,times,avgfit_EV2,'LineWidth',2)
ylabel('\Delta length (\mum)','FontSize',12)
title('Middle Axis','FontSize',14)
subplot(515)
plot(times,avgEV1_sc,times,avgfit_EV1,'LineWidth',2)
ylabel('\Delta length (\mum)','FontSize',12), xlabel('Time (hr)','FontSize',12)
title('Short Axis','FontSize',14)

%% calculate K-S test p-value

[~,pE1] = kstest2(stdEV(1,:),stdEV_s(1,:),'Alpha',0.05/3);
[~,pE2] = kstest2(stdEV(2,:),stdEV_s(2,:),'Alpha',0.05/3);
[~,pE3] = kstest2(stdEV(3,:),stdEV_s(3,:),'Alpha',0.05/3);
[~,pE1c]= kstest2(stdEVc(1,:),stdEV_sc(1,:),'Alpha',0.05/3);
[~,pE2c]= kstest2(stdEVc(2,:),stdEV_sc(2,:),'Alpha',0.05/3);
[~,pE3c]= kstest2(stdEVc(3,:),stdEV_sc(3,:),'Alpha',0.05/3);
pE1, pE2, pE3
pE1c, pE2c, pE3c

%% Plot std

times_long = [times';times';times';times';times';times'];
std_all_EV1 = [squeeze(stdEV_sc(1,1,:)); squeeze(stdEV_sc(1,2,:)); squeeze(stdEV_sc(1,3,:)); ...
    squeeze(stdEV_sc(1,4,:)); squeeze(stdEV_sc(1,5,:)); squeeze(stdEV_sc(1,6,:))];
std_all_EV2 = [squeeze(stdEV_sc(2,1,:)); squeeze(stdEV_sc(2,2,:)); squeeze(stdEV_sc(2,3,:)); ...
    squeeze(stdEV_sc(2,4,:)); squeeze(stdEV_sc(2,5,:)); squeeze(stdEV_sc(2,6,:))];
std_all_EV3 = [squeeze(stdEV_sc(3,1,:)); squeeze(stdEV_sc(3,2,:)); squeeze(stdEV_sc(3,3,:)); ...
    squeeze(stdEV_sc(3,4,:)); squeeze(stdEV_sc(3,5,:)); squeeze(stdEV_sc(3,6,:))];

[p_val_E1,tbl,stats] = anova1(std_all_EV1,times_long);
[p_val_E2,tbl,stats] = anova1(std_all_EV2,times_long);
[p_val_E3,tbl,stats] = anova1(std_all_EV3,times_long);

%f=1; %restart in a new row
%scatter plot with time, mean centered data
f_svc = figure('position',[s(f) 500 1000 600]); f = f+1;
subplot(234)
scatter(times, squeeze(stdEV_sc(1,1,:)),'fill')
hold on
scatter(times, squeeze(stdEV_sc(1,2,:)),'fill')
scatter(times, squeeze(stdEV_sc(1,3,:)),'fill')
scatter(times, squeeze(stdEV_sc(1,4,:)),'fill')
scatter(times, squeeze(stdEV_sc(1,5,:)),'fill')
scatter(times, squeeze(stdEV_sc(1,6,:)),'fill')
xlabel('Time (hr)','FontSize',12); ylabel('Length (\mum)','FontSize',12)
subplot(235)
scatter(times, squeeze(stdEV_sc(2,1,:)),'fill')
hold on
scatter(times, squeeze(stdEV_sc(2,2,:)),'fill')
scatter(times, squeeze(stdEV_sc(2,3,:)),'fill')
scatter(times, squeeze(stdEV_sc(2,4,:)),'fill')
scatter(times, squeeze(stdEV_sc(2,5,:)),'fill')
scatter(times, squeeze(stdEV_sc(2,6,:)),'fill')
axis([0 80 .2 3.1])
xlabel('Time (hr)','FontSize',12); ylabel('Length (\mum)','FontSize',12)
subplot(236)
scatter(times, squeeze(stdEV_sc(3,1,:)),'fill')
hold on
scatter(times, squeeze(stdEV_sc(3,2,:)),'fill')
scatter(times, squeeze(stdEV_sc(3,3,:)),'fill')
scatter(times, squeeze(stdEV_sc(3,4,:)),'fill')
scatter(times, squeeze(stdEV_sc(3,5,:)),'fill')
scatter(times, squeeze(stdEV_sc(3,6,:)),'fill')
axis([0 80 .2 3.1])
xlabel('Time (hr)','FontSize',12); ylabel('Length (\mum)','FontSize',12)
subplot(231)
scatter(times, stdEVc(1,:),'fill')
xlabel('Time (hr)','FontSize',12); ylabel('Length (\mum)','FontSize',12); 
title('Variance of Shortest Axis','FontSize',14)
subplot(232)
scatter(times, stdEVc(2,:),'fill')
axis([0 80 .8 2])
xlabel('Time (hr)','FontSize',12); ylabel('Length (\mum)','FontSize',12); 
title('Variance of Middle Axis','FontSize',14)
subplot(233)
scatter(times, stdEVc(3,:),'fill')
axis([0 80 .8 2])
xlabel('Time (hr)','FontSize',12); ylabel('Length (\mum)','FontSize',12); 
title('Variance of Longest Axis','FontSize',14)

% make histograms variances overlayed
nbins = 20;
stdEV_13c = [stdEVc(1,:),stdEVc(1,:),stdEVc(1,:)];
stdEV_23c = [stdEVc(2,:),stdEVc(2,:), stdEVc(2,:)];
stdEV_33c = [stdEVc(3,:),stdEVc(3,:), stdEVc(3,:)];

f_vch = figure('position',[s(f) 750 1000 350]); f = f+1;
subplot(131)
[hAx,~,~] = plotyy(5:7,9:11,5:7,[20,60,50]);
hold on
[~,centers] = hist(stdEV_sc(1,:),nbins);
hist(stdEV_sc(1,:),nbins)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.5)
hist(stdEV_13c,centers)
h = findobj(gca,'Type','patch');
set(h,'facealpha',0.75);
xlabel('Standard deviation (\mum)','FontSize',12); 
title('Smallest Axis','FontSize',14);
set(hAx(1),'XLim',[0 1.2],'YLim',[0 30],  'YTick',[0, 5, 10, 15, 20, 25, 30],'FontSize',12)
set(hAx(2),'XTick',[0, .4, .8, 1.2],'YLim',[0 30/3],'YTick',[0, 2, 4, 6, 8, 10],'FontSize',12)
set(hAx,{'ycolor'},{'r';'b'})

subplot(132)
[hAx,~,~] = plotyy(5:7,9:11,5:7,[20,60,50]);
hold on
[~,centers] = hist(stdEV_sc(2,:),nbins);
hist(stdEV_sc(2,:),nbins)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.5)
hist(stdEV_23c,centers)
h = findobj(gca,'Type','patch');
set(h,'facealpha',0.75);
xlabel('Standard deviation (\mum)','FontSize',12); 
title('Middle Axis','FontSize',14);
set(hAx(1),'YLim',[0 40],'XLim',[.5 3],'YTick',[0, 10, 20, 30, 40],'FontSize',12)
set(hAx(2),'YLim',[0 40/3],'YTick',[0, 2, 4, 6, 8, 10,12],'XTick',[.5, 1, 1.5, 2, 2.5, 3],'FontSize',12)
set(hAx,{'ycolor'},{'r';'b'})

subplot(133)
[hAx,~,~] = plotyy(5:7,9:11,5:7,[20,60,50]);
hold on
[~,centers] = hist(stdEV_sc(3,:),nbins);
hist(stdEV_sc(3,:),nbins)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.5)
hist(stdEV_33c,centers)
h = findobj(gca,'Type','patch');
set(h,'facealpha',0.75);
xlabel('Standard deviation (\mum)','FontSize',12); 
title('Longest Axis','FontSize',14);
set(hAx(1),'YLim',[0 40],'FontSize',12,'XLim',[.5 3],'YTick',[0, 10, 20, 30, 40,50])
set(hAx(2),'YLim',[0 40/3],'FontSize',12,'YTick',[0, 2, 4, 6, 8, 10,12],'XTick',[.5, 1, 1.5, 2, 2.5, 3])
set(hAx,{'ycolor'},{'r';'b'})

%% Symmetrized Kullback-Leibler Divergence

KLD_Vol = NaN(num_samples);
KLD_Ecc = NaN(num_samples);
for P = 1:num_samples
    for Q = 1:num_samples
        KLD_Vol(P,Q)= SKL( avgVol_sfc(P,:) , avgVol_sfc(Q,:)) ;
        KLD_Ecc(P,Q)= SKL( avgEcc_sfc(P,:) , avgEcc_sfc(Q,:)) ;
    end
end
labels = {'Sample 1','Sample 2','Sample 3','Sample 4','Sample 5','Sample 6'};
labels2 = {'S 1','S 2','S 3','S 4','S 5','S 6'};
f_hm = figure('position',[s(f) 750 1100 350]); f = f+1;
subplot(121)
imagesc(KLD_Ecc);
colorbar; colormap(hot); 
set(gca,'XTickLabel',labels2,'YTickLabel',labels,'FontSize',14);
title('Ecentricity Divergence','FontSize',16)
subplot(122)
imagesc(KLD_Vol);
colorbar; colormap(hot); 
set(gca,'XTickLabel',labels2,'YTickLabel',labels,'FontSize',14);
title('Volume Divergence','FontSize',16)
