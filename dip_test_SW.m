function dip_test(data,pc_a,xlab)
%inputs: 
% data
%pc_a=1 PCA will be perfomed
%pc_a=0 no PCA
%x_lab=labelling of x-axis e.g.: {'PC1','PC2','PC3'}
%Dependencies: 
%Hartigans dip test toolbox needed see paper https://projecteuclid.org/journals/annals-of-statistics/volume-13/issue-1/The-Dip-Test-of-Unimodality/10.1214/aos/1176346577.full
%https://github.com/VH-Lab/vhlab-toolbox-matlab/blob/master/stats/hartigansdipsigniftest.m
%https://github.com/brucegcumming/Matlab/blob/master/HartigansDipTest.m
%% Hartigan's dip test with structure data
%% Do a Hartigan's Dip test with PCA first
%close all;
% define the number of bootstraps
nboot = 500;
if pc_a==1
% define the data to be tested
 [coeff,score,latent,~,explained,mu] = pca([zscore(data)]);
 var_exp(explained,[],[]); 
%  ylim([0 100]);;xlim([0 7]);xticks([1:1:6]);
%  yticks([0:25:100])
  test_data =score(:,1:3);
else  pc_a==0
    test_data =data;
end;
%% 
%test_data = vertcat(str.PCs);
% get the number of dimensions
data_dim = size(test_data,2);

% allocate memory for the results
dip_results = zeros(data_dim,2);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 800, 300]);
% check for all PCs
for i = 1:data_dim
    [dip_results(i,1), dip_results(i,2)] = HartigansDipSignifTest(test_data(:,i), nboot);
    
    % plot the result
    %subplot(round(sqrt(data_dim)),ceil(sqrt(data_dim)),i)
    subplot(1,3,i)
   h1= histogram(test_data(:,i));box off;h1.FaceColor=[0.6 0.6 0.6];h1.EdgeColor=[1 1 1];
   %ylim([0 140]);yticks([0:70:140])
    title(strjoin({'Dip test=',num2str(round(dip_results(i,1),2)),'p=',num2str(round(dip_results(i,2),2))},' '));
    ylabel('Count');
    xlabel(xlab{i})
    set(gca,'FontSize',11);
end

end