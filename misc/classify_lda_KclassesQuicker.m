function [cpred,Sw,L,K] = classify_lda_KclassesQuicker(xtrain, ctrain, xtest,regularization)
% Multi-class classification using linear discriminant analysis 
% without prior!!!!
% INPUT:
%   xtrain : training set, p*ctr matrix with c datapoints in p dimensions
%   ctrain : 1*ctr vector with class labels corresponding to xtrain.
%            Class labels should be in the range 1..N for N classes and
%            each class should occur equally often.
%            To improve speed, these assumptions are not checked.
%   xtest  : test set, p*cte matrix with c2 datapoints in p dimensions 
% OPTIONAL INPUT: 
%       'regularization', default is 0.01 
% OUTPUT:    
%   cpred  : 1*cte vector with predicted class labels for xtest.
%
% Note: 
%
% JD June 2010

if (nargin<4) 
    regularization=0.01;
end;

[P,N]=size(xtrain);     % size of training set 
classes=min(ctrain):max(ctrain); % classses we do classification on
cc=numel(classes);     % class count

muK=zeros([P cc]);      % means
Sw=zeros(P,P);          % Within class variability 


%-------------calculate Parameter-----------------
for i=1:cc;
    j = find(ctrain==classes(i));  % select datapoints in this class
    n = length(j);                 % number of sampels per category 
    muK(:,i) = sum(xtrain(:,j),2)/n; % get the Cluster means 
    res = bsxfun(@minus,xtrain(:,j),muK(:,i));
    Sw = Sw+res*res'; % Estimate common covariance matrix
end;
%-------------Regularisation----------------------
P=size(Sw,1);
Sw=Sw/N; 
Sw=Sw+eye(P)*trace(Sw)*regularization/P;

%-------------classify----------------------------
[dummy N_xtest] = size(xtest);
L=muK'/Sw;                 % Calculate the classifier L(i,:)=muK'*inv(Sw)
K=sum(L.*muK',2);          % constant term for each class muK'*inv(Sw)*muK
G=bsxfun(@plus,-0.5*K,L*xtest);  % Classification function for each class
[gmax,idx]=max(G);         % nearest to class mean
cpred=classes(idx);        % predictor