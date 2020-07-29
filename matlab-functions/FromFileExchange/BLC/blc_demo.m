% Example of training a Bayesian linear classifier and using it to make
% predictions.
%
% Copyright 2013 James Barrett
%
% Date: 21 November 2013
% Contact: james.j.barrett@kcl.ac.uk

figure(10); clf;

sigma = [-1, 1, -1, 1, -1]';
 
 X = [    1.8300   -1.5834
    0.2850    2.1849
    1.5594    2.4828
    0.4276    0.9980
    1.3567    1.0205];

% subset = randperm(numIR+numRed,1000);
X = double(points);
sigma = 2*(species'-1)-1;

% Infer the weights
[w, w0] = blc_train(X,sigma)

% Use the weights to make a prediction
p = blc_predict([-1 -1],w,w0)


%% Plot the data and separating line.
data = [X sigma];
[N, d] = size(data);
Nminus = sum(data(:,d)==-1);
data = sortrows(data,d);

Xminus = data(1:Nminus,1:d-1);
Xplus = data(Nminus+1:N,1:d-1);

plot(Xplus(:,1),Xplus(:,2),'bo')
hold on
plot(Xminus(:,1),Xminus(:,2),'ro')
legend('+1 class','-1 class')

x = -3:1:3;
y=-(w(1)/w(2))*x -(w0/w(2))*ones(1,length(x));
plot(x,y,'-k')

set(gca(),'YLim',[-3 3])
set(gca(),'XLim',[-3 3])
xlabel('x_1')
ylabel('x_2')
grid on