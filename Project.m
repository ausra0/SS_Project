% Statistical Simulation Assignment 
% Ausra Pogozelskyte (nov.-dec. 2017)

% --- OPTIMIZATION ---
% --- Pidx
% Try to optimize this varaible class changing (get rid of it)
% Try to move the test out of the function

% --- general 
% For big redundant values, save as .mat file and load it when needed
% instead of recomputing over and over
% avoid recomputing the same thing over and over again
%% Exercice 1 
clc

%generate x
N = 5; %number of particles, small to have a lot of samples per value
%x = randi(50, N, 1);
x = ones(N,1); 
n = [1e2, 1e3, 1e4, 1e5]; % number of samples wanted
test = true; 

% to cdf plot and Kolmogorov-Smirnov test 
rem = zeros(4, length(n)); 
for i = 1:length(n)
    % test kernel 1
    [~,RawValues, NonZeroIndx, cumweigths] = pidx(1, x, n(i),[]); 
    pidx_test(RawValues, N, NonZeroIndx, cumweigths, 1, n(i));

    % test kernel 2
    [~,RawValues, NonZeroIndx, cumweigths] = pidx(2, x, n(i), [], 1);
    pidx_test(RawValues, N, NonZeroIndx, cumweigths, 2, n(i));

    % test kernel 3 
    [~,RawValues, NonZeroIndx, cumweigths] = pidx(3, x, n(i), []); 
    pidx_test(RawValues, N, NonZeroIndx, cumweigths, 3, n(i));
    
    % test kernel 4
    [~,RawValues, NonZeroIndx, cumweigths] = pidx(4, x, n(i), []); 
    pidx_test(RawValues, N, NonZeroIndx, cumweigths, 4, n(i));    
end

%% Exercice 2.1 
clc
close all 

% test run of Continuous Time/Discrete State MC 
N = 1e3; 
x = ones(N,1); 
T = 10; 

tic;
[JC0, J0] = CTDSMC(1, x, T, true); 
toc;
%%
% plot of the results
figure 
plot(J0, JC0); 
axis([0, T, 0, max(max(JC0))])
title('Continuous Time/Discret State MC'); 

%% Exercice 2.2
clc
close all 

% define parameters
N = 1e3; %1e2, 1e3, 1e4
R = 1e3; 
x = ones(N,1); % x0 = (1, ..., 1)'
T = 10; % define time interval [0,T]
k = [1, 5, 15, 50]; 
p = 1.5; 
kernel = [3]; %, 4]; 

timeline = linspace(0, T, 10*N); 
for i = 1: length(kernel)
    [meanC, sigmaC, meanM, sigmaM] = CandM(kernel(i),timeline, R, x, T, k, p);
    
    figure % plot for the average particle concentration 
    for j = 1: length(k)
        subplot(2, 2, j) 
        plot(timeline, meanC(:, :, j), 'Color', [0, 114, 93]./255)
        hold on 
        plot(timeline, meanC(:, :, j) + sigmaC(:, :, j), 'Color', [0, 198, 162]./255)
        plot(timeline, meanC(:, :, j) - sigmaC(:, :, j), 'Color', [0, 198, 162]./255)
        hold off
        title(sprintf('Plot of c(t, %.0f) for kernel %.0f', k(j), kernel(i)))
        legend('mean', 'upperbound', 'lowerbound', 'Location', 'SouthEast')
    end
    
    figure % plot for the moments
    plot(timeline, meanM(:, :, 1), 'Color', [252, 7, 85]./255)
    hold on 
    plot(timeline, meanM(:, :, 1) + sigmaM(:, :, 1), 'Color', [255, 183, 193]./255)
    plot(timeline, meanM(:, :, 1) - sigmaM(:, :, 1), 'Color', [255, 183, 193]./255)
    hold off
    title(sprintf('Plot of m_{%.1f}(t) for kernel %0.f', p, kernel(i)))
    legend('mean', 'upperbound', 'lowerbound', 'Location', 'SouthEast')
    
end

%% Exercice 2.3 
clc
close all 

% define parameters
a = [0.7, 0.8, 0.9, 1]; 
kernel = 2; 
N = [1e2, 1e3, 1e4]; 
T = 10; 
R = 50; 

for j = 1: length(N) % compute for different values of N
    x = ones(N(j), 1); % define initial state
    timeline = linspace(0, T, 10*N(j)); 
    
    figure
    for i = 1: length(a) % compute for different values of a
        [~ , meanGN, sigmaGN] = Gn(kernel, timeline, x, T, R, a(i)); 
        subplot(2, 2, i)
        plot(timeline, meanGN, 'Color', [66, 134, 244]./255) 
        hold on 
        plot(timeline, meanGN+sigmaGN, 'Color', [164, 197, 252]./255)
        plot(timeline, meanGN-sigmaGN, 'Color', [164, 197, 252]./255)
        hold off
        legend('mean GN value', 'upper bound', 'lower bound', 'Location', 'southeast')
        axis([0, T, 0, 2])
        xlabel('time')
        ylabel('G_{N}(t)')
        title(sprintf('Value of G_{%.0f}(t) for a = %.1f', N(j), a(i)))
    end
end

%% Exercice 3

