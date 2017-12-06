function pidx_test(RawValues, N, NonZeroIndx, cumweigths, kernel, n)
% DOCUMENT 
% --- run test script
   % cumweights = real CDF
   [compte, IndxKhat] = hist(RawValues, unique(RawValues));                 % count collisions
   Khat = zeros(N);                                                         % initialize matrix Khat, our approximate of the kernel
   Khat(IndxKhat) = compte;                                                 % assign count of collisions to the matrix
   Khat= Khat./sum(Khat(:));                                                % normalize Khat to obtain frquencies of collisions
   
   weigthsKhat = Khat(NonZeroIndx);                                         % estimate weights of non-zero values
   cumweigthsKhat = cumsum(weigthsKhat);                                    % build estimate cumsum
   Shat = cumweigthsKhat(length(cumweigthsKhat));                           % extract estimate maximal value of cumweights
   cumweigthsKhat = cumweigthsKhat./Shat;                                   % compute estimate CDF
   
   % plot QQ-plot
   figure
   plot(cumweigths, cumweigthsKhat); 
   xlabel('CDF(Khat)'); 
   ylabel('CDF(K)'); 
   title(sprintf('Q-Q plot for kernel %.0f and n = %.0f', kernel, n))
   
   D = max(abs(cumweigths - cumweigthsKhat));                               % Perform Kolmogorov-Smirnov test
   k = 1.36/sqrt(n);                                                        % K-S critical value
   fprintf('Testing for n=%.0f, kernel %.0f \n', n, kernel);
   if (D<k)
       fprintf('By Kolmogorov-Smirnov test, H0 accepted (alpha = 0.05)\n');
   else
       fprintf('By Kolmogorov-Smirnov test, H0 not accepted (alpha = 0.05)\n');
   end
end