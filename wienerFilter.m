function [yhat, H] = wienerFilter(ideal,observation,noise)
% 
% filtdata = wienerFilter(ideal,observation);
%
% FFT based Wiener filter in one dimension
% 
% Given a ideal of our perfect underlying signal that
% we wish to recover, we estimate the noise from
% noise = observation-ideal;
% The filtering is then performed in the frequency 
% domain by constructing the optimal (Wiener) filter
% for this noise/ideal estimate ... H=Sf2./(Sf2/Nf2)
% See Ref: Numerical Recipes in C, chapter 13.3, Press 
% Ref: http://www.phys.uni.torun.pl/nrbook/c13-3.pdf 
%
% G.D. Clifford 2004 gari AT mit DOT edu


% work out how long to make FFT
% N = length(observation);
% if ~mod(N,2)
%     N = N+1;
% end

N=2*length(observation)-1;

% Wiener filter
%Sf2=real(fft(ideal,N)).^2;   % Smeared ideal
Sf2 = make2sided(pwelch(ideal,4000,[],N))';
%Nf2=smooth(real(fft(noise,N)).^2, N/1000)';   % noise
Nf2 = make2sided(pwelch(noise,4000,[],N))';
Cf=real(fft(observation,N)); % ~= sqrt(Sf2+Nf2); % Corrupted ideal
      
H=Sf2./(Sf2+Nf2);              % Optimal filter
Yhat=H.*Cf;     
yhat=real(ifft(Yhat));  

yhat=yhat(1:length(observation)); 


    function Fout = make2sided(F)
        Fout = [F ; flipud(F(1:end-1))];
    end

end