function [yhat, H] = wienerFilterKP(S2,FFTdata,FFTnoise)
%S2=real(FFTexpected).^2; 
N2=real(FFTnoise).^2;   
D=real(FFTdata);
      
H=S2./(S2+N2);              % Optimal filter
Yhat=H.*D;     
yhat=real(ifft(Yhat));  
end