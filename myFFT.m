function freqDomain = myFFT(x, shift)     

    N = length(x);        
    n = 0:1:N-1;          
    k = 0:1:N-1;
    k = k-shift;
    WN = exp(-1i*2*pi/N);  
    nk = n'*k;            
    WNnk = WN .^ nk;      
    freqDomain = fftshift(WNnk*x);        

end