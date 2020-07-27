# SFND Radar Target Generation and Detection
## [Rubric](https://review.udacity.com/#!/rubrics/2548/view) Points
---
#### 1. FMCW Waveform Design
Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (Tchirp) and slope of the chirp.

```Matlab
%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc= 77e9;         
d_res = 1;
c = 3*10^8;
RMax = 200;

B_sweep = c/(2*d_res);
T_chirp = 5.5*2*RMax/c; 
alpha = Bsweep/Tchirp; 
```

#### 2. Simulation Loop
Simulate Target movement and calculate the beat or mixed signal for every timestamp.

```Matlab
d_0 = 80;         
v_0 = -50;        

Nd=128;          
Nr=1024;         

t=linspace(0,Nd*Tchirp,Nr*Nd); 

Tx=zeros(1,length(t));       
Rx=zeros(1,length(t));       
Mix = zeros(1,length(t));     

r_t=zeros(1,length(t));
td=zeros(1,length(t));

disp(length(t))
 
for i = 1:length(t)         
    r_t(i) = d0 + v0*t(i);
    td(i) = 2*r_t(i)/c; 
    Tx(i) = cos(2*pi*(fc*t(i) + slope*t(i)^2/2));
    Rx(i) = cos(2*pi*(fc*(t(i) - td(i)) + slope*(t(i) - td(i))^2/2));
    Mix(i) = Tx(i) * Rx(i);
end
```

#### 3. Range FFT (1st FFT)

Implement the Range FFT on the Beat or Mixed Signal and plot the result.

```Matlab
sig_fft = fft(Mix,Nr)./Nr;
sig_fft = abs(sig_fft);  

sig_fft = sig_fft(1:(Nr/2));

figure ('Name','Range from First FFT')

plot(sig_fft); grid minor     % plot FFT output 
axis ([0 200 0 1]);
xlabel('measured range');
```
![result](Images/range_FFT.png)

#### 4. doppler FFT (2st FFT)

```Matlab
Mix=reshape(Mix,[Nr,Nd]);

sig_fft2 = fft2(Mix,Nr,Nd);

sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure('Name','Range Doppler Map'),surf(doppler_axis,range_axis,RDM);
```
![result](Images/2D_FFT.png)

#### 5. 2D CFAR
Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map.

Determine the number of Training cells for each dimension. Similarly, pick the number of guard cells.

```Matlab
Tcr = 10;      %the number of Training Cells in both the dimensions.
Tcd = 4;

Gcr = 5;       %the number of Guard Cells in both dimensions around the 
Gcd = 2;       %Cell under test (CUT) for accurate estimation
      
offset = 1.4;  % offset the threshold by SNR value in dB
```

Create a vector to store noise_level for each iteration on training cells，and  get the training Cells Num.

```Matlab
noise_level = zeros(Nr/2-2*(Tcd+Gcd),Nd-2*(Tcr+Gcr));

gridSize = (2*Tcr+2*Gcr+1)*(2*Tcd+2*Gcd+1);
trainingCellsNum = gridSize-(2*Gcr+1)*(2*Gcd+1);   

CFAR_sig = zeros(size(RDM)); 
```

Slide the cell under test across the complete matrix. Make sure the CUT has margin for Training and Guard cells from the edges.

```Matlab
for j=1:Nd-2*(Tcr+Gcr)
    for i=1:Nr/2-2*(Tcd+Gcd)
        trainingCellsPatch = db2pow(RDM(i:i+2*(Tcd+Gcd),j:j+2*(Gcr+Tcr)));
        trainingCellsPatch(Tcd+1:end-Tcd,Tcr+1:end-Tcr) = 0;
        
        noise_level(i,j) = pow2db(sum(sum(trainingCellsPatch))/trainingCellsNum);
        sigThresh = noise_level(i,j)*offset;
        if RDM(i+(Tcd+Gcd),j+(Tcd+Gcr))>sigThresh
            CFAR_sig(i+(Tcd+Gcd),j+(Tcd+Gcr)) = 1;
        else
            CFAR_sig(i+(Tcd+Gcd),j+(Tcd+Gcr)) = 0;
        end
           
    end
end
```

For every iteration,  convert the value from logarithmic to linear using db2pow function, then sum the signal level within all the training cells, and get the mean value. then convert it back to db using pow2db.

```Matlab
trainingCellsPatch = db2pow(RDM(i:i+2*(Tcd+Gcd),j:j+2*(Gcr+Tcr)));
trainingCellsPatch(Tcd+1:end-Tcd,Tcr+1:end-Tcr) = 0;
        
noise_level(i,j) = pow2db(sum(sum(trainingCellsPatch))/trainingCellsNum);
```

Add the offset to it to determine the threshold.

```Matlab
threshold = pow2db(sum_train/Tcell)*offset;
```

Next, compare the signal under CUT against this threshold.
If the CUT level > threshold assign it a value of 1, else equate it to 0.

```Matlab
if RDM(i+(Tcd+Gcd),j+(Tcd+Gcr))>sigThresh
    CFAR_sig(i+(Tcd+Gcd),j+(Tcd+Gcr)) = 1;
else
    CFAR_sig(i+(Tcd+Gcd),j+(Tcd+Gcr)) = 0;
end
```

Training, Guard cells and offset are selected by increasing and decreasing to match the image shared in walkthrough.


![result](Images/CA_CFAR.png)
