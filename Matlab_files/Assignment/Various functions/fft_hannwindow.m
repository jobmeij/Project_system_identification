function [FFTdata,w] = fft_hannwindow(Data,Nwindows)

DataLength = length(Data);
Overlap = 0.5;
WindowLength = round(DataLength/Nwindows);
Window = hann(WindowLength);
nShifts = Nwindows*(1/(1-Overlap))-1;

Ts = 1;
ws = (2*pi)/(Ts);
l = 0:WindowLength/2;
wvect = (l/WindowLength)*ws;

Xdata = zeros(nShifts,floor(WindowLength/2)+1);
for n = 1:nShifts
    %ceil((n-1)*WindowLength*Overlap+1)
    %ceil((n+1)*WindowLength*Overlap)
    DataPart = Window .* Data(ceil((n-1)*WindowLength*Overlap+1):ceil((n+1)*WindowLength*Overlap));
    X = fft(DataPart);
    X = [X(1); X(2:(end/2)+1)*2]/WindowLength;
    Xdata(n,:) = X';
end
FFTdata = mean(Xdata,1);
w = wvect;
end

