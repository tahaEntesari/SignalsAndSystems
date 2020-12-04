function [y] = SoundMaker(gain,note,interval,Fs)
syms t k;
y=0;
thesines=@(k,t)gain(k).*sin(2*pi*440*2.^(note(k)/12)*t);
ts=@(k)0:1/Fs:interval(k);
for l=1:numel(interval)
       y=[y,thesines(l,ts(l))];
end
end