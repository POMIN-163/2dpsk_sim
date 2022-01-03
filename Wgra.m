function  Wgra(x,Y, display)

if display == 1
    i=16;
    j=16*100;
    t=linspace(0,16,j);
    k=j/i;

    y=fft(x);
    n = length(x);
    fshift = (-n/2:n/2-1)*(k/n);
    yshift = fftshift(y);%使零频率分量出现在矩阵的中心处附近

    figure;
    subplot(211);
    plot(t,x);
    axis([0,i,-Y,Y]);
    title('时域图');

    subplot(212);
    plot(fshift,abs(yshift));
    title('对称频域图');
    axis([-15,15,-500,500]);
end