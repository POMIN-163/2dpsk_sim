%% Author: 张海良
%% Date: 2021-12-20 20:53:43
%% Description: 通信原理课设：数字数字调制通信系统设计 -- 2DPSK 相干解调数字传输系统设计
%%

clear all  % 将所有工作空间内的变量、窗口清掉
close all
tic
%{ ---- ---- ---- 参数列表 ---- ---- ---- %}

unit_len = 16;           % 信息长度
fs = 1600;               % 采样率
N = fs / unit_len;       % 单个码元区间内的采样点数
t = linspace(0, 16, fs); % 将区间 [0,16] 分成 16000 份采样点

%{ ---- ---- ---- 读入音频 ---- ---- ---- %}

load 'AudioDataFile.mat' % 将音频数据Aud_data变量载入工作空间，得到的数据进行调制和传输

data_len = 160001 % 音频文件长度

%{ ---- ---- ---- 数字调制 ---- ---- ---- %}

fc = 2 * pi * 5;           % 载波频率 = 5HZ
wave_0 = sin(t * fc);      % 载波
wave_1 = sin(t * fc + pi); % 移相后载波
noise_n = 1;               % 高斯信道信噪比
data_rec = zeros(data_len, unit_len); % 解调接收

for times = 1 : data_len
    unit_tran = Aud_data(times, :); % 一个码元

    W1 = zeros(1, fs);

    for n = 1:unit_len

        if unit_tran(n) < 1

            for m = N * (n - 1) + 1:N * n
                W1(m) = 0;
            end

        else

            for m = N * (n - 1) + 1:N * n
                W1(m) = 1;
            end

        end

    end

    if times == 1
        figure;
        subplot(211);
        plot(t, W1);
        title('绝对码');
        axis([0, unit_len, -1, 2]);
    end

    % 产生相对码: 以 0 为基准，相对码和绝对码取亦或，得到相对码的下一位。即当前相对码等于当前绝对码和前一个相对码的异或。
    compare_tran = zeros(1, unit_len);
    compare_tran(1) = unit_tran(1);

    for n = 2:unit_len

        if unit_tran(n) == compare_tran(n - 1)
            compare_tran(n) = 0;
        else
            compare_tran(n) = 1;
        end

    end

    W2 = zeros(1, fs);

    for n = 1:unit_len

        if compare_tran(n) == 0

            for m = N * (n - 1) + 1:N * n
                W2(m) = 0;
            end

        else

            for m = N * (n - 1) + 1:N * n
                W2(m) = 1;
            end

        end

    end

    if times == 1
        subplot(212);
        plot(t, W2);
        title('相对码');
        axis([0, unit_len, -1, 2]);

        % 载波信息
        figure;
        subplot(211);
        plot(t, wave_1);
        title('载波');
        axis([0, 1, -1.5, 1.5]);

        subplot(212);
        plot(t, wave_0);
        title('移相半个周期后载波');
        axis([0, 1, -1.5, 1.5]);
    end

    % 2pdsk键控法调制：若相对码为0，则加载载波，若相对码为1，则加载移相后载波，实现 2dpsk 信号调制。
    W3 = zeros(1, fs);

    for n = 1:unit_len

        if compare_tran(n) == 0

            for m = N * (n - 1) + 1:N * n
                W3(m) = wave_1(m);
            end

        else

            for m = N * (n - 1) + 1:N * n
                W3(m) = wave_0(m);
            end

        end

    end

    Wgra(W3, 2, times);

    %{ ---- ---- ---- 信道传输 ---- ---- ---- %}

    % 高斯信道加噪：添加信噪比为 noise_n dB的加权高斯白噪声
    W4 = awgn(W3, noise_n);
    Wgra(W4, 5, times);

    %{ ---- ---- ---- 数字解调 ---- ---- ---- %}

    % 相干解调：2dpsk信号与本地载波相乘
    W5 = W4 .* wave_1;

    Wgra(W5, 5, times);

    % (弃用，效率低) 低通滤波器：根据上一步的频域图分析得，应采用截止频率为4的低通滤波器。
    % (弃用，效率低) 低通滤波器的设计利用fdatool工具实现，导出源码见dtlbq.m文件。
    % Hd = dtlbq;
    % lp = filter(Hd, W5);
    % W6 = lp;
    % [f,af] = T2F(t,W5);
    % [t,W6] = lpf(f,af,4);

    fp = 50;
    fs = 120;
    rp = 3;
    rs = 20;
    fn = 11025;
    ws = fs / (fn / 2);
    wp = fp / (fn / 2);            % 计算归一化角频率
    [n,wn] = buttord(wp,ws,rp,rs); % 计算阶数和截止频率
    [b,a] = butter(n, wn);         % 计算H(z)
    W6 = filter(b, a, W5);
    Wgra(W6, 1, times);
    % 抽样判决:低通滤波器有延迟，且低通滤波后的波形不平整。
    % 由于每个码元区间的中点信号最稳定，于是选择每个码元区间的中点作为抽样点，根据正负值进行判决。

    fs = 1600;
    samp = zeros(1, unit_len);
    W7 = zeros(1, fs);

    if W6(N / 2) > 0
        samp(1) = 1;

        for m = 1:N
            W7(m) = 1;
        end

    else
        samp(1) = 0;

        for m = 1:N
            W7(m) = 0;
        end

    end

    for n = 1:unit_len - 1

        if W6(n * N + N / 2) > 0
            samp(n + 1) = 1;

            for m = n * N:n * N + N - 1
                W7(m) = 1;
            end

        else
            samp(n + 1) = 0;

            for m = n * N:n * N + N - 1
                W7(m) = 0;
            end

        end

    end

    if times == 1
        figure;
        subplot(211);
        plot(t, W7);
        axis([0, unit_len, -1, 2]);
        title('抽样判决波形');
    end

    % 相对码转绝对码：对相对码数组逐项判断，根据相对码的0，1选择是否逻辑取反。
    % 相对码中，1表示绝对码当前位与相对码前一位不同，因此当前相对码值为1时，当前绝对码就等于上一位相对码的逻辑取反。
    compare_rec = zeros(1, unit_len);

    if (samp(1) == 1)
        compare_rec(1) = 0;
    else
        compare_rec(1) = 1;
    end

    for n = 2:unit_len

        if (samp(n) == 1)

            if (samp(n - 1) == 1)
                compare_rec(n) = 0;
            else
                compare_rec(n) = 1;
            end

        else
            compare_rec(n) = samp(n - 1);
        end

    end

    W8 = zeros(1, fs);

    for n = 1:unit_len

        if compare_rec(n) == 1

            for m = N * (n - 1) + 1:N * n
                W8(m) = 1;
            end

        else

            for m = N * (n - 1) + 1:N * n
                W8(m) = 0;
            end

        end

    end

    if times == 1
        subplot(212);
        plot(t, W8);
        axis([0, unit_len, -1, 2]);
        title('相对码转绝对码');
    end
    % 信号频域分析
    if times < 100
        figure(8);
        hold on;
        L = 100; %信号长度
        N = 2^nextpow2(L); %取最接近1000的2的幂次,实现快速运算
        sf = fft(W3, N) / 100; %由于时域抽样会有一个 1/Ts的衰减，所以必须乘以Ts也即除以fs
        df = 100 / N; %频率分辨率
        f = [0:df:df * (N - 1)] - 100 / 2; %使频域轴对称
        N = 100;
        plot(f, fftshift(abs(sf))); grid on; %把数据‘循环倒转’
        axis([-50 50 0 1]);
        title('信号功率谱');
    end

    x1 = zeros(1, 16);
    for m1 = 0 : 16 - 1 % 提取样点值恢复原序列
        data_rec(times, m1 + 1) = W8(1, 100 * m1 + 50);
    end

    % 显示计算进度
    if rem(times, 100) == 0
        disp('已完成: ' + string(times) + ' 总数: ' + string(data_len));
    end
end

xx = data_rec; % 假设直接接收到的就是发送端发送的无误码数据，课设设计时应是传输到接收端恢复出的数据。

for ii = 1:size(xx, 1) % 将16位二进制码（一行16个比特构成的向量，最高位为符号位），转换为十进制的采样值
    xx_r(ii, 1) = sum(xx(ii, 2:16) .* (2.^(14:-1:0)));

    if xx(ii, 1) == 1
        xx_r(ii, 1) = -xx_r(ii, 1); % 得到样值序列，存入xx_r变量，如果最高位为1，则样值为负。
    end

end

Fs = 16000; % 采样速率为16000，采样率由音频文件决定
t = (0:1:length(xx_r) - 1) * (1 / Fs); % 采样时刻向量
figure;
plot(t, xx_r); % 画出音频信号波形
sound(xx_r, Fs); % 将样值送到声卡播放，采样率为Fs

N = 16;     % 产生信码长度
Ns = 100;   % 一个码元内的信号点数
Ts = 1;     % 码元时间为1s
t = 0:Ts / Ns:(Ts*N-Ts / Ns); % 产生横坐标时间向量
A = 1;                        % 信号幅度为1

% 后面自己编写代码：
% 完成系统有效性（画出传输信号功率谱密度，分析带宽）分析

pause(10);
% close all;
% 完成系统可靠性（计算误码率）分析

unit_len = 16;           % 信息长度
fs = 16000*1.778;            % 采样率
N = fs / unit_len;       % 单个码元区间内的采样点数
t = linspace(0, 16, fs); % 将区间 [0,16] 分成 1600 份采样点

%{ ---- ---- ---- 读入音频 ---- ---- ---- %}

load 'AudioDataFile.mat' % 将音频数据Aud_data变量载入工作空间，得到的数据进行调制和传输

fc = 2 * pi * 5;           % 载波频率 = 5HZ
wave_0 = sin(t * fc);      % 载波
wave_1 = sin(t * fc + pi); % 移相后载波
noise_n = 1;               % 高斯信道信噪比
data_rec = zeros(data_len, unit_len); % 解调接收

SNR_db = -10:0.5:6;
for i_SNR = 1:length(SNR_db)
    Nerr = 0;
    for times = 1:1000

        x = Aud_data(times, :); % 一个码元
        W1 = zeros(1, fs);

        for n = 1:unit_len
            if x(n) < 1
                for m = N * (n - 1) + 1:N * n
                    W1(m) = 0;
                end
            else
                for m = N * (n - 1) + 1:N * n
                    W1(m) = 1;
                end
            end
        end

        % 产生相对码: 以 0 为基准，相对码和绝对码取亦或，得到相对码的下一位。即当前相对码等于当前绝对码和前一个相对码的异或。
        compare_tran = zeros(1, unit_len);
        compare_tran(1) = x(1);

        for n = 2:unit_len
            if x(n) == compare_tran(n - 1)
                compare_tran(n) = 0;
            else
                compare_tran(n) = 1;
            end
        end

        W2 = zeros(1, fs);

        for n = 1:unit_len
            if compare_tran(n) == 0
                for m = N * (n - 1) + 1:N * n
                    W2(m) = 0;
                end
            else
                for m = N * (n - 1) + 1:N * n
                    W2(m) = 1;
                end
            end
        end

        % 2pdsk键控法调制：若相对码为0，则加载载波，若相对码为1，则加载移相后载波，实现 2dpsk 信号调制。
        W3 = zeros(1, fs);

        for n = 1:unit_len
            if compare_tran(n) == 0
                for m = N * (n - 1) + 1:N * n
                    W3(m) = wave_1(m);
                end
            else
                for m = N * (n - 1) + 1:N * n
                    W3(m) = wave_0(m);
                end
            end
        end

        %{ ---- ---- ---- 信道传输 ---- ---- ---- %}

        % 高斯信道加噪：添加信噪比为 noise_n dB的加权高斯白噪声
        W4 = awgn(W3, SNR_db(i_SNR), 'measured');

        %{ ---- ---- ---- 数字解调 ---- ---- ---- %}

        % 相干解调：2dpsk信号与本地载波相乘
        W5 = W4 .* wave_1;

        % (弃用，效率低) 低通滤波器：根据上一步的频域图分析得，应采用截止频率为4的低通滤波器。
        % (弃用，效率低) 低通滤波器的设计利用fdatool工具实现，导出源码见dtlbq.m文件。
        % Hd = dtlbq;
        % lp = filter(Hd, W5);
        % W6 = lp;
        % [f,af] = T2F(t,W5);
        % [t,W6] = lpf(f,af,4);
        % fp = 500;
        % fs = 700;
        % rp = 4;
        % rs = 20;
        % fn = 11025;
        % ws = fs / (fn / 2);
        % wp = fp / (fn / 2);            % 计算归一化角频率
        % [n,wn] = buttord(wp,ws,rp,rs); % 计算阶数和截止频率
        % [b,a] = butter(n,wn);          % 计算H(z)
        % W6 = filter(b,a,W5);
        wp = 5 / 100; ws = 20 / 100;
        [n2, wn] = buttord(wp, ws, 1, 50); %阻带衰减大于50dB，通带波纹小于1dB
        [b, a] = butter(n2, wn);
        W6 = filter(b, a, W5);

        % 抽样判决:低通滤波器有延迟，且低通滤波后的波形不平整。
        % 由于每个码元区间的中点信号最稳定，于是选择每个码元区间的中点作为抽样点，根据正负值进行判决。
        samp = zeros(1, unit_len);
        W7 = zeros(1, fs);

        if W6(N / 2) > 0
            samp(1) = 1;

            for m = 1:N
                W7(m) = 1;
            end

        else
            samp(1) = 0;

            for m = 1:N
                W7(m) = 0;
            end

        end

        for n = 1:unit_len - 1
            if W6(n * N + N / 2) > 0
                samp(n + 1) = 1;
                for m = n * N:n * N + N - 1
                    W7(m) = 1;
                end
            else
                samp(n + 1) = 0;
                for m = n * N:n * N + N - 1
                    W7(m) = 0;
                end
            end
        end

        % 相对码转绝对码：对相对码数组逐项判断，根据相对码的0，1选择是否逻辑取反。
        % 相对码中，1表示绝对码当前位与相对码前一位不同，因此当前相对码值为1时，当前绝对码就等于上一位相对码的逻辑取反。
        compare_rec = zeros(1, unit_len);

        if (samp(1) == 1)
            compare_rec(1) = 0;
        else
            compare_rec(1) = 1;
        end

        for n = 2:unit_len
            if (samp(n) == 1)
                if (samp(n - 1) == 1)
                    compare_rec(n) = 0;
                else
                    compare_rec(n) = 1;
                end
            else
                compare_rec(n) = samp(n - 1);
            end
        end

        W8 = zeros(1, fs);

        for n = 1:unit_len
            if compare_rec(n) == 1
                for m = N * (n - 1) + 1:N * n
                    W8(m) = 1;
                end
            else
                for m = N * (n - 1) + 1:N * n
                    W8(m) = 0;
                end
            end
        end

        for m1 = 1 : 16 % 提取样点值恢复原序列
            data_rec(times, m1) = compare_rec(m1);
        end

        Nerr0 = sum(abs(Aud_data(times, :) - data_rec(times, :))); % 累计误码
        % if Nerr0 > 0
        %     disp('ERR');
        % end
        Nerr = Nerr + Nerr0;
    end
    Pe(1, i_SNR) = Nerr / (16 * 1000); % 抽取码元数

    disp('误码率计算已完成: ' + string(i_SNR) + ' / 33');
end
figure;
semilogy(SNR_db, Pe, 'b'); % 画出误码率曲线，第三个参数为颜色选择
hold on                    % 保持窗口

% 计算理论误码率
SNR = 10.^(SNR_db / 10);          % 将信噪比分贝值换算成一般值
Pe_theory = 1 / 2*exp(-SNR);      % 理论计算误码率
semilogy(SNR_db, Pe_theory, 'r'); % 画出理论误码率曲线，第三个参数为颜色选择

grid on
xlabel('SNR (dB)'); % 横坐标标注
ylabel('Pe');       % 纵坐标标注
title('误码率曲线'); % 设置标题
toc