% 从AudioDataFile.mat文件中读取音频的PCM数据，完成设计系统仿真

clear all  % 将所有工作空间内的变量清掉

load AudioDataFile.mat % 将音频数据Aud_data变量载入工作空间，得到的数据进行调制和传输

data_len = length(Aud_data)
global data_len;

% 自己编写以下程序（中间可以画过程波形，丰富课设内容）：
% 数字调制
% 信道传输
% 数字解调
% 完成以上程序，最后接收端要得到接收到的Aud_data数据，后面程序是将数据恢复成为样值序列，送给声卡播放。

xx = Aud_data; % 假设直接接收到的就是发送端发送的无误码数据，课设设计时应是传输到接收端恢复出的数据。

for ii = 1:size(xx, 1) % 将16位二进制码（一行16个比特构成的向量，最高位为符号位），转换为十进制的采样值
    xx_r(ii, 1) = sum(xx(ii, 2:16) .* (2.^(14:-1:0)));

    if xx(ii, 1) == 1
        xx_r(ii, 1) = -xx_r(ii, 1); % 得到样值序列，存入xx_r变量，如果最高位为1，则样值为负。
    end

end

Fs = 16000; % 采样速率为16000，采样率由音频文件决定
t = (0:1:length(xx_r) - 1) * (1 / Fs); % 采样时刻向量
plot(t, xx_r); % 画出音频信号波形
sound(xx_r, Fs); % 将样值送到声卡播放，采样率为Fs



% 后面自己编写代码：
% 完成系统有效性（画出传输信号功率谱密度，分析带宽）分析
% 完成系统可靠性（计算误码率）分析
