close all;
% load('Phugoid.mat', '-mat')
load("Roll-Subs.mat")
% load("Dutch-Roll.mat")
% load("SPPO.mat")
% load('Spiral.mat')
% load('SPPO_deltae_to_theta.mat')

s = tf('s');
% num_tf = tf(1, 1);


long_poles = [
                1.836*exp(1j*(pi - acos(0.441))),
                1.836*exp(1j*(pi + acos(0.441))),
                0.123*exp(1j*(pi + acos(0.066))),
                0.123*exp(1j*(pi - acos(0.066))),
             ];

lat_poles = [
                -1/0.2495,
                1.519*exp(1j*(pi + acos(0.15))),
                1.519*exp(1j*(pi - acos(0.15))),
                1/42.6,
             ];

% figure
% plot(long_poles, "x")
% hold on;
% plot(lat_poles, "x")

long_num = tf(1, real(poly(long_poles)));
lat_num  = tf(1, real(poly(lat_poles)));

q = Rollrt;
tp = Time;
ele = Aileron;

start = 1.9;
% start = ;

% q = q(tp>start);
% ele= ele(tp>start);
% t = tp(tp>start);

q =q(:) - mean(q(tp<start));
ele = ele(:) - mean(ele(tp<start));
t = tp;
ndat = length(t);
delt = (t(ndat)-t(1))/(ndat-1);
dterrs = diff(t)/delt - 1;

% figure
% plot(t(2:end), dterrs)

figure
plot(t, q, 'DisplayName', 'Roll Rate')
hold on;
plot(t, ele, 'DisplayName', 'Aileron Rangle')
grid on
legend()
ylabel("Normalised Amplitude")
xlabel("Time")
saveas(gcf,'figs/AtoP_IO','epsc')

%%
freq_lim = 3;
[omega,xfun] = xfer(tp, ele, q);
[coef, xfit] = fitxf(omega, xfun, real(poly(lat_poles)), freq_lim, 3);
coef

coef = coef(1:end-1);
num_tf = tf(real(coef), 1);

figure
subplot(2, 1, 1);
semilogy(omega, abs(xfun), 'DisplayName','FFT')
hold on
semilogy(omega, abs(xfit), 'DisplayName','Fit')
xlim([0, freq_lim*1.5])
xlabel("Frequency (\omega) [rads/s]")
ylabel("Magnitude")
xline(freq_lim, '-', ...
    ['\omega_{lim} = ', num2str(freq_lim), ' rads/s'], ...
    'FontSize',16, ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','top', ...
    'LabelOrientation','horizontal', ...
    'Color', 	"#D95319")
grid on

subplot(2, 1, 2);
plot(omega, angle(xfun))
hold on
plot(omega, angle(xfit))
xlim([0, freq_lim*1.5])
xlabel("Frequency (\omega) [rads/s]")
ylabel("Phase")
xline(freq_lim, '-', ...
    ['\omega_{lim} = ', num2str(freq_lim), ' rads/s'], ...
    'FontSize',16, ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','top', ...
    'LabelOrientation','horizontal', ...
    'Color', 	"#D95319")
grid on
saveas(gcf,'figs/AtoP_FFT','epsc')


% num = [1.133, 0.6522, 0.01197];
% EtoAlpha = tf(coef, [1, 1.052, 1.768, 0.01714, 0.01419]);
% EtoAlpha747 = tf(num, [1, 1.052, 1.768, 0.01714, 0.01419]);

% figure
% step(EtoAlpha)
% step(EtoAlpha747)


EtoPtch = num_tf * lat_num;
EtoQ = s*EtoPtch;

pole(EtoQ)
zero(EtoQ)

%%
% figure
% [y, time] = step(EtoQ);
% subplot(2, 1, 1);
% plot(time, abs(y))
% subplot(2, 1, 2);
% plot(time, angle(y))

load Spiral.mat
q = Rollrt;
tp = Time;
ele = Aileron;

start = 1.9;

q =q(:) - mean(q(tp<start));
ele = ele(:) - mean(ele(tp<start));


% figure
% bode(EtoQ)tf(1, real(poly(long_poles)), 1);

% figure
% step(EtoQ)

figure
plot(tp, q, 'DisplayName', 'Roll Rate')
hold on
% 
% fixed_ele = Elevator - mean(Elevator(Time < start));
[yy, tt] = lsim(EtoQ, ele, tp);
[n, ~] = size(yy);
plot(tt(1:n-1, 1), yy(1:n-1, 1),'DisplayName', 'Simulated')
% hold on
% plot(tp, ele)
hold on
grid on
legend

ylabel("Amplitude")
xlabel("Time")

saveas(gcf,'figs/AtoP_Sim','epsc')

figure
plot(tp, Rollang - mean(Rollang(tp<start)), 'DisplayName', 'Roll Angle')
hold on
[yy, tt] = lsim(EtoPtch, ele, tp);
[n, ~] = size(yy);
plot(tt(1:n-1, 1), yy(1:n-1, 1), 'DisplayName', 'Simulated')

% hold on
% plot(tp, ele, 'DisplayName', 'Aileron')
hold on

legend
ylabel("Amplitude")
xlabel("Time")
grid on
saveas(gcf,'figs/AtoPhi_Sim','epsc')

figure
rlocus(-EtoPtch);
grid on
saveas(gcf,'figs/AtoP_rlocus','epsc')

figure
bode(EtoPtch);
grid("minor")
saveas(gcf,'figs/AtoP_bode','epsc')



