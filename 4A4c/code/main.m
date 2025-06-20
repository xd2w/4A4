close all;
load('Phugoid.mat', '-mat')
% load("Roll-Subs.mat")
% load("Dutch-Roll.mat")
% load("SPPO.mat")
% load('Spiral.mat')
% load('SPPO_deltae_to_theta.mat')

q = Nz;
tp = Time;
ele = Elevator;

start = 16.56;

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
% 
% figure
% plot(t, q)
% hold on;
% plot(t, ele)

[omega,xfun] = xfer(tp, ele, q);
[coef, xfit] = fitxf(omega, xfun, 4, 3, 3);

coef = coef(1, end-1);

figure
subplot(2, 1, 1);
semilogy(omega, abs(xfun))
hold on
semilogy(omega, abs(xfit))
% xlim([0, 100])

subplot(2, 1, 2);
plot(omega, angle(xfun))
hold on
plot(omega, angle(xfit))
% xlim([0, 100])


num = [1.133, 0.6522, 0.01197];
% EtoAlpha = tf(coef, [1, 1.052, 1.768, 0.01714, 0.01419]);
% EtoAlpha747 = tf(num, [1, 1.052, 1.768, 0.01714, 0.01419]);

% figure
% step(EtoAlpha)
% step(EtoAlpha747)

s = tf('s');
num_tf = tf(coef, 1);
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

EtoQ = -s*num_tf * long_num;

pole(EtoQ)
zero(EtoQ)

%%%
figure
[y, time] = step(EtoQ);
subplot(2, 1, 1);
semilogy(time, abs(y))
subplot(2, 1, 2);
plot(time, angle(y))

% figure
% bode(EtoQ)tf(1, real(poly(long_poles)), 1);

figure
step(EtoQ)

ele = ele';
tp = tp';
figure
[yy, tt] = lsim(EtoQ, Elevator(:, 1), Time(:, 1));
% plot(tt(:, 1), yy());



