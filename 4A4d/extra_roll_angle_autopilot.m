cla; close all
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

% my vals
E_theta = -tf([2.6855, 2.4735, 0.0903], real(poly(long_poles)));
E_Nz = -tf([0.0139, 0.0694, 0.4066, 0.0244, 0], real(poly(long_poles)));
R_yaw = -tf([1.8899, 6.5343, 2.6795, -0.5732], real(poly(lat_poles)));
A_Roll = -tf([5.084, 2.9062, 12.0011], real(poly(lat_poles)));

lag_e = tf(1, [0.1, 1]);
lag_a = tf(1, [0.05, 1]);
lag_r = tf(1, [0.25, 1]);

% % given vals
% E_theta = -tf([2.637, 2.475, 0.0607], real(poly(long_poles)));
% E_Nz = -tf([0.0139, 0.0693, 0.4071, 0.0242, 0], real(poly(long_poles)));
% R_yaw = -tf([2.0959, 6.0649, 3.565, -1.1943], real(poly(lat_poles)));
A_Roll = -tf([4.748, 2.957, 11.16], real(poly(lat_poles)));

% figure
% hold on;
% step(E_theta);
% rlocus(-lag_e*E_theta);
% plot_req(0.4, 1);
% hold off
% step(E_Nz);
% step(R_r);
% step(A_Ptch);
figure(3)

sys = A_Roll;

step(A_Roll, 100)

% rlocus(-lag_a*sys);

%%



s = tf('s');
kp = 0.75;
% kp=0
ki = 0;
kd = 0;
omega = 0;
omega = 1.519*sqrt(1-0.15^2); % my vals
dw = 1;

% K = kp + ki/s + kd*s;
K = 100*(s+0.5)/(s+30)
% band stop filter
% K = (s^2 + dw*s+omega^2)/(s^2 + omega^2)
washout_filter = (s/(s+omega));
band_pass = (s^2 + dw*s+omega^2)/(s^2 + omega^2);
% band_pass = washout_filter
% rlocus(-K*sys)


step(feedback(-K*lag_a*sys, band_pass), 100)
stepinfo(feedback(-K*lag_a*sys, band_pass))
% plot_c_star(0.95)
%
plot_lat_rlocus(-band_pass*lag_a*sys,K)

figure
bode(feedback(-K*lag_a*sys, band_pass))

function plot_c_star(k)
[cst,cslo,csup] = csenv(k);

% figure
hold on;
plot(cst, cslo, "Color","b");
plot(cst, csup, "Color","b");
end

function plot_long_rlocus(sys, K)
figure(2)
hold on
% phugoid
subplot(1, 2, 1);
rlocus(K*sys)
hold on
r = rlocus(K*sys, 1);
plot(real(r), imag(r), 'k+')
hold on
plot_req(0.05, 0)
xlim([-0.3, 0.3]);
ylim([-0.3, 0.3]);

% spo
subplot(1, 2, 2);
rlocus(K*sys)
hold on
r = rlocus(K*sys, 1);
plot(real(r), imag(r), 'k+')
hold on
plot_req(0.4, 1)
xlim([-1.3, 0]);
ylim([-4, 4]);

end


function plot_lat_rlocus(sys, K)
figure(2)
hold on
% dutch-roll
% subplot(1, 3, 1);
rlocus(K*sys)
hold on
r = rlocus(K*sys, 1);
plot(real(r), imag(r), 'k+')
hold on
plot_req_wc(0.4, 0.4, 'b')
% xlim([-0.3, 0.3]);
% ylim([-0.3, 0.3]);

% Roll-subsidance
% subplot(1, 3, 2);
rlocus(K*sys)
hold on
r = rlocus(K*sys, 1);
plot(real(r), imag(r), 'k+')
hold on
plot_req_wc(0.8, 1, 'g')
% xlim([-1.3, 0]);
% ylim([-4, 4]);

% spiral
% subplot(1, 3, 3);
rlocus(K*sys)
hold on
r = rlocus(K*sys, 1);
plot(real(r), imag(r), 'k+')
hold on
plot_req_wc(0, log(2)/20, 'r')
% xlim([-3, 0.5]);
% ylim([-2, 2]);

end

function plot_req_wc(zeta, wn, color)
theta = linspace(0, 2*pi, 100);
x = wn*cos(theta);
y = wn*sin(theta);
hold on;
plot(x, y, "Color",color);
r = linspace(0, 5, 2);
rx = -r*zeta;
ry = r/sqrt(1-zeta^2);
plot(rx, ry, "Color",color);
end

function plot_req(zeta, wn)
theta = linspace(0, 2*pi, 100);
x = wn*cos(theta);
y = wn*sin(theta);
hold on;
plot(x, y, "Color","r");
r = linspace(0, 5, 2);
rx = -r*zeta;
ry = r/sqrt(1-zeta^2);
plot(rx, ry, "Color","b");
end




