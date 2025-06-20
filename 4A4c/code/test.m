poles = [-0.5234+1j*1.217, -0.5234-1j*1.217, -0.00247+1j*0.08988, -0.00247-1j*0.08988];
zeros = [0, -11.32, -1.370, -0.01054];

EtoC = tf(poly(zeros), poly(poles));

figure
rlocus(EtoC)