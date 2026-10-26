clearvars
n = 100;
Y = zeros(n,2);
U = zeros(n,2);
t = zeros(n,1);
for i=1:n
    t(i) = i*0.01;
    Y(i,:) = v(t(i));
    U(i,:) = [ 2*sin(t(i)) 1 ];
end

Vr = Y(:,1);
Vc = Y(:,2);
%plot
figure(1);
plot(t,Vr);

title('Vr')
figure(2);
plot(t,Vc);
title('Vc')
 
H = least_squares_rlc(Y(:,1),Y(:,2),U,t,n)
out1 = lsim(H(1,1),U(:,1),t);
out2 = lsim(H(1,2),U(:,2),t);
VR = out1 + out2;
out3 = lsim(H(2,1),U(:,1),t);
out4 = lsim(H(2,2),U(:,2),t);
VC = out3 + out4;

%plot
figure(3);
plot(t,VR);
title('System Response Vr')
figure(4);
plot(t,VC);
title('System Response Vc')

% Part (b): three measurements are corrupted with random values of a much larger magnitude.

rng(1);    % fixed seed so that the run can be repeated
for i=1:3
    Y(randi([1,10],1,1),:) = randi([1,1000],1,2);
end

H2 = least_squares_rlc(Y(:,1),Y(:,2),U,t,n)
out1 = lsim(H2(1,1),U(:,1),t);
out2 = lsim(H2(1,2),U(:,2),t);
VR2 = out1 + out2;
out3 = lsim(H2(2,1),U(:,1),t);
out4 = lsim(H2(2,2),U(:,2),t);
VC2 = out3 + out4;

figure(5);
plot(t,VR);
hold on;
plot(t,VR2);
legend({'initial estimation','false estimation'})
figure(6);
plot(t,VC);
hold on;
plot(t,VC2);
legend({'initial estimation','false estimation'})
