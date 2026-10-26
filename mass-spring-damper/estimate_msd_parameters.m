clearvars

%parameters
m=15;
k=2;
b=0.2;

%initial state
y0(1) = 0;
y0(2) = 0;
u=@(t)5*sin(2*t)+10.5;
yy = zeros(100,2);
t = zeros(100,1);

% Simulate the system through its state equations with ode15s and take samples of the
% displacement for 10 s with a 0.1 s step.

%simulation
for i=0:99
    [t_ode,y] = ode15s(@(t,state)msd_dynamics(t,state,m,b,k,u),[i*0.1 (i+1)*0.1],y0);
    y0(1)= y(end,1);
    y0(2) = y(end,2);
    yy(i+1,:) = y(end,:);
    t(i+1) = t_ode(end);
end

%plot
figure(1);
plot(t,yy(:,1));
xlabel('Time')
ylabel('y(t)')
title('Position')
figure(2)
plot(t,yy(:,2));
xlabel('Time')
ylabel('dy(t)/dt')
title('Velocity')

%linear parametization and least squares method

% Vector U: samples of the external force at the same time instants as the displacement.

U = zeros(100,1);
for i = 1:100
    U(i) = u(t(i));
end

Y = yy(:,1);

%plot
figure(3);
plot(t,U);
hold on;
plot(t,Y);
title('Input - Output samples')

% Build the matrix phi with lsim and choose the coefficients of the filter polynomial Lambda(s).

phi = zeros(100,3);
L = [ 1 10 10 ];
phi(:,1) = lsim(tf([-1 0],L),Y,t);
phi(:,2) = lsim(tf(-1,L),Y,t);
phi(:,3) = lsim(tf(1,L),U,t);

%plot
figure(4);
plot(t,phi(:,1));
title('\phi_1')
figure(5);
plot(t,phi(:,2));
title('\phi_2')
figure(6);
plot(t,phi(:,3));
title('\phi_3')

theta0 = inv(phi'*phi)*phi'*Y;

theta = [b/m-L(2) ; k/m-L(3) ; 1/m ];
error = theta -theta0;
error
