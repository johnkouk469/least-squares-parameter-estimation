function H = least_squares_rlc(Vr,Vc,U,t,n)
% Least-squares estimate of the 2x2 transfer matrix from the two inputs
% (U(:,1) a sine, U(:,2) a constant) to the measured voltages Vr and Vc.
% H(1,:) maps the inputs to Vr, H(2,:) maps them to Vc. Both rows get the
% same denominator, as the two outputs share the dynamics.

L = [ 1 10 10 ];    % filter polynomial of the linear parametrisation

H = [ estimate_row(Vr,U,t,L) ; estimate_row(Vc,U,t,L) ];
end

function Hrow = estimate_row(y,U,t,L)
% Linear parametrisation y = theta'*phi with filtered signals phi.
% The s^2 term of the sine input is left out: for a sinusoid s^2*u = -u, so
% that regressor is linearly dependent on the one for u itself, and only the
% difference of the two coefficients can be identified from this input.
phi = [ lsim(tf([-1 0],L),y,t) ...
        lsim(tf(-1,L),y,t) ...
        lsim(tf([1 0],L),U(:,1),t) ...
        lsim(tf(1,L),U(:,1),t) ...
        lsim(tf([1 0 0],L),U(:,2),t) ...
        lsim(tf([1 0],L),U(:,2),t) ...
        lsim(tf(1,L),U(:,2),t) ];

theta = phi\y;

den = [1 theta(1)+L(2) theta(2)+L(3)];
Hrow = [ tf([0 theta(3:4).'],den) tf(theta(5:7).',den) ];
end
