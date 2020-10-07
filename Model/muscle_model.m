function [out, adot, adotdot, ldot, ldotdot, vdot, vdotdot] = muscle_model(a, L, V, R, F0)
if nargin == 0  % Plot length-velocity-tension curve of the muscle model
    % Muscle parameters
    a = 1.0;
    R = 1; F0 =1;
    N = 20;
    Lmax = 1.5;
    Lmin = 0.7;
    Vmax = 2;
    Vmin = -2;
    
    L = repmat(Lmin + (Lmax - Lmin)/N * (0:N-1), [N, 1]);
    V = repmat(Vmin + (Vmax - Vmin)/N* (0:N-1)', [1, N]);
    
    Tension = muscle_model(a, L, V, R, F0);
    TensionX = muscle_model(a, ones(1, size(V,1)), V(:,1)', R, F0);
    TensionY = muscle_model(a, L(1,:), zeros(1, size(L,2)), R, F0);
    Tmax = max(TensionX);
    
    % Muscle activation dynamics
    t_act = 0.050;      % [s]
    t_deact = 0.066;   % [s]
    dt = 0.005;  % [s]
    T = 0.70;   % [s]
    duration = 0.3; % Input duration [s]
    u = zeros(3, T/dt+1);   % Neural input
    u(:, (0.1/dt)+1:((0.1+duration)/dt)) = repmat([0.3; 0.7; 1.0], [1 , duration./dt]);
    [a, time] = muscle_dyn(u, t_act, t_deact, dt);
    
    % Plot functions
    figure(1);
    subplot(1,2,1);    % Muscle activation dynamics
    plot(time-0.1, a);
    box off;
    ylim([0 1.1]);
    xlabel('Time [s]');
    ylabel('Activation level');
    xlim([-0.1 0.6]);
    legend('u = 0.3','u = 0.7','u = 1.0');
    hold off;
    subplot(1,2,2); % Length-velocity-tension curve
    mesh(V, L, Tension/Tmax);
    hold on;
    plot3(zeros(1,N), L(1,:),TensionY/Tmax,'linewidth', 2);
    plot3(V(:,1), ones(1,N) ,TensionX/Tmax,'linewidth', 2);
    xlabel('Fascicle velocity');
    ylabel('Fascicle length');
    zlabel('Normalized Tension');
    hold off;
        
else % Call to compute muscle tension
    a(a<0) = 0;
    a(a>1) = 1;
    
    szM = size(a, 1);
    szM2 = size(a, 2);
    
    if nargout>1
        [af, afdot_a, afdotdot_a, afdot_l, afdotdot_l] = Af(a, L);
        [fl, fldot_l, fldotdot_l] = FL(L);
        [fv, fvdot_l, fvdotdot_l, fvdot_v, fvdotdot_v] = FV(V, L);
        [fpe1, fpe1dot_l, fpe1dotdot_l] = FPE1(L);
        [fpe2, fpe2dot_l, fpe2dotdot_l] = FPE2(L);      
        
        adot = repmat(F0*R, [1 szM2]).*afdot_a .* (fl .* fv + fpe2);
        adotdot = repmat(F0*R, [1 szM2]).*afdotdot_a .* (fl .* fv + fpe2);        
        ldot = repmat(F0*R, [1 szM2]).*(afdot_l .* (fl .* fv + fpe2) ...
            + af .*(fldot_l.*fv  + fl.*fvdot_l + fpe2dot_l) + fpe1dot_l) ;
        ldotdot = repmat(F0*R, [1 szM2]).*(afdotdot_l .* (fl .* fv + fpe2)  ...
            + afdot_l .* (fldot_l.*fv  + fl.*fvdot_l + fpe2dot_l) ...
            + af .*(fldotdot_l.*fv  + fldot_l.*fvdot_l + fldot_l.*fvdot_l + ...
            fldot_l.*fvdot_l  +  fl.*fvdotdot_l + fpe2dotdot_l)+ fpe1dotdot_l);        
        vdot = repmat(F0*R, [1 szM2]).*(af .* (fl .* fvdot_v + fpe2)+ fpe1);
        vdotdot = repmat(F0*R, [1 szM2]).*(af .* (fl .* fvdotdot_v + fpe2)+ fpe1);
        
        adot = FixNan(adot);
        adotdot = FixNan(adotdot);
        ldot = FixNan(ldot);
        ldotdot = FixNan(ldotdot);
        vdot = FixNan(vdot);
        vdotdot = FixNan(vdotdot);
    else
        af = Af(a, L);
        fl = FL(L);
        fv = FV(V, L);
        fpe2 = FPE2(L);
        fpe1 = FPE1(L);
        adot = NaN;
        adotdot = NaN;
        ldot = NaN;
        ldotdot = NaN;
        vdot = NaN;
        vdotdot = NaN;
    end
    
    out = repmat(F0*R, [1 szM2]).*(af .* (fl .* fv + fpe2)  + fpe1);
    
end
end

function [out, adot, adotdot, ldot, ldotdot] = Af(a, L)
af = 0.56;
nf0 = 2.11;
nf1 = 3.31;
nf = nf0 + nf1*(1./L -1);

out = 1 - exp(-(a./(af.*nf)).^nf);

if  nargout>1
    adot = (exp(-(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1))) ...
        .*(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1) - 1))./af;
    adotdot =  (exp(-(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1))) ...
        .*(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1) - 2) ...
        .*(nf0 + nf1.*(1./L - 1) - 1))./(af.^2.*(nf0 + nf1.*(1./L - 1))) ...
        - (exp(-(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1))) ...
        .*(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(2.*nf0 + 2.*nf1.*(1./L - 1) - 2))./af.^2 ;    
    ldot = -exp(-(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1))) ...
        .*((nf1.*log(a./(af.*(nf0 + nf1.*(1./L - 1)))).*(a./(af.*(nf0 + nf1.*(1./L - 1)))) ...
        .^(nf0 + nf1.*(1./L - 1)))./L.^2 - (a.*nf1.*(a./(af.*(nf0 + nf1.*(1./L - 1)))) ...
        .^(nf0 + nf1.*(1./L - 1) - 1))./(L.^2.*af.*(nf0 + nf1.*(1./L - 1))));
    ldotdot = - exp(-(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1))) ...
        .*((nf1.*log(a./(af.*(nf0 + nf1.*(1./L - 1)))).*(a./(af.*(nf0 + nf1.*(1./L - 1)))) ...
        .^(nf0 + nf1.*(1./L - 1)))./L.^2 - (a.*nf1.*(a./(af.*(nf0 + nf1.*(1./L - 1)))) ...
        .^(nf0 + nf1.*(1./L - 1) - 1))./(L.^2.*af.*(nf0 + nf1.*(1./L - 1)))).^2 ...
        - exp(-(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1))) ...
        .*((nf1.^2.*(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1)))./(L.^4.*(nf0 + nf1.*(1./L - 1))) ...
        - (nf1.*log(a./(af.*(nf0 + nf1.*(1./L - 1)))).*((nf1.*log(a./(af.*(nf0 + nf1.*(1./L - 1)))) ...
        .*(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1)))./L.^2 ...
        - (a.*nf1.*(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1) - 1)) ...
        ./(L.^2.*af.*(nf0 + nf1.*(1./L - 1)))))./L.^2 - (2*nf1.*log(a./(af.*(nf0 + nf1.*(1./L - 1)))) ...
        .*(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1)))./L.^3 + (2*a.*nf1.*(a./(af.*(nf0 + nf1.*(1./L - 1)))) ...
        .^(nf0 + nf1.*(1./L - 1) - 1))./(L.^3.*af.*(nf0 + nf1.*(1./L - 1))) - (a.*nf1.^2.*(a./(af.*(nf0 + nf1.*(1./L - 1)))) ...
        .^(nf0 + nf1.*(1./L - 1) - 1))./(L.^4.*af.*(nf0 + nf1.*(1./L - 1)).^2) + (a.*nf1.*((nf1.*log(a./(af.*(nf0 + nf1.*(1./L - 1)))) ...
        .*(a./(af.*(nf0 + nf1.*(1./L - 1)))).^(nf0 + nf1.*(1./L - 1) - 1))./L.^2 - (a.*nf1.*(a./(af.*(nf0 + nf1.*(1./L - 1)))) ...
        .^(nf0 + nf1.*(1./L - 1) - 2).*(nf0 + nf1.*(1./L - 1) - 1))./(L.^2.*af.*(nf0 + nf1.*(1./L - 1)).^2))) ...
        ./(L.^2.*af.*(nf0 + nf1.*(1./L - 1))));
else
    adot = NaN;
    adotdot = NaN;
    ldot = NaN;
    ldotdot = NaN;
end
end

function [out, ldot, ldotdot] = FL(L)
w = 0.81;
b = 1.55;
p = 2.12;

out = exp(-(abs((L.^b - 1)./w)).^p);

if  nargout>1
    ldot = -(L.^(b - 1).*b.*p.*exp(-(abs(L.^b - 1)./abs(w)).^p) ...
        .*sign(L.^b - 1).*(abs(L.^b - 1)./abs(w)).^(p - 1))./abs(w);
    ldotdot = (L.^(2.*b - 2).*b.^2.*p.^2 ...
        .*exp(-(abs(L.^b - 1)./abs(w)).^p).*sign(L.^b - 1).^2 ...
        .*(abs(L.^b - 1)./abs(w)).^(2.*p - 2))./abs(w).^2 ...
        - (2.*L.^(2.*b - 2).*b.^2.*p.*exp(-(abs(L.^b - 1)./abs(w)).^p) ...
        .*dirac(L.^b - 1).*(abs(L.^b - 1)./abs(w)).^(p - 1))./abs(w) ...
        - (L.^(b - 2).*b.*p.*exp(-(abs(L.^b - 1)./abs(w)).^p).*sign(L.^b - 1) ...
        .*(abs(L.^b - 1)./abs(w)).^(p - 1).*(b - 1))./abs(w) ...
        - (L.^(2.*b - 2).*b.^2.*p.*exp(-(abs(L.^b - 1)./abs(w)).^p) ...
        .*sign(L.^b - 1).^2.*(abs(L.^b - 1)./abs(w)).^(p - 2).*(p - 1))./abs(w).^2;
else
    ldot = NaN;
    ldotdot = NaN;
end
end

function [out, ldot, ldotdot, vdot, vdotdot] = FV(V, L)
Vmax = -7.39;
Cv0 = -3.21;
Cv1 = 4.17;
av0 = -1.53;
av1 = 0;
av2 = 0;
bv = 1.05;

tmp = V < 0;
out = (bv-(av0+av1.*L+av2.*(L.^2)).*V)./(bv+V);
out(tmp) = (Vmax-V(tmp))./(Vmax+(Cv0+Cv1.*L(tmp)).*V(tmp));

if  nargout>1
    ldot = -(V.*(av1 + 2*L.*av2))./(V + bv);
    ldot(tmp) = (Cv1.*V(tmp).*(V(tmp) - Vmax))./(Vmax + V(tmp).*(Cv0 + Cv1.*L(tmp))).^2;
    ldotdot = -(2*V.*av2)./(V + bv);
    ldotdot(tmp) = -(2.*Cv1.^2.*V(tmp).^2.*(V(tmp) - Vmax))./(Vmax + V(tmp).*(Cv0 + Cv1.*L(tmp))).^3;    
    vdot = -(av2.*L.^2 + av1.*L + av0)./(V + bv) ...
        - (bv - V.*(av2.*L.^2 + av1.*L + av0))./(V + bv).^2;
    vdot(tmp) = ((V(tmp) - Vmax).*(Cv0 + Cv1.*L(tmp)))./(Vmax ...
        + V(tmp).*(Cv0 + Cv1.*L(tmp))).^2 - 1./(Vmax + V(tmp).*(Cv0 + Cv1.*L(tmp)));
    vdotdot = (2.*(av2.*L.^2 + av1.*L + av0))./(V + bv).^2 ...
        + (2.*(bv - V.*(av2.*L.^2 + av1.*L + av0)))./(V + bv).^3;
    vdotdot(tmp) = (2*(Cv0 + Cv1.*L(tmp)))./(Vmax + V(tmp).*(Cv0 + Cv1.*L(tmp))).^2 ...
        - (2.*(V(tmp) - Vmax).*(Cv0 + Cv1.*L(tmp)).^2)./(Vmax + V(tmp).*(Cv0 + Cv1.*L(tmp))).^3;
else
    ldot = NaN;
    vdot = NaN;
    ldotdot = NaN;
    vdotdot = NaN;
end
end

function [out, ldot, ldotdot] = FPE1(L)
c1 = 25.6;
k1 = 0.059;
Lr1 = 1.54;

out = c1.*k1.*log(exp((L-Lr1)./k1) + 1);

if  nargout>1
    ldot = (c1.*exp((L - Lr1)./k1))./(exp((L - Lr1)./k1) + 1);
    ldotdot = (c1.*exp((L - Lr1)./k1))./(k1.*(exp((L - Lr1)./k1) + 1)) ...
        - (c1.*exp((2*(L - Lr1))./k1))./(k1.*(exp((L - Lr1)./k1) + 1).^2);
else
    ldot = NaN;
    ldotdot = NaN;
end
end

function [out, ldot, ldotdot] = FPE2(L)
c2 = -0.020;
k2 = -18.7;
Lr2 = 0.79;

out = c2*(exp(k2*(L - Lr2)) - 1);
out(out>0) = 0;

if  nargout>1
    ldot = c2*k2*exp(k2*(L - Lr2));
    ldotdot = c2*k2^2*exp(k2*(L - Lr2));
else
    ldot = NaN;
    ldotdot = NaN;
end
end

function dot = FixNan(dot)
dot(isnan(dot)) = 0;
end

%% Muscle activation dynamics
function [a, time] = muscle_dyn(u, t_act, t_deact, dt)
szU = size(u, 1);
N =size(u, 2);
a = zeros(szU, N);
time = 0:dt:(dt*(N-1));

for k=1:N-1
    t = repmat(t_deact, [szU, 1]);
    tmp = u(:, k) > a(:, k);
    t(tmp) =  t_deact + (t_act-t_deact)*u(tmp,k);
    
    a(:,k+1) = a(:,k) + dt.*(u(:,k) - a(:,k))./t;
end
end