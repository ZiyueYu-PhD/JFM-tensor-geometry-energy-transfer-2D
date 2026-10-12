function [Q,x,y,Q_L,Q_C,Q_S] = fst_geo_efficiency_2D...
    (X,Y,velocitydatax,velocitydatay,R,sigma)

d0=(X(end)-X(1))/length(X);
Wp=2*d0/R;
lowpass=1;

% filtered velocity
u=sharpfilt2(velocitydatax,Wp,sigma,lowpass,0);
v=sharpfilt2(velocitydatay,Wp,sigma,lowpass,0);

% fluctuations
du = velocitydatax - u;
dv = velocitydatay - v;

% SGS stress built from the fluctuations du = u - filter(u)
t11=sharpfilt2(du.^2,Wp,sigma,lowpass,0);
t12=sharpfilt2(du.*dv,Wp,sigma,lowpass,0);
t21=t12;
t22=sharpfilt2(dv.^2,Wp,sigma,lowpass,0);

x=X(1,:);
y=Y(:,1);

dx = x(2)-x(1);
dy = y(2)-y(1);

[u_x,u_y] = gradient(u,dx,dy);
[v_x,v_y] = gradient(v,dx,dy);

s11=u_x;
s12=1/2*(u_y+v_x);
s21=s12;
s22=v_y;

% total energy flux Pi = -tau_ij * S_ij
Q = -(s11.*t11 + s12.*t12 + s21.*t21 + s22.*t22);

% return early if only Q,x,y are requested
if nargout<=3
    return
end


%% -------- Leonard / Cross / SGS decomposition --------

% Leonard
u2_bar = sharpfilt2(u.^2,Wp,sigma,lowpass,0);
uv_bar = sharpfilt2(u.*v,Wp,sigma,lowpass,0);
v2_bar = sharpfilt2(v.^2,Wp,sigma,lowpass,0);

t11_L = u2_bar - u.^2;
t12_L = uv_bar - u.*v;
t21_L = t12_L;
t22_L = v2_bar - v.^2;

% Cross
t11_C = sharpfilt2(u.*du + du.*u,Wp,sigma,lowpass,0);
t12_C = sharpfilt2(u.*dv + du.*v,Wp,sigma,lowpass,0);
t21_C = t12_C;
t22_C = sharpfilt2(v.*dv + dv.*v,Wp,sigma,lowpass,0);

% SGS (already computed)
t11_S = t11;
t12_S = t12;
t21_S = t21;
t22_S = t22;

% energy flux components
Q_L = -(s11.*t11_L + s12.*t12_L + s21.*t21_L + s22.*t22_L);
Q_C = -(s11.*t11_C + s12.*t12_C + s21.*t21_C + s22.*t22_C);
Q_S = -(s11.*t11_S + s12.*t12_S + s21.*t21_S + s22.*t22_S);

end