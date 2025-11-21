function [nreal, nimag] = nrange_no_plot(A, n)

% 
if nargin == 1
    n = 500;
end

T = linspace(0,2*pi,n);
wA = zeros(1,n);
k = 1;
for t=T
    At = exp(-1i*t)*A; 
    Ht = (At+At')/2;
    Kt = (At-At')/(2*1i);
    [Q,E] = eig(Ht);
    e = diag(E);
    m = max(e);
    ms = find(e==m);
    if length(ms)==1
        vtp = Q(:,ms);
        wA(k)=vtp'*A*vtp;
        k = k+1;
    else
        Qt = Q(:,ms);
        KtQ = Qt'*Kt*Q;
        [Qk,Ek] = eig(KtQ);
        e = diag(Ek);
        % Max
        m = max(e);
        ms = find(e==m);
        vtp = Qk(:,ms)'*Qt;
        wA(k)=vtp'*A*vtp;
        k = k+1;
        % Min
        m = min(e);
        ms = find(e==m);
        vtp = Qk(:,ms)'*Qt;
        wA(k)=vtp'*A*vtp;
        k = k+1;
    end    
end
  plot(wA,'x');
  nreal = real(wA);
  nimag = imag(wA);
end