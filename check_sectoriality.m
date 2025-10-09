function check = check_sectoriality(A)

e = eig(A*inv(A'));
phases = angle(e)/2; 

e = eig(A+A');
if e>=-1e-4
    check = true;
    return;
end

AA = exp(1i*(pi/2-(0*pi+phases(1))))*A;
e = eig(AA+AA');
if e>=-1e-6
    check = true;
    return;
end

AA = exp(1i*(pi/2-(pi+phases(1))))*A;
e = eig(AA+AA');
if e>=-1e-6
    check = true;
    return;
end

check = false;
end
