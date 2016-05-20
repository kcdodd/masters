function [no, vro, vzo] = sumRadialSol(ni, vri, vzi)


no = zeros(65, 6);
vro = zeros(65, 6);
vzo = zeros(65, 6);

for p = 1:65
    for r = 1:11
        no(p, r) = (ni(p, r, 1)*0.5 + ni(p, r, 2)*0.75 + ni(p, r, 3)*0.375 + ni(p, r, 4)*0.1875 + ni(p, r, 5)*0.1875)/2;
        vro(p, r) = (vri(p, r, 1)*0.5 + vri(p, r, 2)*0.75 + vri(p, r, 3)*0.375 + vri(p, r, 4)*0.1875 + vri(p, r, 5)*0.1875)/2;
        vzo(p, r) = (vzi(p, r, 1)*0.5 + vzi(p, r, 2)*0.75 + vzi(p, r, 3)*0.375 + vzi(p, r, 4)*0.1875 + vzi(p, r, 5)*0.1875)/2;
    end
end