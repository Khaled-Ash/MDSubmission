%PA 2 Part A
function AddEllipticalAtomicArray(X,Y, X0, Y0, VX0, VY0, IniDist, Temp, Type)
global C
global x y AtomSpacing
global nAtoms
global AtomType Vx Vy Mass0 Mass1

if Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end

L = (2*X - 1) * AtomSpacing;
W = (2*Y - 1) * AtomSpacing;

%Create array of atoms
xp(1, :) = linspace(-L/2, L/2, 2*X);
yp(1, :) = linspace(-W/2, W/2, 2*Y);


numAtoms = 0;
for i = 1:2*X
    for j = 1:2*Y
        if (xp(i)/(X*AtomSpacing))^2 + (yp(j)/(Y*AtomSpacing))^2 <= 1
            numAtoms = numAtoms+1;
            x(nAtoms + numAtoms) = xp(i);
            y(nAtoms  + numAtoms) = yp(j);
        else
            i
            j
        end
    end
end

%Adds noise to the initial positions
x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * IniDist + X0;
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * IniDist + Y0;

AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;

%Adds initial velocities to atoms
if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp / Mass);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end