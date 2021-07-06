%Main
% Load input data
run Sample_Input.m
% Initialization data
run MCMLini.m

if InputCheck
    [Photon,Input,Output] = MCMLGO(Photon,Layer,Input,Output);
end


