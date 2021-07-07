%MCML initial file

%% User-modifiable data
%   SourceType = 1 :infinitely narrow photon beam perpendicularly incident
%              = 2 :convergent photon beam
%              = 3 :Divergent photon beam
SourceType = 1;
%-----------------SourceType 2 parameter-------------
BeamR = 0.8;   %Beam radius [cm]
BeamDepth = 0.2; %beam focal depth [cm]

%------------------SourceType 3 parameter-------------
BeamAngle = pi/4; % Divergent half angle


%% Validation and initialization of input data
%------------------------Tianxiang 21/07/01 -------------------------------
%Add 1 extra grid to save out-of-bounds photon
nz = nz +1;
nr = nr +1;
na = na +1;
%After calculation, redundant grid will be deleted by function 'DeleteXY' 
%in file MCMLGO.m
%--------------------------------------------------------------------------


%------------------------Tianxiang 21/07/02 -------------------------------
%Verify  input data
Warningchar{1} = 'Number of layer(nl) must be positive!';
Warningchar{2} = ['Refractive index (''n'') don''t match the number of layers! ',...
    '''n''  shoule be a 1 by ',num2str(nl+2),' matrix.'];
Warningchar{3} = 'Tissue thickness(''d'') cannot be negative!';
Warningchar{4} = ['absorption( or scattering) coefficient (''mua'' or ''mus'') ',...
    'don''t match the number of layers! ',...
    'mua(or mus)  shoule be a 1 by ',num2str(nl),' matrix.'];
Warningchar{5} = '''mua'' or ''mus'' must be positive!';

InputCheck = 0;

if nl <=0
    errordlg(Warningchar{1},'''nl'' error');

elseif length(n) ~= nl +2
    errordlg(Warningchar{2},'''n'' error');

elseif d <=0
    errordlg(Warningchar{3},'''d'' error');

elseif length(mua) ~= nl ||length(mus) ~=nl
    errordlg(Warningchar{4},'''mua(mus)'' error');

elseif find((mua<0 | mus<0))
     errordlg(Warningchar{5},'''mua(mus)'' error');
else
     InputCheck = 1;
end
%--------------------------------------------------------------------------

%% Structure used to describe tissue layer
Layer = struct('n',[],... % refractive index
               'mua',[],... % absorption coefficient
               'mus',[],... % scattering coefficient
               'g',[],... % anisotropy
               'z0',[],... % z coordinate of upper boundary
               'z1',[],... % z coordinate of lower boundary
               'cos_crit0',[],... % cosine of critical angle for upper boundary
               'cos_crit1',[]... % cosine of critical angle for lower boundary
               );
           
Layer = repmat(Layer,1,nl);
for ilayer = 1:nl
    z1 = sum(d(1:ilayer));
    z0 = z1-d(ilayer);
    n0 = n(ilayer);
    n1 = n(ilayer+1); % This is actually the current layer's n. 
                            % Remember there are (nl+2) elements in n.
    n2 = n(ilayer+2);
    if (n1 > n0)
        cos_crit0 = sqrt(1-(n0^2)/(n1^2));
    else
        cos_crit0 = 0;
    end
    if (n1 > n2)
        cos_crit1 = sqrt(1-(n2^2)/(n1^2));
    else
        cos_crit1 = 0;
    end
    Layer(ilayer).n = n1;
    Layer(ilayer).mua = mua(ilayer);
    Layer(ilayer).mus = mus(ilayer);
    Layer(ilayer).g = g(ilayer);
    Layer(ilayer).z0 = z0;
    Layer(ilayer).z1 = z1;
    Layer(ilayer).cos_crit0 = cos_crit0;
    Layer(ilayer).cos_crit1 = cos_crit1;
end       
% detector layer
Layer(nl+1).n = n(end) ;
Layer(nl+1).mua = 0;
Layer(nl+1).mus = 0;
Layer(nl+1).g = 1;
Layer(nl+1).z0 = Layer(nl).z1;
% Layer(nl+1).z1 = 1e10;

%% Simulation input parameters
Input.Photon_num = np;% number of photons to be traced.
Input.Layer_num = nl;% number of layers
Input.LayerUp_n = n(1);% Refractive index deeper in the tissue
Input.LayerDown_n = n(end);% Refractive index outside the tissue
Input.wth = 0.001; %play roulette if photon weight < Wth.
Input.dz = sum(d)/nz;% z grid separation.[cm]
Input.dr = rlim/nr;% r grid separation.[cm]
Input.da = 0.5*pi/na;% a grid separation.[cm]
Input.nz = nz;% Number of grids considered in z coordinate, array ranges from 1 to nz.
Input.nr = nr;% Number of grids considered in r coordinate, array ranges from 1 to nr.
Input.na = na;% Number of grids considered in alpha coordinate, array ranges from 1 to na.
Input.rlim = rlim;
%Tianxiang 21/07/04  Light source type settings
%
%   SourceType = 1 :infinitely narrow photon beam perpendicularly incident
%              = 2 :convergent photon beam
%              = 3 :Divergent photon beam
Input.SourceType = SourceType;

%-----------Tianxiang 21/07/04 constant for convergent beam----------------
% Depth is not true depth,refraction needs to be considered
% make sure Depth is positive.
if Input.SourceType == 2
%----------------------------Modifiable data-------------------------------  
    Input.BeamR = BeamR;   %Beam radius [cm]
    Input.BeamDepth = BeamDepth; %beam focal depth [cm]
 %-------------------------------------------------------------------------   
    
    real_depth = Input.BeamDepth ; % copy of real depth.
    LayerNo = find([Layer.z0]<Input.BeamDepth,1,'last');
    if LayerNo >=2 % Refraction needs to be considered
        loc_z = Input.BeamDepth;    
        for i =LayerNo:-1:2
            u = loc_z - Layer(i).z0;
            v = Layer(i-1).n/Layer(i).n *u;
            loc_z = v - u + loc_z;
        end
        Input.BeamDepth = loc_z ; 
    end
end
%--------------------------------------------------------------------------

%-----------Tianxiang 21/07/05  Divergent beam----------------
if Input.SourceType ==3 %
    Input.BeamAngle = BeamAngle; % Divergent half angle
end
    
    
    
%% Structure used to describe a photon packet
Photon.x = 0;Photon.y = 0;Photon.z = 0;%Cartesian coordinates.[cm]
Photon.mux =0; Photon.muy =0;Photon.muz =1;%directional cosines of a photon.
Photon.w = 1;% weight.
Photon.dead = 0;% 1 if photon is terminated.
Photon.layer = 0;%index to layer where the photon  packet resides.
Photon.s = 0;%current step size[cm]
Photon.sleft = 0;

%% Structure used to describe imaging lens
% *1 magnification for imaging 
Input.Lens_f = 4; %focal length of a imaging thin lens [cm]
Input.Lens_D = rlim * 2; %aperture size of imaging lens [cm]
Input.Lens_n = 1.5;% refractive index of lens

% ---------------------Tianxiang 21/07/03---------------------------------
%Calculate the depth of Photons after passing through different  layers

 if Input.SourceType ==2 
     starti = LayerNo;
     loc_z = real_depth;
 else
     starti = 1;
     loc_z = 0;
 end
for i =starti:nl
    u = Layer(i).z1 - loc_z;
    v = Layer(i+1).n/Layer(i).n *u;
    loc_z = u - v + loc_z;
end
%Tianxiang 21/07/03  Move the lens according to loc_z.
Input.Lens_z = 2 * Input.Lens_f + loc_z; %z coordinates of imaging lens [cm]
%--------------------------------------------------------------------------
%% Structure used to describe the detector(outside the tissue)
Input.Detector_z = Input.Lens_z + 2 * Input.Lens_f;%z coordinates of detector[cm]
Input.Detector_x = rlim;%size of detector in x [cm]
Input.Detector_y = rlim;%size of detector in y [cm]
Input.Detector_xnum = 512;%resolution of detector in x 
Input.Detector_ynum = 512;%resolution of detector in y
Input.Detector_dx = Input.Detector_x / Input.Detector_xnum; %Pixel size in x[cm]
Input.Detector_dy = Input.Detector_y / Input.Detector_ynum; %Pixel size in y[cm]

%% Structure used to describe output data
Output = struct('Rsp',[],...%specular reflectance.
    'Rd_ra',zeros(nr,na),...% 2D distribution of diffuse  reflectance. [1/(cm2 sr)]
    'Rd_r',zeros(1,nr),...%1D radial distribution of diffuse  reflectance. [1/cm2]
    'Rd_a',zeros(1,nr),...%1D angular distribution of diffuse reflectance. [1/sr]
    'Rd',[],...%total diffuse reflectance. [-]
    'A_rz',zeros(nr,nz),...% 2D probability density in turbid  media over r & z. [1/cm3]
    'A_z',zeros(1,nz),...%1D probability density over z.[1/cm]
    'A_l',zeros(1,nl),...%each layer's absorption probability. [-]
    'A',[],...%total absorption probability. [-]
    'Tt_ra',zeros(nr,na),...%2D distribution of total transmittance. [1/(cm2 sr)]
    'Tt_xy',zeros(nr+1,nr+1),...%2D distribution of total transmittance. 
    'Tt_r',zeros(1,nr),...%1D radial distribution of transmittance. [1/cm2]
    'Tt_a',zeros(1,na),...%1D angular distribution of transmittance. [1/sr]
    'Tt',[],...%total transmittance.
    'Im',zeros(Input.Detector_xnum,Input.Detector_ynum),...%imaging pattern on detector
    'E_rz',zeros(nr,nz)... %energy distribution of tissue.
);

% clear unnecessary vars.
clearvars -except Input Output Layer Photon InputCheck