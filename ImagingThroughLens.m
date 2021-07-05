%Matlab functions for imaging process of light passing through a ideal thin lens
%
% Tianxiang Wu 2021/07/01
% wtx@zju.edu.cn

function Output = ImagingThroughLens(Photon,Input,Output)

if Photon.muz <= 0  %Photon escape
    return;
end

Photon = MoveToLens(Photon,Input);

if HitShutter(Photon,Input)  %Photon blocked
    return;    
end

Photon = RefractOnLens(Photon,Input); %Photon refraction
Photon = MoveToDetector(Photon,Input);%photon moves to detector
Output = RecordImaging(Photon,Input,Output);%Detector processes photon information

end


%Move the photon to the lens surface.
function Photon = MoveToLens(Photon,Input)

    s = (Input.Lens_z - Photon.z) / Photon.muz;
    
    Photon.x = Photon.x + s * Photon.mux;
    Photon.y = Photon.y + s * Photon.muy;
    Photon.z = Input.Lens_z;

end

%Move the photon to the detector plane.
function Photon = MoveToDetector(Photon,Input)

s = (Input.Detector_z - Input.Lens_z)/ Photon.muz;

Photon.x = Photon.x + s * Photon.mux;
Photon.y = Photon.y + s * Photon.muy;
Photon.z = Input.Detector_z;

end

% Calculate the exit angle of the photon after passing through the lens
% make sure the photon has been moved to the lens surface 
% This function only changes cosine of the photon
function Photon = RefractOnLens(Photon,Input)

FocalPoint = FindFcalPoint(Photon,Input);
% d:distance from the reference point to the center of the lens
d = sqrt(FocalPoint.x^2 + FocalPoint.y^2 +...
    (FocalPoint.z -Input.Lens_z)^2);

Photon.mux = - FocalPoint.x / d;
Photon.muy = - FocalPoint.y / d;
Photon.muz = (Input.Lens_z - FocalPoint.z) / d;

end


% find a  reference point on the front focal plane
% FocalPoint.x,.y,.z:location of point
function FocalPoint = FindFcalPoint(Photon,Input)

mux = Photon.mux;muy = Photon.muy;muz = Photon.muz;
x = Photon.x;y = Photon.y;

FocalPoint.x = 0;
FocalPoint.y = 0;
FocalPoint.z = Input.Lens_z - Input.Lens_f;

s = Input.Lens_f / muz;
FocalPoint.x = x - mux *s;
FocalPoint.y = y - muy *s;

end

% Determine whether the photon hits the shutter
% return 1 if Photon blocked
% return 0 otherwise
function HitOrNot = HitShutter(Photon,Input)

  Dr = sqrt(Photon.x^2 + Photon.y^2);
  if Dr < Input.Lens_D / 2
      HitOrNot = 0;
  else
      HitOrNot = 1;
  end
      
end

% Rasterize the photons on the detector and...
% Generate image
function Output = RecordImaging(Photon,Input,Output)
   
ixd = round((Photon.x / Input.Detector_dx + Input.Detector_xnum)/2);
iyd = round((Photon.y / Input.Detector_dy + Input.Detector_ynum)/2);

if ixd > Input.Detector_xnum
    ixd = Input.Detector_xnum;
elseif ixd <=0
    ixd = 1;
end

if iyd >Input.Detector_ynum
    iyd = Input.Detector_ynum;
elseif iyd <= 0
    iyd = 1;
end

Output.Im(ixd,iyd) = Output.Im(ixd,iyd) + Photon.w;
end
