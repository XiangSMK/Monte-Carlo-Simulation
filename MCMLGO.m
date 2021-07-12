% Matlab Functions for Monte Carlo Simulation 
% 

% Tianxiang Wu 2021/06/01 
% wtx@zju.edu.cn

function [Photon,Input,Output] = MCMLGO(Photon,Layer,Input,Output)

    [Photon,Output] = DoOneRun(Photon,Input,Output,Layer);
    
    [Input,Output] = SumScaleResult(Layer,Input,Output);
    
end

%------
function Photon = LaunchPhoton(Photon,Layer,Input) 
Photon.w = 1 ; %
Photon.dead = 0;
Photon.layer = 1;% in tissue

if Input.SourceType == 2  % convergent photon beam
    [Photon] = ConvergBeam(Input.BeamR,Input.BeamDepth,Photon);
elseif Input.SourceType == 3 % Divergent photon beam
    [Photon] = RandomDirect(Photon,Input);
else    % normal incident
    Photon.x = 0;
    Photon.y = 0;
    Photon.z = 0;%Cartesian coordinates.[mm]
    Photon.mux =0;
    Photon.muy =0;
    Photon.muz =1;%directional cosines of a photon.
end



if (Layer(1).mua == 0 && Layer(1).mus ==0) %glass layer
    Photon.layer = 2;
    Photon.z = Layer(2).z0;
%------Tianxiang 21/6/30 Consider non-normally incident photons
    s = (Layer(2).z0 - Layer(1).z0)/Photon.muz; %Attention: muz !=0
    Photon.x = Photon.x +s * Photon.mux;
    Photon.y = Photon.x +s * Photon.muy;
end

end

% Moving photons and change its cosine to focus the beam
function [Photon] = ConvergBeam(R,Depth,Photon)

    rdn =   abs(randn()/2) ;
    while rdn >= 1
       rdn = abs(randn()/2); 
    end
    
    r = R*rdn; %(0,R) approximate normal distribution
    phi = 2*pi*rand;
    x = r*cos(phi);
    y = r*sin(phi);

    s = sqrt(x^2 + y^2 + Depth^2);
    Photon.x = x;
    Photon.y = y;
    Photon.z = 0;
    Photon.mux = -x/s;
    Photon.muy = -y/s;
    Photon.muz = Depth/s;
    
end



% ***********************************************************
%  	Choose (sample) a new theta angle for photon propagation
%  	according to the anisotropy.
%  
%  	If anisotropy g is 0, then
%  		cos(theta) = 2*rand-1.
%  	otherwise
%  		sample according to the Henyey-Greenstein function.
%  
%  	Returns the cosine of the polar deflection angle theta.
%  
function  cost = SpinTheta(g)
    if g == 0
        cost = 2*rand()-1;
    else
        temp = (1-g*g)/(1-g+2*g*rand());
        cost = (1+g*g - temp*temp)/(2*g);
        if cost < -1
            cost = -1;
        elseif cost >1
            cost = 1;
        end
        
    end
end

% ***********************************************************
%  	Choose a new direction for photon propagation by 
%  	sampling the polar deflection angle theta and the 
%  	azimuthal angle psi.
% 
%  	Note:
%    	theta: 0 - pi so sin(theta) is always positive 
%    	feel free to use sqrt() for cos(theta).
%   
%    	phi:   0 - 2pi 
%   	for 0-pi  sin(psi) is + 
%    	for pi-2pi sin(psi) is - 
%  
function Photon = Spin(g,Photon)
mux = Photon.mux;muy = Photon.muy;muz = Photon.muz;

cost = SpinTheta(g);
sint = sqrt(1 - cost^2);	%sqrt() is faster than sin()

psi = 2 * pi * rand;% spin psi 0-2pi
cosp = cos(psi);
sinp = sin(psi);

if abs(muz)>0.99999 %normal incident
    Photon.mux = sint*cosp;
    Photon.muy = sint*sinp;
    Photon.muz = cost*sign(muz);
else
    d = sqrt(1-muz^2);
    Photon.mux = (sint * (((mux*muz*cosp)-(muy*sinp))/d) ) + (mux*cost);
    Photon.muy = (sint * (((muy*muz*cosp)+(mux*sinp))/d) ) + (muy*cost);
    Photon.muz = -sint*cosp*d + muz*cost;
end

end

% Move the photon s away in the current layer of medium.  
function Photon = Hop(Photon)
s = Photon.s;
Photon.x = Photon.x + s* Photon.mux;
Photon.y = Photon.y + s* Photon.muy;
Photon.z = Photon.z + s* Photon.muz;
end

% ***********************************************************
% 	If muz ~= 0, return the photon step size in glass, 
% 	Otherwise, return 0.
% 
% 	The step size is the distance between the current 
% 	position and the boundary in the photon direction.
%  
% 	Make sure uz !=0 before calling this function.
%  
function Photon = StepSizeInGlass(Photon,Layer)
muz = Photon.muz;
%dl_b :step size to boundary
if Photon.layer ==0 %Tianxiang 21/06/30 photon in layer 0.
    if muz <=0
        Photon.dead = 1;
    elseif muz >0
        dl_b = (Layer(1).z0 - Photon.z)/muz;
    end
else
    if muz >0
        dl_b = (Layer(Photon.layer).z1 - Photon.z)/muz;
    elseif muz <0
        dl_b = (Layer(Photon.layer).z0 - Photon.z)/muz;
    else
        dl_b = 0;
    end
end

    Photon.s = dl_b;
end

% ***********************************************************
%  	Pick a step size for a photon packet when it is in 
%  	tissue.
%  	If the member sleft is zero, make a new step size 
%  	with: -log(rnd)/(mua+mus).
%  	Otherwise, pick up the leftover in sleft.
%  
%  	Layer is the index to layer.
%  	In_Ptr is the input parameters.
%  

function Photon=StepSizeInTissue(Photon,Layer)
mua = Layer(Photon.layer).mua;
mus = Layer(Photon.layer).mus;

if Photon.sleft == 0 % make a new step
   Photon.s = -log(rand)/(mua+mus);
else    %take the leftover.
    Photon.s = Photon.sleft/(mua+mus);
    Photon.sleft = 0;
end
end


% ***********************************************************
%  	Check if the step will hit the boundary.
%  	Return 1 if hit boundary.
%  	Return 0 otherwise.
%  
%   If the projected step hits the boundary, the members
%  	s and sleft of Photon are updated.
  
function [HitOrNot,Photon] = HitBoundary(Photon,Layer)

muz = Photon.muz;
if Photon.layer ==0
    mua = 0;
    mus = 0;
else
mua = Layer(Photon.layer).mua;
mus = Layer(Photon.layer).mus;
end
%distance to the boundary
if muz>0
    dl_b = (Layer(Photon.layer).z1-Photon.z)/muz;
elseif muz<0
    dl_b = (Layer(Photon.layer).z0-Photon.z)/muz;
end

if (muz ~=0 && Photon.s >dl_b) %not horizontal & crossing.
    Photon.sleft = (Photon.s - dl_b)*(mua+mus);
    Photon.s = dl_b;
    HitOrNot = 1;
else
    HitOrNot = 0;
end

end

%  ***********************************************************
%  *	Drop photon weight inside the tissue (not glass).
%  *
%  *  The photon is assumed not dead. 
%  *
%  *	The weight drop is dw = w*mua/(mua+mus).
%  *
%  *	The dropped weight is assigned to the absorption array 
%  *	elements.
function [Photon,Output] = Drop(Photon,Layer,Input,Output)
x = Photon.x;
y = Photon.y;
z = Photon.z;
%compute array indices.
izd = floor(z/Input.dz);
if izd >= Input.nz
    iz = Input.nz ;
else
    iz = izd + 1;
end

ird = floor(sqrt(x^2 +y^2)/Input.dr);
if ird >=Input.nr 
    ir = Input.nr;
else
    ir = ird +1;
end

%update photon weight.
mua = Layer(Photon.layer).mua;
mus = Layer(Photon.layer).mus;
w = Photon.w;
dwa = w *mua/(mua+mus);
Photon.w = w - dwa;
% assign dwa to the absorption array element.
Output.A_rz(ir,iz) = Output.A_rz(ir,iz)+ dwa;
Output.E_rz(ir,iz) = Output.E_rz(ir,iz)+ Photon.w;
end


%The photon weight is small, and the photon packet tries to survive a roulette.
function Photon = Roulette(Photon)
if Photon.w <=0
    Photon.dead = 1;
elseif rand()<0.1
    Photon.w = Photon.w * 10;
else
    Photon.dead = 1;
end
end

%  Compute the Fresnel reflectance.
%  Make sure that the cosine of the incident angle a1 is positive, 
%and the case when the angle is greaterthan the critical angle is ruled out.
%  Avoid trigonometric function operations as much as possible, 
%because they are computation-intensive.
% n1  incident refractive index.
% n2  transmit refractive index.
% ca1 cosine of the incident  angle. 0<a1<90 degrees. 
% ca2 cosine of the transmission angle. a2>0. 
function [r,ca2] = RFresnel(n1,n2,ca1)
if n1==n2
    ca2 = ca1;
    r = 0;
elseif ca1 > 0.99999  % normal incident.
    ca2 = ca1;
    r = (n2-n1)/(n2+n1);
    r = r^2;
elseif ca1<0.00001 % very slant
     ca2 = 0;
     r = 1;
else %general.
    sa1 = sqrt(1-ca1*ca1);
    sa2 = n1*sa1/n2;
    %sine of the incident and transmission angles
    if sa2 >1 %double check for total internal reflection.
        ca2 = 0;
        r = 1;  
    else
        ca2 = sqrt(1-sa2^2);
        
      cap = ca1*ca2 - sa1*sa2; % c+ = cc - ss. 
      cam = ca1*ca2 + sa1*sa2; % c- = cc + ss. 
      sap = sa1*ca2 + ca1*sa2; % s+ = sc + cs. 
      sam = sa1*ca2 - ca1*sa2; % s- = sc - cs. 
      
      r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam); 
		% rearranged for speed. 
    end
end
end

% ***********************************************************
%  	Record the photon weight exiting the first layer(uz<0), 
%  	no matter whether the layer is glass or not, to the 
%  	reflection array.
%  
%  	Update the photon weight as well.

function [Photon,Output]= RecordR(refl,Photon,Input,Output)
x = Photon.x;y = Photon.y;
ird = floor(sqrt(x*x+y*y)/Input.dr);
if ird >=Input.nr
    ir = Input.nr;
else
    ir = ird+1;
end
iad = floor(acos(-Photon.muz)/Input.da);
if iad>=Input.na
    ia = Input.na;
else
    ia = iad+1;
end
Output.Rd_ra(ir,ia) = Output.Rd_ra(ir,ia) +Photon.w*(1-refl);
Photon.w = Photon.w *refl;
end

% ***********************************************************
%  	Record the photon weight exiting the last layer(uz>0), 
%  	no matter whether the layer is glass or not, to the 
%  	transmittance array.
%  
%  	Update the photon weight as well.

function [Photon,Output]=RecordT(refl,Photon,Input,Output)
x = Photon.x;
y = Photon.y;
ird = floor(sqrt(x*x+y*y)/Input.dr);
if ird >=Input.nr
    ir = Input.nr;
else
    ir = ird+1;
end

iad = floor(acos(Photon.muz)/Input.da);
if iad >=Input.na
    ia = Input.na;
else
    ia = iad+1;
end
Output.Tt_ra(ir,ia) = Output.Tt_ra(ir,ia) +Photon.w*(1-refl);

%Tianxiang 21/06/30 Transmission imaging needs to keep photons alive
% Photon.w = Photon.w *refl; 

end

function Output = RecordXY(refl,Photon,Input,Output)
x = Photon.x;y = Photon.y;
ixd = round( (x/Input.dr+Input.nr)/2);
iyd = round( (y/Input.dr +Input.nr)/2);

if ixd >= Input.nr
    ixd = Input.nr;
elseif ixd <=0
    ixd = 1;
end

if iyd >= Input.nr
    iyd = Input.nr;
elseif iyd <=0
    iyd = 1;
end

% try
Output.Tt_xy(ixd,iyd) = Output.Tt_xy(ixd,iyd) + Photon.w*(1-refl);
% catch
%     disp(ixd),disp(iyd)
% end
end
% ***********************************************************
%  	Decide whether the photon will be transmitted or 
%  	reflected on the upper boundary (uz<0) of the current 
%  	layer.
%  
%  	If "layer" is the first layer, the photon packet will 
%  	be partially transmitted and partially reflected if 
%  	PARTIALREFLECTION is set to 1,
%  	or the photon packet will be either transmitted or 
%  	reflected determined statistically if PARTIALREFLECTION 
%  	is set to 0.
%  
%  	Record the transmitted photon weight as reflection.  
%  
%  	If the "layer" is not the first layer and the photon 
%  	packet is transmitted, move the photon to "layer-1".
%  
%  	Update the photon parmameters.
%  
function [Photon,Output]=CrossUpOrNot(Photon,Input,Output,Layer)
muz = Photon.muz;
if Photon.layer == 0 %Generally, photons in layer 0 will not enter this function.
    ni = Input.LayerUp_n; % Tianxiang 21/6/30 make sure photon in No.0 layer move downward.
else
    ni = Layer(Photon.layer).n;
end

if Photon.layer == 1 % Photon.layer-1 =0,cannot be the index of 'Layer'
    nt = Input.LayerUp_n;
else
    nt = Layer(Photon.layer-1).n;
end


if -muz <= Layer(Photon.layer).cos_crit0
    r = 1;
else
    [r,muz1] = RFresnel(ni,nt,-muz);
end


%transmitted to layer-1.
if rand() >r 
    if Photon.layer == 1 % Photons escape from tissue(layer 1)
        Photon.muz = -muz1;
        [Photon,Output]=RecordR(0,Photon,Input,Output); %record reflection
        Photon.dead = 1;
    else
        Photon.layer = Photon.layer -1;
        Photon.mux = Photon.mux * ni/nt;
        Photon.muy = Photon.muy * ni/nt;
        Photon.muz = -muz1;
    end
else %reflected.
    Photon.muz = - muz;
end
end

% ***********************************************************
%  	Decide whether the photon will be transmitted  or be 
%  	reflected on the bottom boundary (muz>0) of the current 
%  	layer.
%  
%  	If the photon is transmitted, move the photon to 
%  	"layer+1". If "layer" is the last layer, record the 
%  	transmitted weight as transmittance. See comments for 
%  	CrossUpOrNot.
%  
%  	Update the photon parmameters.
function [Photon,Output]=CrossDnOrNot(Photon,Input,Output,Layer)
muz = Photon.muz;

if Photon.layer==0 % Tianxiang 21/6/30 Photon.layer=0, cannot be the index of Layer.
    ni = Input.LayerUp_n;
else
    ni = Layer(Photon.layer).n;
end

 nt = Layer(Photon.layer+1).n;


if Photon.layer == 0
   [r,muz1] = RFresnel(ni,nt,muz);
else
    if muz <= Layer(Photon.layer).cos_crit1
        r = 1.0;
    else
        [r,muz1] = RFresnel(ni,nt,muz);
    end
end
if rand() >r  %transmitted to layer+1.
    if Photon.layer == Input.Layer_num %last layer
       Photon.muz = muz1; 
       Photon.mux = Photon.mux * ni/nt;
       Photon.muy = Photon.muy * ni/nt;
       [Photon,Output]=RecordT(0,Photon,Input,Output); 
%--------------
%code for Photon passing through the imaging lens.
    Output = RecordXY(0,Photon,Input,Output); 
    Output = ImagingThroughLens(Photon,Input,Output);
%--------------
       Photon.dead = 1;
    else
        Photon.layer = Photon.layer +1;
        Photon.mux = Photon.mux * ni/nt;
        Photon.muy = Photon.muy * ni/nt;
        Photon.muz = muz1;
    end
else %reflected.
     Photon.muz = -muz;
end
end
%-------------------------------------------------------
function [Photon,Output]=CrossOrNot(Photon,Input,Output,Layer)
if Photon.muz <0
    [Photon,Output]=CrossUpOrNot(Photon,Input,Output,Layer);
else 
    [Photon,Output]=CrossDnOrNot(Photon,Input,Output,Layer);
end
end
%------------------------------------------------------

% ***********************************************************
%  	Move the photon packet in glass layer.
%  	Horizontal photons are killed because they will
%  	never interact with tissue again.
% 
function [Photon,Output]=HopInGlass(Photon,Input,Output,Layer)
if Photon.muz ==0
    Photon.dead = 1;
else
    Photon = StepSizeInGlass(Photon,Layer);
    Photon = Hop(Photon);
    [Photon,Output]=CrossOrNot(Photon,Input,Output,Layer);
end
end

% ***********************************************************
%  	Set a step size, move the photon, drop some weight, 
%  	choose a new photon direction for propagation.  
%  
%  	When a step size is long enough for the photon to 
%  	hit an interface, this step is divided into two steps. 
%  	First, move the photon to the boundary free of 
%  	absorption or scattering, then decide whether the 
%  	photon is reflected or transmitted.
%  	Then move the photon in the current or transmission 
%  	medium with the unfinished stepsize to interaction 
%  	site.  If the unfinished stepsize is still too long, 
%  	repeat the above process.  
function [Photon,Output]= HopDropSpinInTissue(Photon,Input,Output,Layer)
    Photon=StepSizeInTissue(Photon,Layer);
   [HitOrNot,Photon] = HitBoundary(Photon,Layer);
   if HitOrNot
       Photon = Hop(Photon); % move to boundary plane.
       [Photon,Output]=CrossOrNot(Photon,Input,Output,Layer);
   else
       Photon = Hop(Photon);
       [Photon,Output] = Drop(Photon,Layer,Input,Output);
       Photon = Spin(Layer(Photon.layer).g,Photon);
   end
end

function [Photon,Output]= HopDropSpin(Photon,Input,Output,Layer)
%----
if Photon.layer==0 % Tianxiang 21/6/30 photon in layer 0 
    [Photon,Output]=HopInGlass(Photon,Input,Output,Layer);
else
if (Layer(Photon.layer).mua==0 && Layer(Photon.layer).mus==0) % glass layer
    [Photon,Output]=HopInGlass(Photon,Input,Output,Layer);
else
    [Photon,Output]= HopDropSpinInTissue(Photon,Input,Output,Layer);
end
end
if Photon.w < Input.wth && Photon.dead ~=0
    Photon = Roulette(Photon);
end
end



%*************************************************
%Execute Monte Carlo simulation for one independent run.
function [Photon,Output] = DoOneRun(Photon,Input,Output,Layer)

%-----------------------Timer parameters------------------
PhotonNum = 0;
inter1 = floor(Input.Photon_num / 10); %interval for timer
inter2 = floor(inter1/10); % Number of photons used in timer
NumFlag = -inter2;
%--------------------------------------------------------

while PhotonNum < Input.Photon_num
    
    %--------------------timer start----------------------------
    ticflag = mod(PhotonNum,inter1) ==0; 
    if ticflag
        NumFlag = PhotonNum;
        tic
    end
    %------------------------------------------------------
    
    %---------------------launch & move--------------------
    Photon = LaunchPhoton(Photon,Layer,Input);
    while Photon.dead == 0 && ...
            (Photon.x <= 2*Input.rlim && Photon.y <= 2*Input.rlim)
        [Photon,Output]= HopDropSpin(Photon,Input,Output,Layer);
    end
    %--------------------------------------------------------
    
    %------------------------Time ends-------------------------------
    if PhotonNum == NumFlag + inter2  
        t = toc;
        try
            DoneTime = t*(Input.Photon_num - PhotonNum)/inter2;
            disp(['Collecting Data...',...
                num2str(PhotonNum/Input.Photon_num *100),'% Done. Simulation will done in ',...
                num2str(DoneTime),' second(s).']);
        catch
            disp(['Collecting Data...',...
                num2str(PhotonNum/Input.Photon_num *100),'% Done.']);
        end
        
    end
    %---------------------------------------------------------------
    
    PhotonNum = PhotonNum +1;
end
    disp('Simulation done.');
end

%---------------------------------------------------------------------
% ***********************************************************
%  	Get 1D array elements by summing the 2D array elements.
%  
function Output = Sum2DRd(Input,Output)
    nr = Input.nr;
    na = Input.na;
for ir = 1:nr
    sum = 0;
    for ia = 1:na
        sum = sum + Output.Rd_ra(ir,ia);
    end
    Output.Rd_r(ir) = sum;
end

for ia = 1:na
    sum = 0;
    for ir = 1:nr
        sum = sum + Output.Rd_ra(ir,ia);
    end
    Output.Rd_a(ia) = sum;
end

sum = 0;
for ir = 1:nr
    sum = sum + Output.Rd_r(ir);
end
Output.Rd = sum;
end

% ***********************************************************
%  	Return the index to the layer according to the index
%  	to the grid line system in z direction (Iz).
%  
%  	Use the center of box.
%  
function i = IzToLayer(Iz,Layer,Input)
 i = 1;dz = Input.dz;
%Tianxiang 21/06/30 Iz +0.5 equal to rounding in C++ ,it's unnecessary in
%Matlab.
while (Iz *dz > Layer(i).z1 && i<Input.Layer_num)
    i = i +1;
end

end
% /***********************************************************
%  *	Get 1D array elements by summing the 2D array elements.
%  ****/
function  Output = Sum2DA(Layer,Input,Output)
    nr = Input.nr;
    nz = Input.nz;
for iz = 1:nz
    sum = 0;
    for ir = 1:nr
        sum = sum + Output.A_rz(ir,iz);
    end
    Output.A_z(iz) = sum;
end


sum = 0;
for iz =1:nz
   sum = sum + Output.A_z(iz);
   Output.A_l(IzToLayer(iz,Layer,Input)) =...
       Output.A_l(IzToLayer(iz,Layer,Input)) +Output.A_z(iz);
end
Output.A = sum;
end

% /***********************************************************
%  *	Get 1D array elements by summing the 2D array elements.
%  ****/
function Output = Sum2DTt(Input,Output)
    nr = Input.nr;
    na = Input.na;
    
    for ir = 1:nr
        sum = 0;
        for ia = 1:na
            sum = sum +Output.Tt_ra(ir,ia);
        end
        Output.Tt_r(ir) = sum;
    end
    
    for ia = 1:na
       sum = 0;
       for ir = 1:nr
           sum = sum + Output.Tt_ra(ir,ia);
       end
       Output.Tt_a(ia) = sum;
    end
    
 sum = 0;
 for ir = 1:nr
     sum = sum + Output.Tt_r(ir);
 end
 Output.Tt = sum;
end



%  **********************************************************
%  *	Scale Rd and Tt properly.
%  *
%  *	"a" stands for angle alpha.
%  ****
%  *	Scale Rd(r,a) and Tt(r,a) by
%  *      (area perpendicular to photon direction)
%  *		x(solid angle)x(No. of photons).
%  *	or
%  *		[2*PI*r*dr*cos(a)]x[2*PI*sin(a)*da]x[No. of photons]
%  *	or
%  *		[2*PI*PI*dr*da*r*sin(2a)]x[No. of photons]
%  ****
%  *	Scale Rd(r) and Tt(r) by
%  *		(area on the surface)x(No. of photons).
%  ****
%  *	Scale Rd(a) and Tt(a) by
%  *		(solid angle)x(No. of photons).
function [Output] = ScaleRdTt(Input,Output)
    nr = Input.nr;
    na = Input.na;
    dr = Input.dr;
    da = Input.da;
    
scale1 = 4 *pi *pi *dr *sin(da/2)*dr *Input.Photon_num;
%The factor (ir+0.5)*sin(2a) to be added.
for ir = 1:nr
    for ia = 1:na
%         if Input.SourceType == 2 
%             scale2 = 1 / sin(2.0*ia*da)*scale1; %Waiting to be modified
%         else
            scale2 = 1.0/((ir-0.5)*sin(2.0*ia*da)*scale1);
%         end
       Output.Rd_ra(ir,ia) = Output.Rd_ra(ir,ia) *scale2;
       Output.Tt_ra(ir,ia) = Output.Tt_ra(ir,ia) *scale2;
    end
end




scale1 = 2 *pi *dr^2*Input.Photon_num;
% area is 2*PI*[(ir+0.5)*dr]*dr.

for ir = 1:nr
    %-----------Tianxiang 21/07/06 -------
    %When non-normal incidence, A_rz and E_rz need not to be scaled
    %
%     if Input.SourceType == 2 
%         scale2 = 1/scale1; %Waiting to be modified
%     else
        scale2 = 1.0/((ir-0.5)*scale1);
%     end
    Output.Rd_r(ir) = Output.Rd_r(ir) * scale2;
    Output.Tt_r(ir) = Output.Tt_r(ir) * scale2;
end

scale1 = 2 *pi *da *Input.Photon_num;
% solid angle is 2*PI*sin(a)*da. sin(a) to be added. 

for ia = 1:na
%     if Input.SourceType == 2 
%         scale2 = 1/scale1; %Waiting to be modified
%     else
        scale2 = 1.0/(sin((ia-0.5)*da)*scale1);
%     end
    Output.Rd_a(ia) = Output.Rd_a(ia) * scale2;
    Output.Tt_a(ia) = Output.Tt_a(ia) * scale2;
end

 scale2 = 1.0/Input.Photon_num;
 Output.Rd =  Output.Rd * scale2;
 Output.Tt =  Output.Tt * scale2;
end

% /***********************************************************
%  *	Scale absorption arrays properly.
%  ****/
function [Output] = ScaleA(Input,Output)
    nr = Input.nr;
    nz = Input.nz;
    dr = Input.dr;
    dz = Input.dz;
    nl = Input.Layer_num;
    
    %Scale A_rz.
    scale1 = 2.0*pi*dr*dr*dz*Input.Photon_num;	
   % volume is 2*pi*ir*dr*dr*dz.
   for iz = 1:nz
       for ir = 1:nr
           %-----------Tianxiang 21/07/06 -------
           %When non-normal incidence, A_rz and E_rz need not to be scaled
           if Input.SourceType == 2 
               scale2 = scale1; %Waiting to be modified
           else
               scale2 = ((ir-0.5)*scale1);
           end
           Output.A_rz(ir,iz) = Output.A_rz(ir,iz) / scale2;%
           %Tianxiang 21/07/02  Energy distribution scaling
           Output.E_rz(ir,iz) = Output.E_rz(ir,iz) / scale2;%
       end
   end
   

   %Scale A_z.
     scale1 = 1.0/(dz*Input.Photon_num);
     for iz = 1:nz
         Output.A_z(iz) = Output.A_z(iz) * scale1;
     end
     
     scale1 = 1.0/Input.Photon_num;	
     for il = 1:nl
         Output.A_l(il) = Output.A_l(il) * scale1; 
     end
   
     Output.A = Output.A * scale1;
end


%Tianxiang 21/07/01 Delete out-of-bounds data
function [Input,Output] = DeleteXY(Input,Output)
    nr = Input.nr;
    nz = Input.nz;
    na = Input.na;
    
Output.A_rz(nr,:)=[];
Output.A_rz(:,nz)=[];
Output.A_z(nz)=[];

Output.E_rz(nr,:)=[];
Output.E_rz(:,nz)=[];

Output.Rd_r(nr)=[];
Output.Rd_a(na)=[];

Output.Rd_ra(nr,:)=[];
Output.Rd_ra(:,na)=[];

Output.Tt_r(nr)=[];
Output.Tt_a(na)=[];
Output.Tt_ra(nr,:)=[];
Output.Tt_ra(:,na)=[];

Output.Tt_xy(nr,:)=[];
Output.Tt_xy(1,:)=[];
Output.Tt_xy(:,nr)=[];
Output.Tt_xy(:,1)=[];

%Tianxiang 21/07/01 Restore the number of grids
Input.da = Input.da * Input.na/(Input.na -1);
Input.dr = Input.dr * Input.nr/(Input.nr -1);
Input.dz = Input.dz * Input.nz/(Input.nz -1);

Input.Detector_dx = Input.Detector_dx *Input.Detector_xnum/(Input.Detector_xnum-1);
Input.Detector_dy = Input.Detector_dy *Input.Detector_ynum/(Input.Detector_ynum-1);

Input.na = Input.na -1;
Input.nr = Input.nr -1;
Input.nz = Input.nz -1;

end
% ***********************************************************
%  	Sum and scale results of current run.
% 
function [Input,Output] = SumScaleResult(Layer,Input,Output)
    Output = Sum2DRd(Input,Output);
    Output = Sum2DA(Layer,Input,Output);
    Output = Sum2DTt(Input,Output);
    
    Output = ScaleRdTt(Input,Output);
    Output = ScaleA(Input,Output);
    % Tianxiang 21/07/01
    [Input,Output] = DeleteXY(Input,Output);
end
%---------------------------Tianxiang 21/07/04-----------------------------
%-----------------------Functions for creating suorce----------------------

%make sure Input.SourceType ==3
function [Photon] = RandomDirect(Photon,Input)
%Setting light angle with angle as half angle
    angle = Input.BeamAngle;
    rdn =   randn()/2 ;
    while abs(rdn) >= 1
       rdn =   randn()/2 ; 
    end
    phi1 =   angle * rdn ;%(0,angle)
    phi2 = 2 * pi * rand;
    cos_phi1 = cos(phi1);
    sin_phi1 = sqrt(1-cos_phi1^2); %Tianxiang 21/07/06 sqrt() is fast than sin()
    cos_phi2 = cos(phi2);
    sin_phi2 = sin(phi2);  
    
  Photon.mux = sin_phi1 * sin_phi2;
  Photon.muy = sin_phi1 * cos_phi2;
  Photon.muz = cos_phi1;
  
  Photon.x = 0;
  Photon.y = 0;
  Photon.z = 0;
end



% move photon to creat a source
function [Photon] = SetPhotonLoc(Photon,loc)
%SetPhotonLoc 
  Photon.x = loc.x;
  Photon.y = loc.y;
  Photon.z = loc.z;
end

%CreatShape creat a  rectangle source 
function [loc_mat] = CreatShape()
    width = 20;% total length[mm] col
    length = 0.1;%
    
    w_num = width*10; %here 10 is sample rate.
    l_num = length*10;
    
    x0=0;y0=0; %central position
    
    loc_w = linspace(-width/2,width/2,w_num);
    loc_l = linspace(-length/2,length/2,l_num);
    loc_mat=zeros(w_num * l_num,2);
    for i = 1:w_num
        for j = 1:l_num
            loc_mat((i-1)*l_num+j,1) = x0 + loc_w(i);
            loc_mat((i-1)*l_num+j,2) = y0 + loc_l(j);
        end
    end
end

