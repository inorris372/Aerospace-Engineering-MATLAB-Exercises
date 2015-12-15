function oev = eci2orb_gooding (cbmu, r, v)

% conversion of three-dimensional position and
% velocity vectors to classical orbital elements

% input

%  cbmu = central body gravitational constant (km**3/sec**2)
%  r    = eci position vector (kilometers)
%  v    = eci velocity vector (kilometers/second)

% output

%  oev(1) = semimajor axis (kilometers)
%  oev(2) = orbital eccentricity (non-dimensional)
%           (0 <= eccentricity < 1)
%  oev(3) = orbital inclination (radians)
%           (0 <= inclination <= pi)
%  oev(4) = argument of perigee (radians)
%           (0 <= argument of perigee <= 2 pi)
%  oev(5) = right ascension of ascending node (radians)
%           (0 <= raan <= 2 pi)
%  oev(6) = true anomaly (radians)
%           (0 <= true anomaly <= 2 pi)

% R. H. Gooding, "On Universal Elements, and Conversion
% Procedures To and From Position and Velocity",
% Celestial Mechanics 44 (1988), 283-298

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = r(1);
y = r(2);
z = r(3);

xdot = v(1);
ydot = v(2);
zdot = v(3);

xsqysq = x * x + y * y;

rsq = xsqysq + z * z;

rmag = sqrt(rsq);

vr = (x * xdot + y * ydot + z * zdot) / rmag;

hx = y * zdot - z * ydot;

hy = z * xdot - x * zdot;

hz = x * ydot - y * xdot;

hsq = hx * hx + hy * hy + hz * hz;

if (hsq == 0.0)
    
    % rectilinear orbit
    
    ei = 0.5 * pi;
    
    if (xsqysq == 0.0)
        
        % axial orbit
        
        bom = 0.0;
        
    else
        
        % general rectilinear orbit
        
        bom = atan3(y, x);
        
    end
    
    u = atan3(z, sqrt(xsqysq));
    
    vt = 0.0;
    
else
    
    % non-degenerate orbit
    
    bx = hy * z - hz * y;
    by = hz * x - hx * z;
    bz = hx * y - hy * x;
    
    hx = y * bz - z * by;
    hy = z * bx - x * bz;
    hz = x * by - y * bx;
    
    w = hx * hx + hy * hy;
    
    h = sqrt(w + hz * hz);
    
    ei = atan2(sqrt(w), hz);
    
    if (w == 0.0)
        
        % orbit in reference plane
        
        bom = 0.0;
        
        u = atan3(y * sign(hz), x);
        
    else
        
        % general orbit
        
        bom = atan3(hx, -hy);
        
        u = atan3(h * z, rsq * bz);
        
    end
    
    vt = h / (rmag * rsq);
    
end

[al, q, om, tau] = pv2els(cbmu, rmag, u, vr, vt);

% semimajor axis

oev(1) = cbmu / al;

% orbital eccentricity (non-dimensional)

oev(2) = 1.0 - q / oev(1);

% orbital inclination (radians)

oev(3) = ei;

% argument of perigee (radians)

oev(4) = om;

% raan (radians)

oev(5) = bom;

% true anomaly (radians)

oev(6) = mod(u - om, 2.0 * pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [al, q, om, tau] = pv2els (gm, r, u, vr, vt)

% two-dimensional conversion from position and velocity
% to classical orbital elements (Gooding's method)

% input

%  gm = gravitational constant
%  r  = radial distance
%  u  = angle from assumed reference direction
%  vr = radial velocity
%  vt = transverse velocity

% output

%  al  = alpha (gm / sma)
%  q   = periapsis distance
%  om  = argument of periapsis relative to reference direction
%  tau = time from periapsis

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vsq = vr * vr + vt * vt;

al = 2.0 * gm / r - vsq;

alp = abs(al);

rtal = sqrt(alp);

d = r * vr;

h = r * vt;

p = h * h;

esq1 = p * al;

es = d * rtal;

eses = es * es;

ec = r * vsq - gm;

ecec = ec * ec;

if (al > 0.0)
    
    esq = ecec + eses;
    
else
    
    esq = gm * gm - esq1;
    
end

e = sqrt(esq);

q = p / (gm + e);

if (al == 0.0)
    
    tau = d * (2.0 * q + r) / (3.0 * gm);
    
    v = 2.0 * atan2(vr, vt);
    
elseif (e == 0.0)
    
    tau = 0.0;
    
    v = 0.0;
    
else
    
    e1 = al * q;
    
    if (al > 0.0)
        
        eh = atan2(es, ec);
        
        if (gm * eh * eh / 6.0 + e1 >= 0.25 * gm)
            
            em = gm * eh - esq;
            
            ecesq = gm * ec - esq;
            
        else
            
            em = gm * emkep(e1 / gm, eh);
            
            ecesq = (esq1 * ecec - esq * eses) / (esq + gm * ec);
            
        end
        
    else
        
        eh = asinh(es / e);
        
        if (gm * eh * eh / 6.0 - e1 >= 0.25 * gm)
            
            em = es - gm * eh;
            
            ecesq = esq - gm * ec;
            
        else
            
            em = e * shmkep(-e1 / e, es/ e);
            
            ecesq = -(esq1 * ecec + esq * eses) / (esq + gm * ec);
            
        end
        
    end
    
    en = alp * rtal;
    
    tau = em / en;
    
    v = atan2(es * h * rtal, ecesq);
    
end

om = u - v;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = emkep (e1, ee)

% solves ee - (1 - e1) * sin(ee)
% when (e1, ee) is close to (1, 0)

% Gooding support function

%%%%%%%%%%%%%%%%%%%%%%%%%%

y = e1 * ee;

ee2 = -ee * ee;

term = ee * (1.0 - e1);

d = 0.0;

while (1)
    
    d = d + 2.0;
    
    term = term * ee2 / (d * (d + 1.0));
    
    y0 = y;
    
    y = y - term;
    
    if (y0 == y)
        break;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = shmkep(g1, s)

% accurate computation of s - (1 - g1) * asinh(s)
% when (g1, s) is close to (0, 0)

% Gooding support function

%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 1.0 - g1;

t = s / (1.0 + sqrt(1.0 + s * s));

tsq = t * t;

x = s * (g1 + g * tsq);

term = 2.0 * g * t;

twoi1 = 1.0;

while (1)
    
    twoi1 = twoi1 + 2.0;
    
    term = term * tsq;
    
    x0 = x;
    
    x = x - term * tsq;
    
    if (x == x0)
        
        break;
        
    end
    
end



