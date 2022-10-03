config const describe = false;
const doc = "\
=================================================================================\
==> staeble-c: incorporates stability-dependent transfer coefficients based on   \
    Charnock's equation for momentum roughness                                   \
                                                                                 \
./staeble-c                                                                      \
                                                                                 \
=================================================================================\
";
if describe then {
  writeln(doc);
  exit(0);
}
// -----------------------------------------------------------------------------
// 2021-05-27T14:15:36 beginning the new coding of staeble
// 2022-05-17T09:42:49 renaming to staeble-n: neutral only, does not
// take into account stability effects
// 2022-05-18T09:10:21 introducing iterative search for z0
// 2022-05-18T10:04:42 worked! now going for charnock
// -----------------------------------------------------------------------------
// /============================================================================
// modules and initializations
// /============================================================================
use IO, Time;                  // system modules
use dgrow;
use sunearth;
use atmgas;
use angles;
use evap;
use nstat;
use ssr only sum;
use water;
IniPar(374.6,0.97);                // Lake Meads's altitude, water emissivity
var rlat = dec2rad(36.146084);     // Lake Mead's latitude in radians
Prescott(a=0.3,b=0.575);           // so say Lake Mead's measurements
Brutsaert(a=0.7140,b=0.0687);
// -----------------------------------------------------------------------------
// measurement levels
// -----------------------------------------------------------------------------
const za = 2;                 // ea and Ta are nominally measured at 2 m
const zb = 10;                // wind is measured at 10 m
// -----------------------------------------------------------------------------
// for simplicity, rho, cp, gamma and L0 will all be held constant
// -----------------------------------------------------------------------------
const g = 9.81;                    // the acceleration of gravity
const L0 = 2.464e6 ;               // latent heat of evaporation
const cp = 1005.0;                 // specific heat of air at const pressure
const gamma = cp*Ph/(0.622*L0);    // psychrometric constant
const rho = rho_air(Ph,Th,0.0);    // density of (dry) air
const nu = viscair(Th)/rho;        // kinematic viscosity of air
// -----------------------------------------------------------------------------
// --> Psi_tau is the momentum integral function
// -----------------------------------------------------------------------------
private proc Psi_tau(const in zeta: real): real {
  var psi: real;
  if zeta < 0.0 then {
    var b = (1 - 16*zeta)**0.25;
    psi = log((b**2 + 1.0)/2.0);
    psi += 2*log((b+1.0)/2.0);
    psi -= 2*atan(b);
    psi += half_pi;
  }
  else if zeta > 1.0 then {
    psi = -5.0;
  }
  else {
    psi = -5.0*zeta;
  }
  return psi;
}
// -----------------------------------------------------------------------------
// --> Psi_E: water vapor integral function
// -----------------------------------------------------------------------------
private proc Psi_E(const in zeta: real): real {
  var psi: real;
  if zeta < 0.0 then {
    var b = (1 - 16*zeta)**0.25;
    psi = 2*log((b**2 + 1.0)/2.0);
  }
  else if zeta > 1.0 then {
    psi = -5.0;
  }
  else {
    psi = -5.0*zeta;
  }
  return psi;
}
// -----------------------------------------------------------------------------
// --> Flux: iterative flux solution using MOST similarity functions
// -----------------------------------------------------------------------------
private proc Flux(
    in i: int,
    in aM: real,              // effective momentum roughness
    in T0: real,              // water surface temperature
    in Ta: real,              // air temperature
    in e0: real,              // water vapor pressure at the surface
    in ea: real,              // water vapor pressure in air
    in uu: real               // wind speed
): (real,real,real,real) {
  const kappa = 0.4;              // von Kármán's constant
  const q0 = 0.622*e0/Ph;         // convert to specific humidity
  const qa = 0.622*ea/Ph;         // convert to specific humidity
  const delq = q0 - qa;           // sp humidity difference
  const delT = T0 - Ta;           // temperature difference
  const Tv = (1 + 0.61*qa)*Ta;    // virtual temperature
  var zeta_a: real;               // stability at za
  var zeta_b: real;               // stability at zb
  var
    CEold,             // various iterations of the transfer coeff
  CEnew,             //
  Cm,                // partial transfer coeff for momentum
  Ce,                // partial transfer coeff for water vapor
  // & sensible heat
  ustar,             // the friction velocity
  z0,                // the effective roughness
  z0plus,            // The roughness Reynolds number
  z0E                // the scalar roughness length
  : real;
  // -----------------------------------------------------------------------------
  // always start at neutral
  // -----------------------------------------------------------------------------
  z0 = 0.0002;                         // initial value: typical over water
  Cm = kappa/(log(zb/z0));             // neutral Cm
  ustar = Cm*uu;                       // "neutral" ustar
  z0plus = ustar*z0/nu;                // roughness Reynolds number
  z0E = 7.4*z0*exp(-2.25*z0plus**0.25) ; 
  Ce = kappa/(log(za/z0E));            // neutral Ce
  CEnew = Cm*Ce;                       // first guess is anything
  CEold = 2*CEnew;                     // to ensure that we enter the while
  var kount = 0;
  // if i == 817 then writeln("-"*20);
  // if i == 817 then writeln("aM = ",aM);
  // if i == 817 then writeln("Ta = ",Ta);
  // if i == 817 then writeln("Tv = ",Tv);
  // if i == 817 then writeln("delT = ",delT);
  // if i == 817 then writeln("delq = ",delq);
  // if i == 817 then writeln("u = ",uu);
  while abs(CEnew - CEold) >= 1.0e-5 do {
    CEold = CEnew;
    var qstar = Ce*delq;
    var tstar = Ce*delT;
    var tvstar = (1 + 0.61*qa)*tstar + 0.61*Ta*qstar;
    zeta_a = -kappa*g*za*tvstar/(Tv*ustar**2);
    zeta_b = -kappa*g*zb*tvstar/(Tv*ustar**2);
    // if i == 817 then writeln("CE = ",CEnew);
    // if i == 817 then writeln("u* = ",ustar);
    // if i == 817 then writeln("z0 = ",z0);
    Ce = kappa/(log(za/z0E) - Psi_E(zeta_a));
    Cm = kappa/(log(zb/z0) - Psi_tau(zeta_b));
    ustar = Cm*uu;
    z0 = aM*ustar**2/g;          // Charnock
    z0plus = ustar*z0/nu;        // roughness Reynolds number
    z0E = z0*exp(-2.25*z0plus**0.25) ; 
    CEnew = Cm*Ce;
    //    if i == 817 then writeln("CEnew = ",CEnew);
    kount += 1;
    if kount > 100 then {
      halt("did not converge in 100 iterations");
    }
  }
  var HH = rho*cp*CEnew*uu*delT;
  var LE = rho*L0*CEnew*uu*delq;
  return(HH,LE,ustar,zeta_a);
}
// =============================================================================
// start reading data
// =============================================================================
var cin = openreader("lakemead-mod.dat");
var ddata = {1..100};
var idata: domain(string);    // an associative domain
var line: string;
var sdate: [ddata] string;    // the table of dates
var idate: [idata] int;       // the inverse table of sdate
var n: int = 0;               // count the days
var
  yy,                        // year,
mm,                        // month,
dd                         // day
: [ddata] int ;
var
  T0,              // water surface temperature, MODIS
Ta,              // air temperature, ERA5
e0,              // sat vapor pressure at T0
ea,              // water vapor pressure, ERA5
uu,              // wind speed at 10 m, ERA5
Rs,              // solar radiation, ERA5
Ra,              // estimated downwelling atmospheric radiation
Re,              // emitted radiation
Rn,              // estimated net radiation
HH,              // sensible heat flux
LE,              // latent heat flux
ustar,           // friction velocity
zeta,            // Obukhov stabilit variable
DD,              // rate of change of enthalpy
XX:              // cumulative mass and heat transfer term
[ddata] real;
// /============================================================================
// reading the data is always a pain in the butt
// /============================================================================
while cin.readLine(line) do {
  if line[0] == '#' || line[0] == '\n' then {
    continue;
  }
  line = line.strip();
  var field = line.split();
  n += 1;
  dgrow(n,ddata);
  sdate[n] = field[0];            // the date table
  idata += field[0];
  idate[field[0]] = n;               // the inverse
  yy[n] = sdate[n][0..3]: int;
  mm[n] = sdate[n][5..6]: int;
  dd[n] = sdate[n][8..9]: int;
  Ta[n] = field[1]:real;
  ea[n] = field[2]:real;
  uu[n] = field[3]:real;
  Rs[n] = field[4]:real;
  T0[n] = field[5]:real;
  e0[n] = svp(T0[n]+273.15,"Tetens");
}
writeln(n," days");
ddata = {1..n};     // resize ddata
// =============================================================================
// now we try to nail it
// =============================================================================
var
  laMl,                 // log of right roughness
laMr,                 // log of left roughness
laMc:                 // middle roughness
real;
laMl = log(2.0e-7);      // very small aM
laMr = log(20.0);        // very large aM
laMc = (laMl + laMr)/2.0;
var aM = exp(laMc);
const
  dacc = 0.01;         // we want to close energy budget within 1%
var
  dclose = 0.02;
while abs(dclose) > dacc do {
  dclose = Enclose(aM);
  writeln(laMl," ",laMr," ",laMc," aM = ",aM," dclose = ",dclose);
  if dclose < 0.0 then {
    laMr = laMc;
  }
  else {
    laMl = laMc;
  }
  laMc = (laMl + laMr)/2.0;
  aM = exp(laMc);
}
// =============================================================================
// --> print results for daily data
// =============================================================================
var cod = openwriter("lakemead-ch.out");
const dheader: string =
      "#     date   Rs(W/m2)  Ra(W/m2)  Re(W/m2)  Rn(W/m2)  D (W/m2)    H(W/m2)   LE(W/m2) uu(m/s)    T0(oC)    Ta(oC)    e0(Pa)    ea(Pa)\n";
cod.writef(dheader);
for i in 1..n do {
  cod.writef("%s "+"%10.2dr"*12+"\n",
             sdate[i],Rs[i], Ra[i], Re[i], Rn[i], DD[i],
             HH[i], LE[i], uu[i], T0[i],  Ta[i],e0[i], ea[i]);
}
cod.close();
// -----------------------------------------------------------------------------
// --> EnClose: calculate the difference sum(R) - sum(XX) as a function of z0
// -----------------------------------------------------------------------------
proc Enclose(in aM: real): real {
  // -----------------------------------------------------------------------------
  // 1st pass: calculate net irradiance and accumulate Delta T, Delta e
  // -----------------------------------------------------------------------------
  for i in 1..n do {
    var (delta,rr) = ddse(yy[i],mm[i],dd[i]);
    var (Rsea,dsmax) = rsds(rlat,rr,delta);
    var S = SPrescott(Rsea,Rs[i]);
    var alb = WaterAlbedo(yy[i],mm[i],dd[i],rlat);
    var TaK = Ta[i]+273.15;
    var T0K = T0[i]+273.15;
    Radiation(alb,ea[i],TaK,T0K,S,Rs[i],Ra[i],Re[i],Rn[i]);
    // -----------------------------------------------------------------------------
    // accumulate intermediate results
    // -----------------------------------------------------------------------------
    //      writeln("i = ",i," *** date = ",sdate[i]);
    (HH[i],LE[i],ustar[i],zeta[i]) = Flux(i,aM,T0K,TaK,e0[i],ea[i],uu[i]);
    XX[i] = HH[i] + LE[i];
    DD[i] = Rn[i] - XX[i];
  }
  var SRn = sum(Rn[1..n]);
  var SXX = sum(XX[1..n]);
  // writeln("SRn = ",SRn);
  // writeln("SXX = ",SXX);
  var dclose = (SRn - SXX)/SRn;
  return dclose;
}
