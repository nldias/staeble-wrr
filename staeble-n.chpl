config const describe = false;
const doc = "\
=================================================================================\
==> staeble-n: my world domination evaporation model. This version specifically  \
    for Lake Mead                                                                \
                                                                                 \
./staeble-n --choice a|b|ab                                                      \
                                                                                 \
where:                                                                           \
                                                                                 \
a -- LE = rho cp A (e0 - ea)/gamma,                                              \
b -- LE = rho cp B u(e0 - ea)/gamma,                                             \
ab -- LE = rho cp (A + Bu)(e0-ea)/(2gamma)                                       \
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
// -----------------------------------------------------------------------------
// /============================================================================
// modules and initializations
// /============================================================================
use IO, Time;                  // system modules
use dgrow;
use sunearth only ddse, rsds;
use atmgas;
use angles;
use evap;
use ssr only sum;
IniPar(374.6,0.97);                // Lake Meads's altitude, water emissivity
var rlat = dec2rad(36.146084);     // Lake Mead's latitude in radians
Prescott(a=0.3,b=0.575);           // so say Lake Mead's measurements
Brutsaert(a=0.7140,b=0.0687);
// -----------------------------------------------------------------------------
// for simplicity, rho, cp, gamma and L0 will all be held constant
// -----------------------------------------------------------------------------
const L0 = 2.464e6 ;
const cp = 1005.0;
const gamma = cp*Ph/(0.622*L0);
const rho = rho_air(Ph,Th,0.0);
// /============================================================================
// declare variables that will be used (mostly) throughout
// /============================================================================
config const choice: string;
select choice {
  when "a","b","ab" {
    ;
  }
  otherwise {
    halt("choice may be one of a, b or ab");
  }
}
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
HH,              // sensible heat fluxe
LE,              // latent heat flux
DD,              // rate of change of enthalpy
XX,              // cumulative mass and heat transfer terms
YY:              // cumulative mass and heat transfer terms
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
// -----------------------------------------------------------------------------
// 1st pass: calculate net irradiance and accumulate Delta T, Delta e
// -----------------------------------------------------------------------------
for i in 1..n do {
  var (delta,rr) = ddse(yy[i],mm[i],dd[i]);
  var (Rsea,dsmax) = rsds(rlat,rr,delta);
  var S = SPrescott(Rsea,Rs[i]);
  var alb: real;
  // -----------------------------------------------------------------------------
  // check for the albedo of ice
  // -----------------------------------------------------------------------------   
  if T0[i] > 0.0 then {
    alb = WaterAlbedo(yy[i],mm[i],dd[i],rlat);
  }
  else {
    alb = 0.45;
  }
  Radiation(alb,ea[i],Ta[i]+273.15,T0[i]+273.15,S,Rs[i],Ra[i],Re[i],Rn[i]);
  // -----------------------------------------------------------------------------
  // accumulate intermediate results
  // -----------------------------------------------------------------------------
  XX[i] = rho*cp*((T0[i] - Ta[i]) + (e0[i] - ea[i])/gamma);
  YY[i] = uu[i]*XX[i];
}
// -----------------------------------------------------------------------------
// 2nd pass: calibrate the transfer coefficients A and B
// -----------------------------------------------------------------------------
var A,B: real;
var SRn = sum(Rn[1..n]);
var SXX = sum(XX[1..n]);
var SYY = sum(YY[1..n]);
A = SRn/SXX;
B = SRn/SYY;
writef("A = %8.6dr\n",A);
writef("B = %8.6dr\n",B);
// -----------------------------------------------------------------------------
// mean overall a and b in "hydrological" units
// -----------------------------------------------------------------------------
const Deltat_d = 24*3600.0;   // number of seconds in a day
var aa = Deltat_d*rho*cp*A/(gamma*L0);
var bb = Deltat_d*rho*cp*B/(gamma*L0);
// -----------------------------------------------------------------------------
// 3rd and last pass: with A and/or b adjusted, calculate 1st values of H, LE
// and D
// -----------------------------------------------------------------------------
for i in 1..n do {
  select choice {
    when "a" {
      HH[i] = rho*cp*A*(T0[i] - Ta[i]);        
      LE[i] = rho*cp*A*(e0[i] - ea[i])/gamma;
    }
    when "b" {
      HH[i] = rho*cp*B*uu[i]*(T0[i] - Ta[i]);        
      LE[i] = rho*cp*B*uu[i]*(e0[i] - ea[i])/gamma;
    }
    when "ab" {
      HH[i] = rho*cp*((A + B*uu[i])/2.0)*(T0[i] - Ta[i]);        
      LE[i] = rho*cp*((A + B*uu[i])/2.0)*(e0[i] - ea[i])/gamma;
    }
  }
  DD[i] = Rn[i] - HH[i] - LE[i];
}
// -----------------------------------------------------------------------------
// --> print results for daily data
// -----------------------------------------------------------------------------
var cod = openwriter("lakemead-"+choice+".out");
const dheader: string =
      "#     date   Rs(W/m2)  Ra(W/m2)  Re(W/m2)  Rn(W/m2)  D (W/m2)    H(W/m2)   LE(W/m2) uu(m/s)    T0(oC)    Ta(oC)    e0(Pa)    ea(Pa)\n";
cod.writef(dheader);
for i in 1..n do {
  cod.writef("%s "+"%10.2dr"*12+"\n",
             sdate[i],Rs[i], Ra[i], Re[i], Rn[i], DD[i],
             HH[i], LE[i], uu[i], T0[i],  Ta[i],e0[i], ea[i]);
}
cod.close();
