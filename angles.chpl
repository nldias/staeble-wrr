// -----------------------------------------------------------------------------
// --> todec: converts degrees, minutes and seconds to decimals
//
//
// input:  ggg.mmss  ( degrees (gg), minutes (mm) and seconds (ss) )
// output: ggg.xxxx  ( decimal degrees )
// -----------------------------------------------------------------------------
proc todec(
   in xmsc: real
): real {
   var xint = trunc(xmsc);
   var xdec = xint;
   xmsc = (xmsc - xint)*100.0;
   xint = nearbyint(xmsc);
   xdec += xint/60.0;
   xmsc = (xmsc - xint)*100;
   xdec += xmsc/3600.0;
   return(xdec) ;
}
// -----------------------------------------------------------------------------
// rad2dec: convert from radians to decimal degrees
// -----------------------------------------------------------------------------
inline proc rad2dec(const in r: real): real {
   return r*180.0/pi;
}
inline proc dec2rad(const in r: real): real {
   return r*pi/180.0;
}

