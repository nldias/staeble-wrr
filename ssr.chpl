// =============================================================================
// ssr is a set of search and sort functions
//
// Nelson Luis Dias
// 19900000 (circa)
// 20060421 (today)
// 2021-03-20T14:44:19 now this is today: Chapel!
// 2021-08-28T13:47:27 amin: the minimum of an array
// 2021-09-01T08:59:58 aminz: the minimum of an array not including zeros
// 2022-04-21T20:23:07 quickselect & friends added
// =============================================================================
use nstat only nanstat1;
// -----------------------------------------------------------------------------
// --> amin: the minimum of an array
// -----------------------------------------------------------------------------
proc amin(ref a: [] ?ta): ta where (ta == int || ta == real) {
   var b = max(ta);
   for x in a do {
      if x < b then {
         b = x;
      }
   }
   return b;
}
// -----------------------------------------------------------------------------
// --> amax: the maximum of an array
// -----------------------------------------------------------------------------
proc amax(
   ref a: [] ?ta
): ta where (ta == int || ta == real) {
   var b = min(ta);
   for x in a do {
      if x > b then {
         b = x;
      }
   }
   return b;
}
// -----------------------------------------------------------------------------
// --> aminz: the minimum of an array not including zeros
// -----------------------------------------------------------------------------
proc aminz(
   ref a: [] ?ta
): ta where (ta == int || ta == real) {
   const zero = 0.0:ta;              // convert to the base type zero
   var b = max(ta);
   for x in a do {
      if x == zero then {
         continue;
      }
      if x < b then {
         b = x;
      }
   }
   return b;
}
// -------------------------------------------------------------------
// --> heapsort: sorts an array of floats, using the heap algorithm
//
// adapted from heapsort (hpsort) in Numerical Recipes
// -------------------------------------------------------------------
proc heapsort(
   ref ax: [] ?ta
) where ax.rank == 1 {
// -------------------------------------------------------------------
// painful reindexing
// -------------------------------------------------------------------
   var n = ax.size;
   ref x = ax.reindex(1..n);
   var qq: ta;           // to be truly generic
// -------------------------------------------------------------------
// beginning of algorithm proper
// -------------------------------------------------------------------
   // if countval(NAN,x) > 0 then {
   //    halt("something very fishy still ...\n");
   // }
   var l = n/2 + 1 ;
   var ir = n ;
   while true do {
      if (l > 1) then {
         l = l-1;
         qq = x[l];
      }
      else {
         qq = x[ir];
         x[ir] = x[1];
         ir -= 1;
         if ir == 1 then {
            x[1] = qq;
            break ;
         }
      }
      var i = l;
      var j = l+l;
      while j <= ir do {
         if j < ir && x[j] < x[j+1] then {
            j += 1;
         }
         if (qq < x[j]) then {
            x[i] = x[j];
            i = j;
            j = j+j;
         }
         else {
            j=ir+1;
         }
      }
      x[i] = qq;
   }
//   writeln(x);
   return;
}


// -------------------------------------------------------------------
// --> indxsort: sorts an array of floats, using the heap algorithm,
// and sorting by index, stored in indx
//
// adapted from heapsort (hpsort) in Numerical Recipes
// -------------------------------------------------------------------
proc indxsort(
ref ax: [] ?ta,
ref aindx: [] int
) {
// -------------------------------------------------------------------
// painful reindexing
// -------------------------------------------------------------------
   assert(ax.rank == 1);
   var n = ax.size;
   ref x = ax.reindex(1..n);
   assert(aindx.rank == 1);
   ref indx = aindx.reindex(1..n);
   var qq: ta;           // to be truly generic
// -------------------------------------------------------------------
// index is started from 1 to n
// -------------------------------------------------------------------
   for j in 1..n do {
      indx[j] = j;
   }
// -------------------------------------------------------------------
// beginning of algorithm proper
// -------------------------------------------------------------------
   var indxt: int;
   var l = n/2 + 1 ;
   var ir = n ;
   while true do {
      if (l > 1) then {
         l -= 1;
         indxt = indx[l];
         qq = x[indxt];
      }
      else {
         indxt=indx[ir];
         qq = x[indxt];
         indx[ir] = indx[1];
         ir -= 1;
         if ir == 1 then {
            indx[1] = indxt ;
            break ;
         }
      }
      var i = l;
      var j = l+l;
      while j <= ir do {
         if ((j < ir) && (x[indx[j]] < x[indx[j+1]])) then {
            j += 1;
         }
         if (qq < x[indx[j]]) then {
            indx[i] = indx[j];
            i = j;
            j = j+j;
         }
         else {
            j=ir+1;
         }
      }
      indx[i] = indxt;
   }
// -------------------------------------------------------------------
// we still need to fix the indices
// -------------------------------------------------------------------
   const del = ax.indices.first - 1;
   indx += del;
   return;
}
// -------------------------------------------------------------------
// --> countval: count how many times val occurs in x
// -------------------------------------------------------------------
proc countval(val: ?tv, x: [] ?tx): int {
// --------------------------------------------------------------------
// be careful with empty arrays
// --------------------------------------------------------------------
   assert(x.rank == 1);
   assert(tx == tv);
   var n = x.size;
   if n == 0 then {
      halt("ssr-->countval: empty array");
   }
   var count = 0;
// -------------------------------------------------------------------
// be careful with NANs
// -------------------------------------------------------------------
   if tv == real && isnan(val) then {
      for e in x do {
         if isnan(e) then count += 1;
      }
   }
   else {
      for e in x do {
         if e == val then count +=1 ;
      }
   }
   return count;
}
// -------------------------------------------------------------------
// --> purgeval: purge all occurrences of val inside the array
// -------------------------------------------------------------------
proc purgeval(val: ?tv, ax: [] ?tx): [] tx {
// --------------------------------------------------------------------
// be careful with empty arrays
// --------------------------------------------------------------------
   assert(ax.rank == 1);
   assert(tx == tv);
   var n = ax.size;
   if n == 0 then {
      halt("ssr-->purgeval: empty array");
   }
   var count = countval(val, ax);
   ref x = ax.reindex(0..n-1);
   var y: [0..n-1-count] tx;       // will return this array
// -------------------------------------------------------------------
// be careful with NANs
// -------------------------------------------------------------------
   if tv == real && isnan(val) then {
      var k = 0;
      for i in 0..#n do {
         if !isnan(x[i]) then {
            y[k] = x[i];
            k += 1;
        }
      }
   }
   else {
      var k = 0;
      for i in 0..#n do {
         if x[i] != val then {
            y[k] = x[i];
            k += 1;
         }
      }
   }
   return y;
}

// -------------------------------------------------------------------
// --> sum: sum all elements in array x
// -------------------------------------------------------------------
inline proc sum(x: [] ?tx): tx {
   return (+ reduce x) ;
}
// -------------------------------------------------------------------
// --> butsum: sum all elements but those which are equal to val in
// array x. Returns the number of elements found that are different
// from val, and the overall sum.
// -------------------------------------------------------------------
proc butsum(val: ?tv, x: [] ?tx): (int, tv) {
//   assert (x.rank == 1);
   assert (tx == tv);
   var n = x.size;
   if n == 0 then {
      halt("ssr-->countval: empty array");
   }
   var nb = 0;
   var sum: tv = 0;
   if tv == real && isnan(val) then {
      for e in x do {
         if !isnan(e) then sum += e;
         nb += 1;
      }
   }
   else {
      for e in x do {
         if e != val then sum += e;
         nb += 1;
      }
   }
   return (nb,sum);
}
// -------------------------------------------------------------------
// --> whereval: where val occurs in x
//
// important: whereval *always* returns a 0-based array
// -------------------------------------------------------------------
use dgrow;
proc whereval(val, x: [?dx] ?tx): [] int {
// --------------------------------------------------------------------
// be careful with empty arrays
// --------------------------------------------------------------------
   assert(x.rank == 1);
   assert(tx == val.type);
   var n = x.size;
   if n == 0 then {
      halt("--> whereval: empty array");
   }
   const xf = dx.first;       // try to be agnostic
   const xl = dx.last;        // try to be agnostic
   const m = max(n/10,2);     // guess that 10% of x elements == val
   var dw = {0..#m};          // ... but at least 2
   var ww: [dw] int;          // where they are
   var ct = 0;                // count how many
   if val.type == real && isnan(val) then {  // be careful with NANs
      for i in xf..xl do {
         if isnan(x[i]) then {
            dgrow(ct,dw);
            ww[ct] = i;
            ct += 1;
         }
      }
   }
   else {                                    // just find them
      for i in xf..xl do {
         if x[i] == val then {
            dgrow(ct,dw);
            ww[ct] = i;
            ct += 1;
         }
      }
   }
   dw = {0..#ct};                   // adjust domain
   return ww;
}
// -------------------------------------------------------------------
// --> diff: calculates the 1st discrete difference of a 1D array
// bool--int "magic" is used!
//
// important: diff *always* returns a 0-based array
// -------------------------------------------------------------------
proc diff(ref ax: [] ?tx) {
   assert (ax.rank == 1);     // must be 1D
   var n = ax.size;           // count elements
   ref x = ax.reindex(1..n);  // reindex
// -------------------------------------------------------------------
// return an array with a domain that is compatible with ax's domain
// -------------------------------------------------------------------   
   type td;
   if tx == bool then {
      td = int;
   }
   else {
      td = tx;
   }
   var dd = {0..#(n-1)};      // always return a 0-based array
   var dx: [dd] td;           // the return array
   dx = x[2..n] - x[1..n-1];  // differentiated
   return dx;                 // end of the story
}
// -------------------------------------------------------------------
// --> linspace: my equivalent of a (simple!) numpy linspace. returns
// a 0-based 1D array with n linearly interpolated values including,
// and between, start and stop.
// -------------------------------------------------------------------
proc linspace(
   start: real,               // the first value
   stop: real,                // the last value
   n: int                     // how many do you want?
): [] real {
   var x: [0..#n] real;
   assert( n > 1);
   var dx = (stop-start)/(n-1);
   forall i in 0..#n do {     // this is a parallel algorithm!
      x[i] = start + i*dx;
   }
   return x;
}
// -------------------------------------------------------------------
// --> flip: flips a 1D array
// -------------------------------------------------------------------
proc flip(ref ax: [] ) {
   assert(ax.rank == 1);
   var n = ax.size;
   ref x = ax.reindex(0..#n);
   for i in 0..n/2 do {
      x[i] <=> x[n-1-i];
   }
}
// -----------------------------------------------------------------------------
// --> interp: searches the array x until x[i] < xc ; then linearly
// interpolates. THE ARRAY x MUST BE SORTED!
//
// 2008-05-15T09:44:20 -- version for an array of shorts
//
// 2008-05-15T11:07:07 -- implementing a much faster binary search (hopefully)
//
// 2021-04-16T16:57:14 -- translating from C to Chapel
// -----------------------------------------------------------------------------
proc interp(
   const in xc: ?tc,         // the value being searched  
   const ref ax: [] ?tx,     // the table                 
   const ref ay: [] ?ty      // the table 
) : ty {
// -----------------------------------------------------------------------------
// the painful checks
// -----------------------------------------------------------------------------
   assert(tc == tx);
   assert(ax.rank == 1);
   assert(ay.rank == 1);
   const n = ax.size;
// -------------------------------------------------------------------------------------------
// must have at least 2 points for interpolation
// -------------------------------------------------------------------------------------------
   if ( n < 2 ) then {
      halt("--> interp: x,y size must be >= 2") ;
   }
// -------------------------------------------------------------------------------------------
// error conditions: xc must be within (x[0],x[n-1])
// -------------------------------------------------------------------------------------------
   ref x = ax.reindex(0..n-1);
   ref y = ay.reindex(0..n-1);
   if  (xc < x[0]) || (xc > x[n-1])  then {
      halt("xc = ",xc," out of range:",x[0]," ",x[n-1]) ;
   }
// -------------------------------------------------------------------------------------------
// a simple binary search algol
// -------------------------------------------------------------------------------------------
   var iu = n - 1 ;
   var il = 0 ;
   var im: int;
   while iu - il > 1 do {
      im = (iu + il)/2 ;
      if xc <= x[im] { 
	 iu = im ;
      }
      else {
	 il = im ;
      }
   }
// -------------------------------------------------------------------------------------------
// linear interpolation and return
// -------------------------------------------------------------------------------------------
   var dx = (x[iu] - x[il]) ;
   var dy =  (y[iu] - y[il]) ;
   var m = dy/dx ;
   return (y[il] + m * (xc - x[il])) ;
}
// -----------------------------------------------------------------------------
// --> indx_interp: the indices around the interpolated value
//
// 2021-05-16T11:52:51 -- may be quite useful! but unfinished!
// -----------------------------------------------------------------------------
proc indx_interp(
   const in xc: ?tc,         // the value being searched  
   const ref ax: [] ?tx      // the table                 
   ) : (int,int) {
// -----------------------------------------------------------------------------
// the painful checks
// -----------------------------------------------------------------------------
   assert(tc == tx);
   assert(ax.rank == 1);
   const n = ax.size;
   const xfirst = ax.domain.first;
// -------------------------------------------------------------------------------------------
// must have at least 2 points for interpolation
// -------------------------------------------------------------------------------------------
   if ( n < 2 ) then {
      halt("--> interp: x size must be >= 2") ;
   }
// -------------------------------------------------------------------------------------------
// error conditions: xc must be within (x[0],x[n-1])
// -------------------------------------------------------------------------------------------
   ref x = ax.reindex(0..n-1);
   if  (xc < x[0]) || (xc > x[n-1])  then {
      halt("xc = ",xc," out of range:",x[0]," ",x[n-1]) ;
   }
// -------------------------------------------------------------------------------------------
// a simple binary search algol
// -------------------------------------------------------------------------------------------
   var iu = n - 1 ;
   var il = 0 ;
   var im: int;
   while iu - il > 1 do {
      im = (iu + il)/2 ;
      if xc <= x[im] { 
	 iu = im ;
      }
      else {
	 il = im ;
      }
   }
// -------------------------------------------------------------------------------------------
// the bracketing indices of the *original* (not reindexed!) array
// -------------------------------------------------------------------------------------------
   return (il+xfirst,iu+xfirst);
}
// -----------------------------------------------------------------------------
// The partition, quickselect, indxpartition and indxquickselect are an
// adaptation of two sources:
//
// 1. https://stackoverflow.com/questions/5380568/\
//    algorithm-to-find-k-smallest-numbers-in-array-of-n-items
// 2. https://en.wikipedia.org/wiki/Quickselect
//
// I did not understand the 1. well, but it helped implement 2., with
// adaptations
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// --> partition: private partition for indxquickselect
// -----------------------------------------------------------------------------
private proc partition(ref A: [] real, in left: int, in right: int): int 
   where A.rank == 1 {
   const pivot = A[right];
   var i = left;
   var x: int;
   for x in left..right-1 do { // (x = left; x < right; x++) {
      if (A[x] < pivot) then {
         A[i] <=> A[x];
         i += 1;
      }
   }
   A[i] <=> A[right];
   return i;
}


// -----------------------------------------------------------------------------
// --> quickselect: returns the kth smallest value of A. After it
//     is called with quickselect(indx,A,A.indices.first,A.indices.last,k):
//
// A[A.indices.first] is the smallest
// A[A.indices.first+1] is the second smallest
// ...
// A[A.indices.first+k-1] is the kth smallest
// -----------------------------------------------------------------------------
proc quickselect(
   ref A: [] real,       // partially sorted array
   in left: int,         // left index
   in right: int,        // right index
   in k: int             // k smallest
   ): real {
   if left == right then {
      return A[left];
   }
// -----------------------------------------------------------------------------   
// p is position of pivot in the partitioned array
// -----------------------------------------------------------------------------   
   var p = partition(A, left, right); 
   if k == p then {
      return A[k] ;
   }
   else if k < p then {
      return quickselect(A,left,p-1,k);
   }
   else {
      return quickselect(A,p+1,right,k);
   }
}

// -----------------------------------------------------------------------------
// --> indxpartition: private partition for indxquickselect
// -----------------------------------------------------------------------------
private proc indxpartition(
   ref indx: [] int,     // indexes of partially sorted array
   ref A: [] real,       // A remains intact
   in left: int,         // left index
   in right: int         // right index
   ): int 
   where A.rank == 1 {
   const pivot = A[indx[right]];
   var i = left;
   var x: int;
   for x in left..right-1 do { 
      if (A[indx[x]] < pivot) then {
         indx[i] <=> indx[x];
         i += 1;
      }
   }
   indx[i] <=> indx[right];
   return i;
}
// -----------------------------------------------------------------------------
// --> indxquickselect: returns the index of the kth smallest value. After it
//     is called with indxquickselect(indx,A,A.indices.first,A.indices.last,k):
//
// A[indx[1]] is the smallest
// A[indx[2]] is the second smallest
// ...
// A[indx[k]] is the kth smallest
// -----------------------------------------------------------------------------
proc indxquickselect(
   ref indx: [] int,     // indexes of partially sorted array
   ref A: [] real,       // A remains intact
   in left: int,         // left index
   in right: int,        // right index
   in k: int             // k smallest values
   ) : int where (A.rank == 1) {
   assert (A.shape == indx.shape);
// -----------------------------------------------------------------------------
// is this game over?
// -----------------------------------------------------------------------------
   if left == right then {
      return indx[left];
   }
// -----------------------------------------------------------------------------   
// p is position of pivot in the partitioned array
// -----------------------------------------------------------------------------   
   var p = indxpartition(indx, A, left, right); 
   if k == p then {
      return indx[k] ;
   }
   else if k < p then {
      return indxquickselect(indx,A,left,p-1,k);
   }
   else {
      return indxquickselect(indx,A,p+1,right,k);
   }
}
// ------------------------------------------------------------------------------
// --> fillgaps_li: fill gaps in array ax (flagged with NANs) by linear
// interpolation from adjacent data.
// ------------------------------------------------------------------------------
proc fillgaps_li(ref ax: [] real) where ax.rank == 1 {
   var nx = ax.size;
   ref x = ax.reindex(1..nx);
// -----------------------------------------------------------------------------
// if passed everything but there are still gaps, fills each gap with linear
// interpolation. first, find each gap:
// -----------------------------------------------------------------------------
   var db = {0..nx+1};
   var binv: [db] bool;       // true when NAN in x 
   binv[0] = false;           // place sentinels around binv
   binv[1..nx] = isnan(x);    // and find nans inside
   binv[nx+1] = false;        // place sentinels around binv
   var dinv = diff(binv);     // differentiate binv
// -----------------------------------------------------------------------------
// now +1 marks the beginning of gaps, and -1 their end, in dinv. the number of
// -1s is equal to the number of +1s. ibeg and iend are the indices of the last
// valid datum before each gap and first valid datum after each gap in 1-based
// array x
// -----------------------------------------------------------------------------
   var ibeg = whereval(1,dinv);         // note that ibeg is 0-based
   var iend = whereval(-1,dinv) + 1;    // note that iend is 0-based
   var nruns = ibeg.size;               // how many gaps?
   assert (nruns == iend.size);
// -----------------------------------------------------------------------------
// finally, linear interpolation of gaps: we place two sentinels with the mean
// at both ends of x. Since there are nans in x, we need to use nanmean
// -----------------------------------------------------------------------------
   var (nnans, xmean): 2*real = nanstat1(x);
// -----------------------------------------------------------------------------
// xcat acts like x with two sentinels at the extremeties
// -----------------------------------------------------------------------------
   var xcat: [db] real;
   xcat[0] = xmean;
   xcat[1..nx] = x;
   xcat[nx+1] = xmean;
// -----------------------------------------------------------------------------
// (parallel!) loop over gaps: linear interpolation with linspace
// -----------------------------------------------------------------------------
   forall ir in 0..#nruns do {
      var irb = ibeg[ir];          // last valid position before
      var ire = iend[ir];          // first valid position after
      var xstart = xcat[irb];      // valid datum before
      var xstop  = xcat[ire];      // valid datum after
      var gaplen = ire - irb + 1;  // gap size
      var xfill = linspace(xstart,xstop,gaplen);  // linear interp
      x[irb+1..ire-1] = xfill[1..gaplen-2];       // fill gaps
   }
   return;
}
