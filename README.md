# staeble-wrr


This repository contains 3 Chapel programs which implement the STAEBLE
Lake Evaporation Model.
You will need to install the Chapel programming language: see [https://chapel-lang.org/].
The model is described [here](https://doi.org/10.1002/essoar.10511612.1)
A description of each program (a "help") can be printed with

./staeble-n --describe

(and so forth)

- staeble-n.chpl

  Implements STAEBLE-A, STAEBLE-B and STAEBLE-AB

- staeble-c.chpl

  Implements STAEBLE-C

- staeble-ch.chpl

  Implements STAEBLE-CH

The input file is lakemead-mod.dat. The outputs are staeble-a.out,
staeble-b.out, staeble-ab.out, staeble-c.out and staeble-ch.out, and
will be overwritten when you compile and run the program
