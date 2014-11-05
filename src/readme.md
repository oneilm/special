## Special functions repository

The following library files are described:

#### prini.f
Basic screen/file printing routines.

#### legeexps.f
This is a standard set of Legendre polynomial functions, also available
as a part of the FMM2D and FMM3D libraries on the CMCL website.

#### qlegefuns.f
Routines for compute Legendre functions of the _second_ kind, Q_n. The
important routines in this file are `hilbert_legendre()` and
`zqneval()`. The latter routine evalutes a series of Q_n's in the
complex plane, and the former evaluates the Hilbert tranform of a
function on [-1,1] which is samples at Legendre nodes. These routines
are more or less accurate to full machine precision using a combination
of up/down recurrences and re-scaling.

#### legehalf.f
TBA
