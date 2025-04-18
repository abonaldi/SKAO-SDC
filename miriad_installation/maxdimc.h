/*----------------------------------------------------------------------------
*  maxdimc.h - C include file containing MIRIAD-wide parameters.
*
*  Note that MAXBUF defined here is only used by the routines in xyzio.c
*  which are only used by a handful of tasks.  It is NOT the same as the
*  MAXBUF defined for Fortran use in maxdim.h and only needs to be large
*  enough for efficient I/O.  xyzio always allocates two buffers, one int
*  and the other float, using malloc() for up to the maximum number of
*  elements specified by MAXBUF.  The default value allows 64MiB each for
*  4-byte ints and floats each should be plenty.
*
*  The remaining parameters must match those in maxdim.h, or maxnax.h for
*  MAXNAX.
*
*  MAXDIM is only used by the xyio.c routines.
*
*  MAXNAX is only used by the xyio.c, xyzio.c and xyziowrap.c routines.
*
*  MAXANT and MAXWIN are only used by the uvio.c routines.
*
*  MAXBASE, MAXCHAN, and MAXWIDE are not used by anything.
*
* $Id: maxdimc.h,v 1.11 2020/05/18 05:34:57 mirmgr Exp $
*---------------------------------------------------------------------------*/
/* Maximum buffer size for xyzio (N.B. see comment above). */
#define MAXBUF  16777216

/* Maximum number of pixels on an image axis. */
#define MAXDIM  65536

/* Maximum number of image axes. */
#define MAXNAX  7

/* Maximum number of antennas. */
#define MAXANT  1024

/* Maximum number of baselines, including autocorrelations. */
#define MAXBASE ((MAXANT * (MAXANT + 1)) / 2)

/* Maximum number of channels in spectral data. */
#define MAXCHAN 70000

/* Maximum number of windows in uv data. */
#define MAXWIN  48

/* Maximum number of wideband channels. */
#define MAXWIDE 18

/* Maximum number of mosaic pointings. */
#define MAXPNT 20000

