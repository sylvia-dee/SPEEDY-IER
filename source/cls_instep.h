C--
C--   Length of the integration and time stepping constants (common ISTEPS)

C      NMONTS = 2
C      NMONTS = 12
C      NMONTS = 24
      NMONTS = 36
      NDAYSL = 30
      NSTEPS = 36

      NSTDIA = 36*30
      NSTPPR = 6
      NSTOUT = 36*30
      IDOUT  = 2
      NMONRS = 3

      ISEASC = 1
      IYEAR0 = 1948
      IMONT0 = 1

      NSTRAD = 3
      NSTRDF = 0
      INDRDF = 1

      IALST  = 1
      IASST  = 1
      IAICE  = 1

      ISST0  = (IYEAR0-1854)*12+IMONT0

      IOBSNINO = 1
C      IOBSNINO = 0

C      IFLUXCORR = 3
C      IFLUXCORR = 2
      IFLUXCORR = 0

C--
C--   Logical flags (common LFLAG1)

      LPPRES = .true.
C      LCLSTR = .false.
      LCO2 = .false.


