*****************************************************************************
*                                                                           *
* Copyright 2004 Greg Bryan                                                 *
* Copyright 2004 Laboratory for Computational Astrophysics                  *
* Copyright 2004 Board of Trustees of the University of Illinois            *
* Copyright 2004 Regents of the University of California                    *
*                                                                           *
* This software is released under the terms of the "Enzo Public License"    *
* in the accompanying LICENSE file.                                         *
*                                                                           *
*****************************************************************************

c=======================================================================
c/////////////////////  SUBROUTINE F_ERROR  \\\\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine f_error (sourcefile, linenumber)
c
c     PRINT ERROR MESSAGE AND EXIT PROGRAM
c=======================================================================


      implicit none

      CHARACTER sourcefile*(*)
      INTEGER linenumber

      CALL fc_error (sourcefile // char(0), linenumber)

      return
      end

c=======================================================================
c/////////////////////  SUBROUTINE F_WARNING  \\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine f_warning (sourcefile, linenumber)
c
c     PRINT WARNING MESSAGE AND CONTINUE
c=======================================================================


      implicit none

      CHARACTER sourcefile*(*)
      INTEGER linenumber

      CALL fc_warning (sourcefile // char(0), linenumber)
      return
      end
