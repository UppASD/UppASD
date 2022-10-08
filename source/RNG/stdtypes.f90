module stdtypes
   !---------------------------------------------------------------------
   ! From the Algorithmic Conjurings of Scott Robert Ladd comes...
   !---------------------------------------------------------------------
   !
   !   stdtypes.f90 (a Fortran 95 module)
   !
   !   Definitions of common and standard integer and real types used in
   !   32- and 64-bit architectures.
   !---------------------------------------------------------------------
   !
   !  COPYRIGHT NOTICE, DISCLAIMER, and LICENSE:
   !
   !  This notice applies *only* to this specific expression of this
   !  algorithm, and does not imply ownership or invention of the
   !  implemented algorithm.
   !  
   !  If you modify this file, you may insert additional notices
   !  immediately following this sentence.
   !  
   !  Copyright 2001, 2002, 2004 Scott Robert Ladd.
   !  All rights reserved, except as noted herein.
   !
   !  This computer program source file is supplied "AS IS". Scott Robert
   !  Ladd (hereinafter referred to as "Author") disclaims all warranties,
   !  expressed or implied, including, without limitation, the warranties
   !  of merchantability and of fitness for any purpose. The Author
   !  assumes no liability for direct, indirect, incidental, special,
   !  exemplary, or consequential damages, which may result from the use
   !  of this software, even if advised of the possibility of such damage.
   !  
   !  The Author hereby grants anyone permission to use, copy, modify, and
   !  distribute this source code, or portions hereof, for any purpose,
   !  without fee, subject to the following restrictions:
   !  
   !      1. The origin of this source code must not be misrepresented.
   !  
   !      2. Altered versions must be plainly marked as such and must not
   !         be misrepresented as being the original source.
   !  
   !      3. This Copyright notice may not be removed or altered from any
   !         source or altered source distribution.
   !  
   !  The Author specifically permits (without fee) and encourages the use
   !  of this source code for entertainment, education, or decoration. If
   !  you use this source code in a product, acknowledgment is not required
   !  but would be appreciated.
   !  
   !  Acknowledgement:
   !      This license is based on the wonderful simple license that
   !      accompanies libpng.
   !
   !-----------------------------------------------------------------------
   !
   !  For more information on this software package, please visit
   !  Scott's web site, Coyote Gulch Productions, at:
   !
   !      http://www.coyotegulch.com
   !  
   !-----------------------------------------------------------------------

   ! Kind types for 64-, 32-, 16-, and 8-bit signed integers
   integer, parameter :: INT64 = selected_int_kind(18)
   integer, parameter :: INT32 = selected_int_kind(9)
   integer, parameter :: INT16 = selected_int_kind(4)
   integer, parameter :: INT08 = selected_int_kind(2)

   ! Kind types for IEEE 754/IEC 60559 single- and double-precision reals
   integer, parameter :: IEEE32 = selected_real_kind(  6,  37 )
   integer, parameter :: IEEE64 = selected_real_kind( 15, 307 )

end module stdtypes
