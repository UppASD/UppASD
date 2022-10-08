!> Sorting routines
module Sorting
   use Parameters
   use Profiling

   implicit none
   public
contains
   subroutine Merge(A,NA,B,NB,C,NC)

      integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
      integer, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
      integer, intent(in)     :: B(NB)
      integer, intent(in out) :: C(NC)

      integer :: I,J,K

      I = 1; J = 1; K = 1;
      do while(I <= NA .and. J <= NB)
         if (A(I) <= B(J)) then
            C(K) = A(I)
            I = I+1
         else
            C(K) = B(J)
            J = J+1
         endif
         K = K + 1
      enddo
      do while (I <= NA)
         C(K) = A(I)
         I = I + 1
         K = K + 1
      enddo
      return

   end subroutine merge

   recursive subroutine MergeSort(A,N,T)

      integer, intent(in) :: N
      integer, dimension(N), intent(in out) :: A
      integer, dimension((N+1)/2), intent (out) :: T

      integer :: NA,NB,V

      if (N < 2) return
      if (N == 2) then
         if (A(1) > A(2)) then
            V = A(1)
            A(1) = A(2)
            A(2) = V
         endif
         return
      endif
      NA=(N+1)/2
      NB=N-NA

      call MergeSort(A,NA,T)
      call MergeSort(A(NA+1),NB,T)

      if (A(NA) > A(NA+1)) then
         T(1:NA)=A(1:NA)
         call Merge(T,NA,A(NA+1),NB,A,N)
      endif
      return

   end subroutine MergeSort

   !subroutine Merge2D(A,NA,B,NB,C,NC)
   subroutine Merge2D(A,NA,B,NB,C,NC,IA,IB,IC,IW,RA,RB,RC,RW)
      implicit none

      integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
      integer, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
      integer, intent(in)     :: B(NB)
      integer, intent(in out) :: C(NC)
      integer, intent(in) :: IW
      integer, intent(in) :: RW
      integer, dimension(IW), intent(in out) :: IA(NA,1:IW)
      integer, dimension(IW), intent(in)     :: IB(NB,1:IW)
      integer, dimension(IW), intent(in out) :: IC(NC,1:IW)
      real(dblprec), dimension(RW), intent(in out) :: RA(NA,1:RW)
      real(dblprec), dimension(RW), intent(in)     :: RB(NB,1:RW)
      real(dblprec), dimension(RW), intent(in out) :: RC(NC,1:RW)

      integer :: I,J,K

      write(*,*) "Entered Merge2D"
      I = 1; J = 1; K = 1;
      do while(I <= NA .and. J <= NB)
         write(*,*) "Before C(K) ", C
         write(*,*) "Before IC(K,1:IW) ", IC
         write(*,*) "Before RC(K,1:RW) ", RC
         if (A(I) <= B(J)) then
            C(K) = A(I)
            IC(K,1:IW) = IA(I,1:IW)
            RC(K,1:RW) = RA(I,1:RW)
            I = I+1
         else
            C(K) = B(J)
            IC(K,1:IW) = IB(J,1:IW)
            RC(K,1:RW) = RB(J,1:RW)
            J = J+1
         endif
         write(*,*) "After  C(K) ", C
         write(*,*) "After  IC(K,1:IW) ", IC
         write(*,*) "After  RC(K,1:RW) ", RC
         K = K + 1
      enddo
      write(*,*) "Middle Merge2D"
      do while (I <= NA)
         C(K) = A(I)
         IC(K,1:IW) = IA(I,1:IW)
         RC(K,1:RW) = RA(I,1:RW)
         I = I + 1
         K = K + 1
      enddo
      write(*,*) "End Merge2D"
      return

   end subroutine merge2D

   recursive subroutine MergeSort2D(A,N,T,IA,IW,IT,RA,RW,RT)

      implicit none
      integer, intent(in) :: N
      integer, dimension(N), intent(in out) :: A
      integer, dimension((N+1)/2), intent (out) :: T
      integer, intent(in) :: IW
      integer, dimension(N,IW), intent(in out) :: IA
      integer, dimension((N+1)/2,IW), intent (out) :: IT
      integer, intent(in) :: RW
      real(dblprec), dimension(N,RW), intent(in out) :: RA
      real(dblprec), dimension((N+1)/2,RW), intent (out) :: RT

      integer :: NA,NB,V
      integer, dimension(IW) :: IV
      real(dblprec), dimension(RW) :: RV

      write(*,*) "Entered MergeSort2D"
      if (N < 2) return
      if (N == 2) then
         if (A(1) > A(2)) then
            V = A(1)
            A(1) = A(2)
            A(2) = V
            IV(1:IW) = IA(1,1:IW)
            IA(1,1:IW) = IA(2,1:IW)
            IA(2,1:IW) =IV(1:IW)
            RV(1:RW) = RA(1,1:RW)
            RA(1,1:RW) = RA(2,1:RW)
            RA(2,1:RW) = RV(1:RW)
         endif
         return
      endif
      NA=(N+1)/2
      NB=N-NA
      write(*,*) "Middle MergeSort2D"

      call MergeSort2D(A,NA,T,IA,IW,IT,RA,RW,RT)
      call MergeSort2D(A(NA+1),NB,T,IA(NA+1,1:IW),IW,IT,RA(NA+1,1:RW),RW,RT)

      if (A(NA) > A(NA+1)) then
         T(1:NA)=A(1:NA)
         IT(1:NA,1:IW)=IA(1:NA,1:IW)
         RT(1:NA,1:RW)=RA(1:NA,1:RW)
         write(*,*) "Do merge"
         call Merge2D(T,NA,A(NA+1),NB,A,N,  IT,IA(NA+1,1:IW),IA,IW,  RT,RA(NA+1,1:RW),RA,RW)
      endif
      return

   end subroutine MergeSort2D

   subroutine MergeIR(A,NA,B,NB,C,NC,RA,RB,RC)
      implicit none

      integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
      integer, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
      integer, intent(in)     :: B(NB)
      integer, intent(in out) :: C(NC)
      real(dblprec), intent(in out) :: RA(NA)
      real(dblprec), intent(in)     :: RB(NB)
      real(dblprec), intent(in out) :: RC(NC)

      integer :: I,J,K

      I = 1; J = 1; K = 1;
      do while(I <= NA .and. J <= NB)
         if (A(I) <= B(J)) then
            C(K) = A(I)
            RC(K) = RA(I)
            I = I+1
         else
            C(K) = B(J)
            RC(K) = RB(J)
            J = J+1
         endif
         K = K + 1
      enddo
      do while (I <= NA)
         C(K) = A(I)
         RC(K) = RA(I)
         I = I + 1
         K = K + 1
      enddo
      return

   end subroutine mergeIR


   recursive subroutine MergeSortIR(A,N,T,RA,RT)

      implicit none
      integer, intent(in) :: N
      integer, dimension(N), intent(inout) :: A
      integer, dimension((N+1)/2), intent (out) :: T
      real(dblprec), dimension(N), intent(in out) :: RA
      real(dblprec), dimension((N+1)/2), intent (out) :: RT

      integer :: NA,NB,V
      real(dblprec) :: RV

      if (N < 2) return
      if (N == 2) then
         if (A(1) > A(2)) then
            V = A(1)
            A(1) = A(2)
            A(2) = V
            RV = RA(1)
            RA(1) = RA(2)
            RA(2) = RV
         endif
         return
      endif
      NA=(N+1)/2
      NB=N-NA

      call MergeSortIR(A,NA,T,RA,RT)
      call MergeSortIR(A(NA+1),NB,T,RA(NA+1),RT)

      if (A(NA) > A(NA+1)) then
         T(1:NA)=A(1:NA)
         RT(1:NA)=RA(1:NA)
         call MergeIR(T,NA,A(NA+1),NB,A,N,RT,RA(NA+1),RA)
      endif
      return

   end subroutine MergeSortIR

   subroutine MergeII(A,NA,B,NB,C,NC,RA,RB,RC)
      implicit none

      integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
      integer, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
      integer, intent(in)     :: B(NB)
      integer, intent(in out) :: C(NC)
      integer, intent(in out) :: RA(NA)
      integer, intent(in)     :: RB(NB)
      integer, intent(in out) :: RC(NC)

      integer :: I,J,K

      I = 1; J = 1; K = 1;
      do while(I <= NA .and. J <= NB)
         if (A(I) <= B(J)) then
            C(K) = A(I)
            RC(K) = RA(I)
            I = I+1
         else
            C(K) = B(J)
            RC(K) = RB(J)
            J = J+1
         endif
         K = K + 1
      enddo
      do while (I <= NA)
         C(K) = A(I)
         RC(K) = RA(I)
         I = I + 1
         K = K + 1
      enddo
      return

   end subroutine mergeII


   recursive subroutine MergeSortII(A,N,T,RA,RT)

      implicit none
      integer, intent(in) :: N
      integer, dimension(N), intent(in out) :: A
      integer, dimension((N+1)/2), intent (out) :: T
      integer, dimension(N), intent(in out) :: RA
      integer, dimension((N+1)/2), intent (out) :: RT

      integer :: NA,NB,V
      integer :: RV

      if (N < 2) return
      if (N == 2) then
         if (A(1) > A(2)) then
            V = A(1)
            A(1) = A(2)
            A(2) = V
            RV = RA(1)
            RA(1) = RA(2)
            RA(2) = RV
         endif
         return
      endif
      NA=(N+1)/2
      NB=N-NA

      call MergeSortII(A,NA,T,RA,RT)
      call MergeSortII(A(NA+1),NB,T,RA(NA+1),RT)

      if (A(NA) > A(NA+1)) then
         T(1:NA)=A(1:NA)
         RT(1:NA)=RA(1:NA)
         call MergeII(T,NA,A(NA+1),NB,A,N,RT,RA(NA+1),RA)
      endif
      return

   end subroutine MergeSortII

   subroutine MergeRI(A,NA,B,NB,C,NC,RA,RB,RC)
      implicit none

      integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
      real(dblprec), intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
      real(dblprec), intent(in)     :: B(NB)
      real(dblprec), intent(in out) :: C(NC)
      integer, intent(in out) :: RA(NA)
      integer, intent(in)     :: RB(NB)
      integer, intent(in out) :: RC(NC)

      integer :: I,J,K

      I = 1; J = 1; K = 1;
      do while(I <= NA .and. J <= NB)
         if (A(I) <= B(J)) then
            C(K) = A(I)
            RC(K) = RA(I)
            I = I+1
         else
            C(K) = B(J)
            RC(K) = RB(J)
            J = J+1
         endif
         K = K + 1
      enddo
      do while (I <= NA)
         C(K) = A(I)
         RC(K) = RA(I)
         I = I + 1
         K = K + 1
      enddo
      return

   end subroutine mergeRI


   recursive subroutine MergeSortRI(A,N,T,RA,RT)

      implicit none
      integer, intent(in) :: N
      real(dblprec), dimension(N), intent(in out) :: A
      real(dblprec), dimension((N+1)/2), intent (out) :: T
      integer, dimension(N), intent(in out) :: RA
      integer, dimension((N+1)/2), intent (out) :: RT

      integer :: NA,NB,RV
      real(dblprec) :: V

      if (N < 2) return
      if (N == 2) then
         if (A(1) > A(2)) then
            V = A(1)
            A(1) = A(2)
            A(2) = V
            RV = RA(1)
            RA(1) = RA(2)
            RA(2) = RV
         endif
         return
      endif
      NA=(N+1)/2
      NB=N-NA

      call MergeSortRI(A,NA,T,RA,RT)
      call MergeSortRI(A(NA+1),NB,T,RA(NA+1),RT)

      if (A(NA) > A(NA+1)) then
         T(1:NA)=A(1:NA)
         RT(1:NA)=RA(1:NA)
         call MergeRI(T,NA,A(NA+1),NB,A,N,RT,RA(NA+1),RA)
      endif
      return

   end subroutine MergeSortRI

   ! Quick sort
   subroutine QSORT(A,IA,N)
      !
      implicit none
      integer, intent(in) :: N ! Number of entries that will be sorted
      real(dblprec), dimension(N), intent(inout) :: A ! Square modulus of the entries to be sorted (for vectors)
      integer, dimension(N), intent(out) :: IA  ! Indexing counters
      integer :: I,IL,IL1,IP,IQ,ISAV,ISTACK,IU,IU1,IX,IZ,J
      real(dblprec) :: X,XSAV,XX,Z,ZZ
      integer, dimension(20) :: ILT,IUT

      do I = 1,N
         IA(I) = I
      end do
      if (N == 1) return
      IL1 = 1
      IU1 = N
      ISTACK = 0
      OUTER_LOOP: do
         ISTACK = ISTACK + 1
         IL = IL1
         IU = IU1
         do
            IP = IL
            IQ = IU
            X = A(IP)
            Z = A(IQ)
            I = 0
            J = IQ - IP - 1
            if (X > Z) then
               XSAV = X
               X = Z
               Z = XSAV
               A(IP) = X
               A(IQ) = Z
               ISAV = IA(IP)
               IA(IP) = IA(IQ)
               IA(IQ) = ISAV
            end if
            if (IU-IL > 1) then
               XX = X
               IX = IP
               ZZ = Z
               IZ = IQ
               do
                  IP = IP + 1
                  if (IP >= IQ) exit
                  X = A(IP)
                  if (X >= XX) then
                     do
                        IQ = IQ - 1
                        if (IP >= IQ) then
                           IQ = IP
                           IP = IP - 1
                           Z = X
                           X = A(IP)
                           exit
                        else
                           Z = A(IQ)
                           if (Z <= ZZ) exit
                        end if
                     end do
                     if (X > Z) then
                        XSAV = X
                        X = Z
                        Z = XSAV
                        A(IP) = X
                        A(IQ) = Z
                        ISAV = IA(IP)
                        IA(IP) = IA(IQ)
                        IA(IQ) = ISAV
                     end if
                     if (X > XX) then
                        XX = X
                        I = I + 1
                        IX = IP
                     end if
                     if (Z < ZZ) then
                        ZZ = Z
                        I = I + 1
                        IZ = IQ
                     end if
                  end if
               end do
               IP = IQ - 1
               if (IP/=IX .and. abs(X-XX)<abs(XX)*dbl_tolerance) then
                  A(IP) = XX
                  A(IX) = X
                  ISAV = IA(IX)
                  IA(IX) = IA(IP)
                  IA(IP) = ISAV
               end if
               if (IQ/=IZ .and. abs(Z-ZZ)<abs(ZZ)*dbl_tolerance) then
                  A(IQ) = ZZ
                  A(IZ) = Z
                  ISAV = IA(IZ)
                  IA(IZ) = IA(IQ)
                  IA(IQ) = ISAV
               end if
               if (IU-IQ <= IP-IL) then
                  IU1 = IU
                  IL1 = IQ + 1
                  IU = IP - 1
               else
                  IL1 = IL
                  IU1 = IP - 1
                  IL = IQ + 1
               end if
               if (I /= J) then
                  if (IU1 > IL1) exit
                  if (IU > IL) cycle
               end if
            end if
            do
               ISTACK = ISTACK - 1
               if (ISTACK <= 0) exit OUTER_LOOP
               IL = ILT(ISTACK)
               IU = IUT(ISTACK)
               if (IU > IL) exit
            end do
         end do
         ILT(ISTACK) = IL
         IUT(ISTACK) = IU
      end do OUTER_LOOP
      !
   end subroutine QSORT


end module Sorting
