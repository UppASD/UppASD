!> Drivers for JSON writing
!> @copyright
!> GNU Public License.
!> @author
!> Anders Bergman
module json

   use Parameters
   use Profiling
   implicit none
   public

contains
         subroutine json_key(fileno,key)
            implicit none
            integer, intent(in) :: fileno
            character(*), intent(in) :: key
            !
            write(fileno,'(a20,a)',advance='no') '"'//trim(key)//'"','  :  '
         end subroutine json_key

         subroutine json_string(fileno,val)
            implicit none
            integer, intent(in) :: fileno
            character(*), intent(in) :: val
            !
            write(fileno,'(1x,a,a,a)') ' "',val,'" ,'
         end subroutine json_string

         subroutine json_char(fileno,val,nel,contflag,indflag)
            implicit none
            integer, intent(in) :: fileno
            integer, intent(in) :: nel
            character, dimension(nel), intent(in) :: val
            logical, intent(in), optional :: contflag
            logical, intent(in), optional :: indflag
            !
            integer :: iel
            logical :: cont 
            logical :: ind 

            !
            cont = .false.
            if (present(contflag)) cont=contflag
            ind = .false.
            if (present(indflag)) ind=indflag
            !
            if(nel==1) then
               write(fileno,'(a,a,a)') ' "',val(1),'" ,'
            else
               if(ind) then
                  write(fileno,'(26x,a)',advance='no') ' [ '
               else
                  write(fileno,'(1x,a)',advance='no') ' [ '
               end if
               do iel=1,nel-1
                  write(fileno,'(a,a,a)',advance='no') '"',val(iel),'" , '
               end do
               if (cont) then 
                  write(fileno,'(a,a,a)') '"',val(iel),'" ] '
               else
                  write(fileno,'(a,a,a)') '"',val(iel),'" ] , '
               end if
            end if
         end subroutine json_char

         subroutine json_float(fileno,val,nel,contflag,indflag)
            implicit none
            integer, intent(in) :: fileno
            integer, intent(in) :: nel
            real(dblprec), dimension(nel), intent(in) :: val
            logical, intent(in), optional :: contflag
            logical, intent(in), optional :: indflag
            !
            integer :: iel
            logical :: cont
            logical :: ind
            !
            cont = .false.
            if(present(contflag)) cont=contflag
            ind = .false.
            if(present(indflag)) ind=indflag

            if(nel==1) then
               write(fileno,'(g14.6,a)') val(1),' ,'
            else
               if(ind) then
                  write(fileno,'(26x,a)',advance='no') ' [ '
               else 
                  write(fileno,'(1x,a)',advance='no') ' [ '
               end if
               do iel=1,nel-1
                  write(fileno,'(g14.6,a)',advance='no') val(iel),' , '
               end do
               if(cont) then
                  write(fileno,'(g14.6,a)') val(iel),' ]  '
               else
                  write(fileno,'(g14.6,a)') val(iel),' ] , '
               end if
            end if
         end subroutine json_float

         subroutine json_int(fileno,val,nel,contflag,indflag)
            implicit none
            integer, intent(in) :: fileno
            integer, intent(in) :: nel
            integer, dimension(nel), intent(in) :: val
            logical, intent(in), optional :: contflag
            logical, intent(in), optional :: indflag
            !
            integer :: iel
            logical :: cont
            logical :: ind

            cont = .false.
            if(present(contflag)) cont=contflag
            ind = .false.
            if(present(indflag)) ind=indflag
            !
            if(nel==1) then
               write(fileno,'(i8,a)') val(1),' ,'
            else
               if(ind) then
                  write(fileno,'(26x,a)',advance='no') ' [ '
               else
                  write(fileno,'(1x,a)',advance='no') ' [ '
               end if
               do iel=1,nel-1
                  write(fileno,'(i8,a)',advance='no') val(iel),' , '
               end do
               if(cont) then
                   write(fileno,'(i8,a)') val(iel),' ]  '
               else
                   write(fileno,'(i8,a)') val(iel),' ] , '
                end if
            end if
         end subroutine json_int

!!!          subroutine json_vector(fileno,val,nel)
!!!             implicit none
!!!             integer, intent(in) :: fileno
!!!             integer, intent(in) :: nel
!!!             real(dblprec), dimension(nel), intent(in) :: val
!!!             !
!!!             integer :: iel
!!!             !
!!!             write(fileno,'(a)',advance='no') ' [ '
!!!             do iel=1,nel-1
!!!                write(fileno,'(f12.6,a)',advance='no') val(iel),' , '
!!!             end do
!!!             write(fileno,'(f12.6,a)') val(iel),' ] '
!!!          end subroutine json_vector
end module json
