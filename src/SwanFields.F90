module SwanFields
  !
  !   --|-----------------------------------------------------------|--
  !     | Delft University of Technology                            |
  !     | Faculty of Civil Engineering and Geosciences              |
  !     | Environmental Fluid Mechanics Section                     |
  !     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
  !     |                                                           |
  !     | Programmer: Fedor Baart                                   |
  !   --|-----------------------------------------------------------|--
  !
  !
  !     SWAN (Simulating WAves Nearshore); a third generation wave model
  !     Copyright (C) 1993-2012  Delft University of Technology
  !
  !     This program is free software; you can redistribute it and/or
  !     modify it under the terms of the GNU General Public License as
  !     published by the Free Software Foundation; either version 2 of
  !     the License, or (at your option) any later version.
  !
  !     This program is distributed in the hope that it will be useful,
  !     but WITHOUT ANY WARRANTY; without even the implied warranty of
  !     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  !     GNU General Public License for more details.
  !
  !     A copy of the GNU General Public License is available at
  !     http://www.gnu.org/copyleft/gpl.html#SEC3
  !     or by writing to the Free Software Foundation, Inc.,
  !     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  !
  !
  !   Authors
  !
  !   50.00: Fedor Baart
  !
  !   Updates
  !

  !   Purpose
  !
  !   Module containing field calculations for online data exchange
  !
  !   Method
  !
  !
  !   Modules used
  !
  !   none


  use iso_c_binding

  ! allow access to variables
  use swcomm1
  use swcomm2
  use swcomm3
  use m_genarr
  !
  implicit none
  !
  !   Module parameters
  !
  !
  !
  !   Module variables
  !

  logical(c_bool), allocatable, target :: mask(:,:) ! Land mask (0=land, 1=sea)
  logical(c_bool), parameter :: land =.false., sea=.true.

  !   ---
  !
  !   Source text
  !

contains
  subroutine field_init()
    ! TODO? do something with mxcgl, mycgl, iode to compute local bounds

    allocate(mask(mxc,myc))

  end subroutine field_init

  subroutine field_final()

    deallocate(mask)

  end subroutine field_final

  subroutine field_update(compda)
    real, intent(in) :: compda(mcgrd,mcmvar)
    call make_mask(compda, mask)

  end subroutine field_update

  subroutine make_mask(compda, mask)
    !   Allocate and set land/sea mask.  This mask is dynamic, includes
    !   halos and is based on both excluded points and depth.  Local
    !   mask values are used (1=sea, 0=land).

    real, intent(in) :: compda(mcgrd,mcmvar)
    logical(c_bool), allocatable, target, intent(out) :: mask(:,:)
    ! local
    integer :: ix, iy, i


    do iy = 1, myc
       do ix = 1, mxc
          i = kgrpnt(ix,iy)
          if (i.gt.1) then
             if (compda(i,jdp2).gt.depmin) then
                mask(ix,iy) = sea
             else
                mask(ix,iy) = land
             endif
          else
             mask(ix,iy) = land
          endif
       enddo
    enddo

  end subroutine make_mask


  end module SwanFields
