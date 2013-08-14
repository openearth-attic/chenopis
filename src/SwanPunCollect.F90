subroutine SwanPunCollect ( blkndc )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
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
#ifdef HAVE_PUN
    !
    !   Authors
    !
    !   41.36: Marcel Zijlema
    !
    !   Updates
    !
    !   41.36, July 2012: New subroutine
    !
    !   Purpose
    !
    !   Determines for each grid vertex in global grid the node number
    !
    !   Method
    !
    !   Collect specific data from all nodes
    !
    !   Modules used
    !
        use mpi
        use ocpcomm4
        use m_parall, only: IAMMASTER
        use SwanGriddata, only: ivertg, nverts, nvertsg
        use SIZES, only: MYPROC, MNPROC
    !
        implicit none
    !
    !   Argument variables
    !
        real,    dimension(nvertsg), intent(out) :: blkndc   ! array giving node number in each grid vertex in global grid
    !
    !   Local variables
    !
        integer, save                            :: ient = 0 ! number of entries in this subroutine
        integer                                  :: ierr     ! error value of MPI call
        integer                                  :: j        ! loop counter
        integer                                  :: k        ! loop counter
        integer                                  :: nownv    ! number of vertices in own subdomain (without ghost vertices)
        !
        integer, dimension(:), allocatable       :: iarr     ! auxiliary integer array to gather data
        integer, dimension(:), allocatable       :: icount   ! array specifying array size of data received from each processor
        integer, dimension(:), allocatable       :: idsplc   ! array specifying the starting address of the incoming data from each processor, relative to the global array
        integer, dimension(:), allocatable       :: ivertp   ! vertex index of global grid in own subdomain (without ghost vertices)
        !
        real   , dimension(:), allocatable       :: arr      ! auxiliary real array to gather data
        real   , dimension(:), allocatable       :: blknd    ! node number per subdomain (without ghost vertices)
        !
        character(80)                            :: msgstr   ! string to pass message
    !
    !   Structure
    !
    !   Description of the pseudo code
    !
    !   Source text
    !
        if (ltrace) call strace (ient,'SwanPunCollect')
        !
        nownv = count(ivertg>0)
        !
        allocate(ivertp(nownv))
        allocate( blknd(nownv))
        !
        ! determine node number per subdomain
        !
        k = 0
        do j = 1, nverts
           if ( ivertg(j) > 0 ) then
              k = k + 1
              ivertp(k) = ivertg(j)
              blknd (k) = real(MYPROC+1)
           endif
        enddo
        !
        if ( IAMMASTER ) then
           allocate(icount(0:MNPROC-1))
           allocate(idsplc(0:MNPROC-1))
        endif
        !
        ! gather the array sizes to the master
        !
        call MPI_GATHER( nownv, 1, MPI_INTEGER, icount, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
        if ( ierr /= MPI_SUCCESS ) then
           write (msgstr, '(a,i6)') ' MPI produces some internal error - return code is ',ierr
           call msgerr( 4, trim(msgstr) )
           return
        endif
        !
        ! check consistency with respect to size of gathered data
        !
        if ( IAMMASTER ) then
           if ( sum(icount) /= nvertsg ) then
              call msgerr(4, 'inconsistency found in SwanPunCollect: size of gathered data not correct ')
              return
           endif
        endif
        !
        ! calculate starting address of each local array with respect to the global array
        !
        if ( IAMMASTER ) then
           idsplc(0) = 0
           do j = 1, MNPROC-1
              idsplc(j) = icount(j-1) + idsplc(j-1)
           enddo
        endif
        !
        if ( IAMMASTER ) then
           allocate(iarr(nvertsg))
           allocate( arr(nvertsg))
        endif
        !
        ! gather different amounts of data from each processor to the master
        !
        call MPI_GATHERV( ivertp, nownv, MPI_INTEGER, iarr, icount, idsplc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
        if ( ierr == MPI_SUCCESS ) call MPI_GATHERV( blknd, nownv, MPI_REAL, arr, icount, idsplc, MPI_REAL, 0, MPI_COMM_WORLD, ierr )
        if ( ierr /= MPI_SUCCESS ) then
           write (msgstr, '(a,i6)') ' MPI produces some internal error - return code is ',ierr
           call msgerr( 4, trim(msgstr) )
           return
        endif
        !
        if ( IAMMASTER ) then
           do j = 1, nvertsg
              blkndc(iarr(j)) = arr(j)
           enddo
        endif
        !
        deallocate(ivertp,blknd)
        if ( IAMMASTER ) deallocate(icount,idsplc,iarr,arr)
        !
#endif
end subroutine SwanPunCollect
