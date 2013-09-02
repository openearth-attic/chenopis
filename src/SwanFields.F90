module SwanFields
  !
  !   50.00: Fedor Baart
  !   Purpose
  !   Module containing field calculations for online data exchange
  !
  !   Method
  !   Modules used
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
  !   Module variables
  !

  logical(c_bool), allocatable, target, save :: landmask(:,:) ! Land landmask (0=land, 1=sea)
  real(c_float), allocatable, target, save :: tpsmooth(:)
  real(c_float), allocatable, target, save :: ubot(:)
  real(c_float), allocatable, target, save :: hsig(:)
  real(c_float), allocatable, target, save :: dep(:)
  real(c_float), allocatable, target, save :: k(:,:)
  real(c_float), allocatable, target, save :: cg(:,:)
  real(c_float), allocatable, target, save :: n(:,:)
  real(c_float), allocatable, target, save :: nd(:,:)
  logical(c_bool), parameter :: land =.false., sea=.true.

  ! not in a module....
  real, external ::  SwanIntgratSpc

  !   ---
  !
  !   source text
  !

contains
  subroutine field_init()
    ! TODO? do something with mxcgl, mycgl, iode to compute local bounds

    allocate(landmask(mxc,myc)) ! landmask
    allocate(tpsmooth(mcgrd)) ! peak period
    allocate(ubot(mcgrd)) ! orbital velocity
    allocate(hsig(mcgrd)) ! significant wave height
    allocate(dep(mcgrd)) ! depth
    allocate(k(msc,mcgrd)) ! wave number
    allocate(cg(msc,mcgrd)) ! group velocity
    allocate(n(msc,mcgrd)) ! ratio of group and phase velocity
    allocate(nd(msc,mcgrd)) ! derivative of N with respect to D

  end subroutine field_init

  subroutine field_finalize()

    deallocate(landmask)

  end subroutine field_finalize

  subroutine field_update(compda)
    real, intent(in) :: compda(mcgrd,mcmvar)
    ! Update module variables
    call make_dep(compda) ! needed for make_k
    call make_k_and_cg_and_n_and_nd(compda)
    call make_landmask(compda)
    call make_tpsmooth(compda)
    call make_ubot(compda)
    call make_hsig(compda)

  end subroutine field_update

  subroutine make_dep(compda)
    implicit none
    real, intent(in) :: compda(mcgrd,mcmvar)

    dep = compda(1,jdp2)
  end subroutine make_dep

  subroutine make_k_and_cg_and_n_and_nd(compda)

    implicit none
    real, intent(in) :: compda(mcgrd,mcmvar)

    integer :: ip
    ! variables
    ! mmt, sig, d, k, cg, n, nd
    ! int:: mmt
    ! real:: sig(mmt), d, k(mmt), cg(mmt), n(mmt), nd(mmt)
    ! intent::  in, in, in, out, out, out, out, out

    do ip=1,mcgrd
       call kscip1 (msc, spcsig, dep(ip), k(:,ip), cg(:,ip), n(:,ip), nd(:,ip))
    end do

  end subroutine make_k_and_cg_and_n_and_nd

  subroutine make_hsig(compda)
    implicit none
    real, intent(in) :: compda(mcgrd,mcmvar)
    real :: ds ! delta spectrum
    real :: ead ! energy per bin, per point, per direction
    real :: ehfr, etot, eftail, ecs
    real :: fmin, fmax
    integer :: ip, id, is


    eftail = 1. / (pwtail(1) - 1.)

    do ip=1,mcgrd
       if (outpar(6).eq.0.) then
          ! integration over [0,inf]
          etot = 0.
          ! trapezoidal rule is applied
          do id=1, mdc
             do is=2,msc
                ds = spcsig(is)-spcsig(is-1)
                ead = 0.5*(spcsig(is)*ac2(id,is,ip)+ spcsig(is-1)*ac2(id,is-1,ip))*ds*ddir
                etot = etot + ead
             enddo
             if (msc .gt. 3) then
                ! contribution of tail to total energy density
                ehfr = ac2(id,msc,ip) * spcsig(msc)
                etot = etot + ddir * ehfr * spcsig(msc) * eftail
             endif
          enddo
       else
          ! integration over [fmin,fmax]
          fmin = pi2*outpar(21)
          fmax = pi2*outpar(36)
          ecs  = 1.
          etot = SwanIntgratSpc(0., fmin, fmax, spcsig, spcdir(1,1), k(:,ip), ecs, 0., 0., ac2(:,:,ip), 1)
       endif

       if (etot .ge. 0.) then
          hsig(ip) = 4.*sqrt(etot)
       else
          hsig(ip) = 0.
       endif
    end do

  end subroutine make_hsig

  subroutine make_ubot(compda)

    implicit none
    real, intent(in) :: compda(mcgrd,mcmvar)
    integer :: ip ! counter over grid points
    real :: rr
    rr = sqrt(2.)
    do ip=1,mcgrd
       ubot(:) = compda(1, jubot)*rr
    enddo
  end subroutine make_ubot


  subroutine make_tpsmooth(compda)

    implicit none
    real, intent(in) :: compda(mcgrd,mcmvar)

    real :: emax  ! Maximum energy according to the depth and the breaker parameter ?
    real :: etd ! sum of energy over theta direction?
    real :: ed, emaxd ! TODO what is this?
    real :: emaxu ! TODO what is this?
    real :: e1, e2, e3 ! TODO what is this?
    integer :: isigm ! location of maximum in sigma direction?
    real :: sig1, sig2, sig3 ! neighbours for interpolation?
    real :: sigp ! the result?
    integer :: is ! counter over grid points in sigma-direction of computational grid
    integer :: id ! counter over grid points in theta-direction of computational grid
    integer :: ip ! counter over grid points
    real :: p, q, r, t, a ! some more unspecified terms...

    do ip=1,mcgrd
       emax = 0.
       etd  = 0.
       isigm = -1
       do is = 1, msc
          ed  = etd
          etd = 0.
          do id = 1, mdc
             etd = etd + spcsig(is)*ac2(id,is,ip)*ddir
          end do
          if (etd.gt.emax) then
             emax  = etd
             isigm = is
             emaxd = ed
             emaxu = 0.
             if (is.lt.msc) then
                do id = 1, mdc
                   emaxu = emaxu + spcsig(is+1)*ac2(id,is+1,ip)*ddir
                end do
             else
                emaxu = emax
             end if
          end if
       end do


       if (isigm.gt.1 .and. isigm.lt.msc) then
          sig1 = spcsig(isigm-1)
          sig2 = spcsig(isigm+1)
          sig3 = spcsig(isigm  )
          e1   = emaxd
          e2   = emaxu
          e3   = emax
          p    = sig1+sig2
          q    = (e1-e2)/(sig1-sig2)
          r    = sig1+sig3
          t    = (e1-e3)/(sig1-sig3)
          a    = (t-q)/(r-p)
          if (a.lt.0) then
             sigp = (-q+p*a)/(2.*a)
          else
             sigp = sig3
          end if
          tpsmooth(ip) = 2.*pi/sigp
       else if (isigm.eq.1) then
          tpsmooth(ip) = 2.*pi/spcsig(1)
       else
          tpsmooth(ip) = -9
       end if
    end do

  end subroutine make_tpsmooth

  subroutine make_landmask(compda)
    !   Allocate and set land/sea landmask.  This landmask is dynamic, includes
    !   halos and is based on both excluded points and depth.  Local
    !   landmask values are used (1=sea, 0=land).

    real, intent(in) :: compda(mcgrd,mcmvar)
    ! local
    integer :: ix, iy, i


    do iy = 1, myc
       do ix = 1, mxc
          i = kgrpnt(ix,iy)
          if (i.gt.1) then
             if (compda(i,jdp2).gt.depmin) then
                landmask(ix,iy) = sea
             else
                landmask(ix,iy) = land
             endif
          else
             landmask(ix,iy) = land
          endif
       enddo
    enddo

  end subroutine make_landmask


end module SwanFields
