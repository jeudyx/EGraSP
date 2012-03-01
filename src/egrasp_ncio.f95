module Egrasp_NCIO

	use netcdf
	use Tipos

	private

	integer , parameter :: filelen = 256

	character(len=32) :: reference = 'http://github.com/jeudyx/EGraSP'

	double precision , parameter :: fillvalue = 1D+20

	! NetCDF ID
	integer :: ncid
	! Error Status
	integer :: stat

	! Variable IDs
	integer :: time_id
	integer :: position_X_id
	integer :: position_Y_id
	integer :: position_Z_id
	integer :: velocity_U_id
	integer :: velocity_V_id
	integer :: velocity_W_id
	integer :: mass_id
	integer :: density_id
	integer :: distance_id

	integer :: idrec = 1

	double precision :: xtimefac

	public :: runparameters , update_ncio, init_ncio , writerec , release_ncio , open_ncio

contains

	subroutine check_err
		implicit none
		if (stat /= NF90_NOERR) then
			print *, nf90_strerror(stat)
			stop
		endif
	end subroutine check_err

	subroutine update_ncio(filename,rp)
		implicit none
		character (len=filelen) , intent(in) :: filename
		type(runparameters) , intent(in) :: rp 
		stat = nf90_open(filename, NF90_WRITE, ncid)
		call check_err

		stat = nf90_put_att(ncid, NF90_GLOBAL, 'model_dt', rp%dt)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'model_temperature', rp%temperature)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'model_BH_theta', rp%BH_theta)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'model_N_neighbour', rp%N_neighbour)
		call check_err

		stat = nf90_sync(ncid)
		call check_err		

	end subroutine update_ncio

	subroutine init_ncio(filename,NP_len,rp)
		implicit none
		character (len=filelen) , intent(in) :: filename
		integer , intent(in) :: NP_len
		type(runparameters) , intent(in) :: rp 
		integer :: NP_dim
		integer :: time_dim
		character (len=64) :: history
		integer , dimension(8) :: xtime

		integer , dimension(2) :: idims

		xtimefac = dble(rp%save_every)
		stat = nf90_create(filename, nf90_clobber, ncid);
		call check_err

		stat = nf90_def_dim(ncid, 'NP', NP_len, NP_dim);
		call check_err
		stat = nf90_def_dim(ncid, 'time', nf90_unlimited, time_dim);
		call check_err

		idims(1) = NP_dim
		idims(2) = time_dim

		stat = nf90_def_var(ncid, 'time', nf90_double, idims(2:2), time_id)
		call check_err

		stat = nf90_def_var(ncid, 'position_X', nf90_double, idims, position_X_id)
		call check_err
		stat = nf90_def_var(ncid, 'position_Y', nf90_double, idims, position_Y_id)
		call check_err
		stat = nf90_def_var(ncid, 'position_Z', nf90_double, idims, position_Z_id)
		call check_err

		stat = nf90_def_var(ncid, 'velocity_U', nf90_double, idims, velocity_U_id)
		call check_err
		stat = nf90_def_var(ncid, 'velocity_V', nf90_double, idims, velocity_V_id)
		call check_err
		stat = nf90_def_var(ncid, 'velocity_W', nf90_double, idims, velocity_W_id)
		call check_err

		stat = nf90_def_var(ncid, 'mass', nf90_double, idims, mass_id)
		call check_err
		stat = nf90_def_var(ncid, 'density', nf90_double, idims, density_id)
		call check_err
		stat = nf90_def_var(ncid, 'distance', nf90_double, idims, distance_id)
		call check_err

		stat = nf90_put_att(ncid, NF90_GLOBAL, 'title', rp%title)
		call check_err

		call date_and_time(values=xtime)
		write (history,'(a,i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')	 &
						 "Created on ", xtime(1) , '-' , xtime(2) , '-' , xtime(3) , &
						 ' ' , xtime(5) , ':' , xtime(6) , ':' , xtime(7)
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'history', history)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'institution', rp%institution)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'source', rp%source)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'references', reference)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'comment', rp%comment)
		call check_err

		stat = nf90_put_att(ncid, NF90_GLOBAL, 'model_dt', rp%dt)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'model_totalmass', rp%totalmass)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'model_initial_density', rp%initial_density)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'model_beta', rp%beta)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'model_temperature', rp%temperature)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'model_BH_theta', rp%BH_theta)
		call check_err
		stat = nf90_put_att(ncid, NF90_GLOBAL, 'model_N_neighbour', rp%N_neighbour)
		call check_err

		stat = nf90_put_att(ncid, time_id, 'units', 'years since 0000-00-00')
		call check_err

		stat = nf90_put_att(ncid, position_X_id, 'units', 'parsec')
		call check_err
		stat = nf90_put_att(ncid, position_X_id, '_FillValue', fillvalue)
		call check_err
		stat = nf90_put_att(ncid, position_Y_id, 'units', 'parsec')
		call check_err
		stat = nf90_put_att(ncid, position_Y_id, '_FillValue', fillvalue)
		call check_err
		stat = nf90_put_att(ncid, position_Z_id, 'units', 'parsec')
		call check_err
		stat = nf90_put_att(ncid, position_Z_id, '_FillValue', fillvalue)
		call check_err

		stat = nf90_put_att(ncid, velocity_U_id, 'units', 'm/s')
		call check_err
		stat = nf90_put_att(ncid, velocity_U_id, '_FillValue', fillvalue)
		call check_err
		stat = nf90_put_att(ncid, velocity_V_id, 'units', 'm/s')
		call check_err
		stat = nf90_put_att(ncid, velocity_V_id, '_FillValue', fillvalue)
		call check_err
		stat = nf90_put_att(ncid, velocity_W_id, 'units', 'm/s')
		call check_err
		stat = nf90_put_att(ncid, velocity_W_id, '_FillValue', fillvalue)
		call check_err

		stat = nf90_put_att(ncid, mass_id, 'units', 'solar mass')
		call check_err
		stat = nf90_put_att(ncid, mass_id, '_FillValue', fillvalue)
		call check_err
		stat = nf90_put_att(ncid, density_id, 'units', 'g/cm^3')
		call check_err
		stat = nf90_put_att(ncid, density_id, '_FillValue', fillvalue)
		call check_err
		stat = nf90_put_att(ncid, distance_id, 'units', 'parsec')
		call check_err
		stat = nf90_put_att(ncid, distance_id, '_FillValue', fillvalue)
		call check_err

		stat = nf90_enddef(ncid);
		call check_err

	end subroutine init_ncio

	subroutine release_ncio
		stat = nf90_close(ncid)
		call check_err
		ncid = -1
	end subroutine release_ncio

	subroutine writerec(NP,X,Y,Z,U,V,W,mass,rho,dist)
		implicit none
		integer , intent(in) :: NP
		real*8 X(1:NP)		
		real*8 Y(1:NP)
		real*8 Z(1:NP)
		real*8 U(1:NP)
		real*8 V(1:NP)
		real*8 W(1:NP)
		real*8 mass(1:NP)
		real*8 rho(1:NP)
		real*8 dist(1:NP)
		integer , dimension(2) :: istart , icount
		double precision , dimension(1) :: xtime

		xtime(1) = dble(idrec) * xtimefac
		istart(:) = idrec
		icount(:) = 1
	
		stat = nf90_put_var(ncid, time_id, xtime, istart(1:1), icount(1:1))
		call check_err

		istart(2) = 1
		icount(2) = 1		!Recuerde, la dimension temporal es la segunda, la primera es de los puntos!

		icount(1) = NP

		stat = nf90_put_var(ncid, position_X_id, X, istart, icount)
		call check_err
		stat = nf90_put_var(ncid, position_Y_id, Y, istart, icount)
		call check_err
		stat = nf90_put_var(ncid, position_Z_id, Z, istart, icount)
		call check_err
		stat = nf90_put_var(ncid, velocity_U_id, U, istart, icount)
		call check_err
		stat = nf90_put_var(ncid, velocity_V_id, U, istart, icount)
		call check_err
		stat = nf90_put_var(ncid, velocity_W_id, W, istart, icount)
		call check_err
		stat = nf90_put_var(ncid, mass_id, mass, istart, icount)
		call check_err
		stat = nf90_put_var(ncid, density_id, rho, istart, icount)
		call check_err
		stat = nf90_put_var(ncid, distance_id, dist, istart, icount)
		call check_err
		stat = nf90_sync(ncid)
		call check_err
		idrec = idrec + 1
	end subroutine writerec

	subroutine open_ncio(NP,filename,X,Y,Z,U,V,W,mass,rho,dist)
		implicit none
		integer , intent(in) :: NP
		character (len=filelen) , intent(in) :: filename
		real*8 X(1:NP)		
		real*8 Y(1:NP)
		real*8 Z(1:NP)
		real*8 U(1:NP)
		real*8 V(1:NP)
		real*8 W(1:NP)
		real*8 mass(1:NP)
		real*8 rho(1:NP)
		real*8 dist(1:NP)
		integer , dimension(2) :: istart , icount
		double precision , dimension(2) :: xtime
		integer :: NPlen , NP_dimid , time_dimid

		stat = nf90_open(filename, NF90_WRITE, ncid)
		call check_err

		stat = nf90_inq_dimid(ncid, "NP", NP_dimid)
		call check_err
		stat = nf90_inquire_dimension(ncid, NP_dimid, len=NPlen)
		call check_err
		if ( NP /= NPlen ) then
			print *, 'Error ! Dimension of stored and requested differ !'
			stop
		end if
		stat = nf90_inq_dimid(ncid, "time", time_dimid)
		call check_err
		stat = nf90_inquire_dimension(ncid, NP_dimid, len=idrec)
		call check_err

		stat = nf90_inq_varid(ncid, 'time', time_id)
		call check_err
		istart(1) = 1
		icount(1) = 2
		stat = nf90_get_var(ncid, time_id, xtime, istart(1:1), icount(1:1))
		call check_err
		xtimefac = xtime(2) - xtime(1)
		stat = nf90_inq_varid(ncid, 'position_X', position_X_id)
		call check_err
		stat = nf90_inq_varid(ncid, 'position_Y', position_Y_id)
		call check_err
		stat = nf90_inq_varid(ncid, 'position_Z', position_Z_id)
		call check_err
		stat = nf90_inq_varid(ncid, 'velocity_U', velocity_U_id)
		call check_err
		stat = nf90_inq_varid(ncid, 'velocity_V', velocity_V_id)
		call check_err
		stat = nf90_inq_varid(ncid, 'velocity_W', velocity_W_id)
		call check_err
		stat = nf90_inq_varid(ncid, 'mass', mass_id)
		call check_err
		stat = nf90_inq_varid(ncid, 'density', density_id)
		call check_err
		stat = nf90_inq_varid(ncid, 'distance', distance_id)
		call check_err
		istart(1) = idrec
		istart(2) = 1
		icount(1) = 1
		icount(2) = NP
		stat = nf90_get_var(ncid, position_X_id, X, istart, icount)
		call check_err
		stat = nf90_get_var(ncid, position_Y_id, Y, istart, icount)
		call check_err
		stat = nf90_get_var(ncid, position_Z_id, Z, istart, icount)
		call check_err
		stat = nf90_get_var(ncid, velocity_U_id, U, istart, icount)
		call check_err
		stat = nf90_get_var(ncid, velocity_V_id, V, istart, icount)
		call check_err
		stat = nf90_get_var(ncid, velocity_W_id, W, istart, icount)
		call check_err
		stat = nf90_get_var(ncid, mass_id, mass, istart, icount)
		call check_err
		stat = nf90_get_var(ncid, density_id, rho, istart, icount)
		call check_err
		stat = nf90_get_var(ncid, distance_id, dist, istart, icount)
		call check_err
		idrec = idrec + 1
	end subroutine open_ncio

end module egrasp_ncio
