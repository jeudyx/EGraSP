!--------------------------------------------
!EGraSP: Evolucion GRAvitacional de Sistemas de Particulas
!Maestria en Astrofisica, Universidad de Costa Rica - Jeudy Blanco
!Ultimo cambio: 27/02/2012 -----
!
!Modulo Almacenamiento: se encarga de salvar y recuperar el sistema hacia y desde el disco.
!--------------------------------------------

MODULE Almacenamiento
	CONTAINS


subroutine crearNubeNetCDF(path, N, masa, densidad, variacion, beta, tipo, altura)

	use Egrasp_NCIO

	implicit none
	
	integer N, i, ncid, status, tipo
	real*8	masa, densidad, variacion, beta, altura
	character(len=256) :: path
	type(runparameters) params

	params%title = "Star formation - collapse of Interstellar cloud"
	params%institution = "Centro de Investigaciones Espaciales- Universidad de Costa Rica"
	params%source = "Model simulation output"
	params%comment = ""
	params%totalmass = masa
	params%dt = 0.0D+0
	params%initial_density = densidad
	params%beta = beta
	params%temperature = -1.0D+0
	params%BH_theta = -1.0D+0
	params%N_neighbour = -1
	params%save_every = 0

	call init_ncio(path,N,params)

end subroutine

subroutine CerrarNubeCDF()
	use Egrasp_NCIO
	call release_ncio
end subroutine CerrarNubeCDF

subroutine guardarNube(unidad, path, N, masas, coordenadas_x, coordenadas_y, coordenadas_z, v_x, v_y, v_z, densidades)

	use Egrasp_NCIO

	implicit none
	
	integer N, i, unidad
	real*8 masas(0:N-1), coordenadas_x(0:N-1), coordenadas_y(0:N-1), coordenadas_z(0:N-1), v_x(0:N-1), v_y(0:N-1), v_z(0:N-1), densidades(0:N-1)
	real*8 radio
	real*8 vector_posicion(0:2)
	character(len=256) :: path	

	call guardarNubeNetCDF(N, masas, coordenadas_x, coordenadas_y, coordenadas_z, v_x, v_y, v_z, densidades)
end subroutine

SUBROUTINE CargarNube(path, masas, coordenadas_x, coordenadas_y, coordenadas_z, v_x, v_y, v_z, distancias, densidades, N)

	implicit none

	integer N, i
	real*8 distancias(0:N-1), masas(0:N-1), coordenadas_x(0:N-1), coordenadas_y(0:N-1), coordenadas_z(0:N-1), v_x(0:N-1), v_y(0:N-1), v_z(0:N-1), densidades(0:N-1)
	character(len=100) :: path

	call CargarNubeCSV(path, masas, coordenadas_x, coordenadas_y, coordenadas_z, v_x, v_y, v_z, distancias, densidades, N)
end subroutine

subroutine guardarNubeNetCDF(N, masas, coordenadas_x, coordenadas_y, coordenadas_z, v_x, v_y, v_z, densidades)

	use Egrasp_NCIO
	use Vectores

	implicit none
	
	integer N, i
	real*8 masas(0:N-1), coordenadas_x(0:N-1), coordenadas_y(0:N-1), coordenadas_z(0:N-1), v_x(0:N-1), v_y(0:N-1), v_z(0:N-1), densidades(0:N-1), distancias(0:N-1)
	real*8 radio
	real*8 vector_posicion(0:2)
	character(len=256) :: path

	do i = 0, N - 1, 1
		vector_posicion(0) = coordenadas_x(i)
		vector_posicion(1) = coordenadas_y(i)
		vector_posicion(2) = coordenadas_z(i)
		distancias(i) = magnitudVector3D(vector_posicion)
	enddo

	call writerec(N, coordenadas_x, coordenadas_y, coordenadas_z, v_x, v_y, v_z,masas,densidades,distancias)

end subroutine

subroutine guardarNubeCSV(unidad, path, N, masas, coordenadas_x, coordenadas_y, coordenadas_z, v_x, v_y, v_z, densidades)
	use Vectores
	
	implicit none
	
	integer N, i, unidad
	real*8 masas(0:N-1), coordenadas_x(0:N-1), coordenadas_y(0:N-1), coordenadas_z(0:N-1), v_x(0:N-1), v_y(0:N-1), v_z(0:N-1), densidades(0:N-1)
	real*8 radio
	real*8 vector_posicion(0:2)
	character(len=100) :: path
	
	open(unit=unidad, file=path)

	
	write(unidad,*) 'X,Y,Z,VX,VY,VZ,D,Masa,RHO'
	do i = 0, N - 1, 1
		vector_posicion(0) = coordenadas_x(i)
		vector_posicion(1) = coordenadas_y(i)
		vector_posicion(2) = coordenadas_z(i)
		radio = magnitudVector3D(vector_posicion)
		!No escribo si masa es menor que cero (señal de que la partícula fue destruida)
		
		if(masas(i) > 0 ) then
			write(unidad,*) coordenadas_x(i), ",", coordenadas_y(i), ",", coordenadas_z(i), ",", v_x(i), ",", v_y(i), ",", v_z(i), ",", radio, ",", masas(i), ",", densidades(i)
		endif
	enddo

	close(unit=unidad)

end subroutine


SUBROUTINE CargarNubeCSV(path, masas, coordenadas_x, coordenadas_y, coordenadas_z, v_x, v_y, v_z, distancias, densidades, N)
	
	use Constantes
	
	implicit none

	integer N, i
	real*8 distancias(0:N-1), masas(0:N-1), coordenadas_x(0:N-1), coordenadas_y(0:N-1), coordenadas_z(0:N-1), v_x(0:N-1), v_y(0:N-1), v_z(0:N-1), densidades(0:N-1)
	character(len=100) :: path

	open(unit=UNIT_NUBE, file=path, form='formatted', access='sequential')
	!open(unit=UNIT_NUBE, file=path, status='old', action='read')
	!  X,Y,Z,VX,VY,VZ,D,Masa,rho
	
	!Salto la linea del header
	
	read(UNIT_NUBE, *)
	
	do i = 0, N - 1, 1
		read(UNIT_NUBE, *) coordenadas_x(i), coordenadas_y(i), coordenadas_z(i), v_x(i), v_y(i), v_z(i), distancias(i), masas(i), densidades(i)
!		format (F25.15,F25.15,F25.15,F25.15,F25.15,F25.15,F25.15,F25.15,ES20.5E2)
	enddo

	close(unit=UNIT_NUBE)

END SUBROUTINE

END MODULE Almacenamiento
