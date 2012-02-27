MODULE Almacenamiento
	CONTAINS

subroutine guardarNube(unidad, path, N, masas, coordenadas_x, coordenadas_y, coordenadas_z, v_x, v_y, v_z, densidades)
	use Vectores
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
		!No escribo si masa es menor que cero (se�al de que la part�cula fue destruida)
		
		if(masas(i) > 0 ) then
			write(unidad,*) coordenadas_x(i), ",", coordenadas_y(i), ",", coordenadas_z(i), ",", v_x(i), ",", v_y(i), ",", v_z(i), ",", radio, ",", masas(i), ",", densidades(i)
		endif
	enddo

	close(unit=unidad)

end subroutine


SUBROUTINE CargarNube(path, masas, coordenadas_x, coordenadas_y, coordenadas_z, v_x, v_y, v_z, distancias, densidades, N)
	
	use Constantes
	
	implicit none

	integer N, i
	real*8 distancias(0:N-1), masas(0:N-1), coordenadas_x(0:N-1), coordenadas_y(0:N-1), coordenadas_z(0:N-1), v_x(0:N-1), v_y(0:N-1), v_z(0:N-1), densidades(0:N-1)
	character(len=100) :: path

	open(unit=UNIT_NUBE, file=path, form='formatted', access='sequential')
	!  X,Y,Z,VX,VY,VZ,D,Masa,rho
	
	!Salto la linea del header
	
	read(UNIT_NUBE, *)
	
	do i = 0, N - 1, 1
		read(UNIT_NUBE, *) coordenadas_x(i), coordenadas_y(i), coordenadas_z(i), v_x(i), v_y(i), v_z(i), distancias(i), masas(i), densidades(i)
	enddo

	close(unit=UNIT_NUBE)

END SUBROUTINE

END MODULE Almacenamiento
