program principal
	
	use Almacenamiento

	implicit none
	
	real*8 masas(0:999), pos_x(0:999), pos_y(0:999), pos_z(0:999), distancias(0:999), v_x(0:999), v_y(0:999), v_z(0:999), densidades(0:999)
	character(len=256) :: path

	path = "./datos/sim1k0.1MSbeta0.16.nc"		
	
	write(*,*) "Path: ", path

	call CargarNube(path, masas, pos_x, pos_y, pos_z, v_x, v_y, v_z, distancias, densidades, 1000)

	stop
end
