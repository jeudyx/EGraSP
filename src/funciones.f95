MODULE Funciones
	CONTAINS


function asignarPosicion(Zcurr, rndX, rndY, rndZ) result(posicion)

	implicit none

	real*8 Zcurr, rndX, rndY, rndZ, rndS
	
	real*8 posicion(0:2)
	
	integer signo	

	CALL RANDOM_NUMBER(rndS)

	if(rndS <= 0.5) then
		signo = 1
	else
		signo = -1
	endif

	posicion(0) = frand(0.0D+0, Zcurr, rndX) * signo

	!------------
	
	CALL RANDOM_NUMBER(rndS)

	if(rndS <= 0.5) then
		signo = 1
	else
		signo = -1
	endif

	posicion(1) = frand(0.0D+0, Zcurr, rndY) * signo
	
	!------------
	
	if(rndZ >= 0.0D+0) then
	
		CALL RANDOM_NUMBER(rndS)

		if(rndS <= 0.5) then
			signo = 1
		else
			signo = -1
		endif

		posicion(2) = frand(0.0D+0, Zcurr, rndZ) * signo
	else
		posicion(2) = 0.0D+0
	endif

	return

end function

function asignarPosicion2(iniciox, inicioy, inicioz, finalx, finaly, finalz, rndX, rndY, rndZ) result(posicion)

	implicit none

	real*8 iniciox, inicioy, inicioz, finalx, finaly, finalz, rndX, rndY, rndZ, rndS
	
	real*8 posicion(0:2)
	
	integer signo	

	CALL RANDOM_NUMBER(rndS)

	if(rndS <= 0.5) then
		signo = 1
	else
		signo = -1
	endif

	posicion(0) = frand(iniciox, finalx, rndX) * signo

	!------------
	
	CALL RANDOM_NUMBER(rndS)

	if(rndS <= 0.5) then
		signo = 1
	else
		signo = -1
	endif

	posicion(1) = frand(inicioy, finaly, rndY) * signo
	
	!------------
	
	if(rndZ >= 0.0D+0) then
	
		CALL RANDOM_NUMBER(rndS)

		if(rndS <= 0.5) then
			signo = 1
		else
			signo = -1
		endif

		posicion(2) = frand(inicioz, finalz, rndZ) * signo
	else
		posicion(2) = 0.0D+0
	endif

	return

end function

subroutine imprimirNube(N, masas, coordenadas_x, coordenadas_y, coordenadas_z, v_x, v_y, v_z)
	use Vectores
	integer N, i
	real*8 masas(0:N-1), coordenadas_x(0:N-1), coordenadas_y(0:N-1), coordenadas_z(0:N-1), v_x(0:N-1), v_y(0:N-1), v_z(0:N-1)
	real*8 radio
	real*8 vector_posicion(0:2)

	!type(Particula) nube(N)
	write(*,*) 'X,Y,Z,VX,VY,VZ,D,Masa'
	do i = 0, N - 1, 1
		vector_posicion(0) = coordenadas_x(i)
		vector_posicion(1) = coordenadas_y(i)
		vector_posicion(2) = coordenadas_z(i)
		radio = magnitudVector3D(vector_posicion)
		!write(*,*) i, ",", coordenadas_x(i), ",", coordenadas_y(i), ",", coordenadas_z(i), ",",radio, ",", masas(i)
		write(*,*) coordenadas_x(i), ",", coordenadas_y(i), ",", coordenadas_z(i), ",", v_x(i), ",", v_y(i), ",", v_z(i), ",", radio, ",", masas(i)
	enddo

end subroutine

!-----------------------------------------------
!
!Distancia entre 2 particulas en metros
!-----------------------------------------------

real*8 function distanciaParticulas(particula1, particula2)
	use Constantes
	use Tipos

	implicit none

	type(Particula) particula1, particula2

	!write(*,*) "En distanciaParticulas1", particula1%posicion(0), particula2%posicion(0), particula1%posicion(0) - particula2%posicion(0), (particula1%posicion(0) - particula2%posicion(0))**2
	!write(*,*) "En distanciaParticulas2", particula1%posicion(1), particula2%posicion(1), particula1%posicion(1) - particula2%posicion(1), ((particula1%posicion(1) - particula2%posicion(1))**2)
	!write(*,*) "En distanciaParticulas3", particula1%posicion(2), particula2%posicion(2), particula1%posicion(2) - particula2%posicion(2), ((particula1%posicion(2) - particula2%posicion(2))**2)

	distanciaParticulas = PARSEC_MTS * SQRT( ((particula1%posicion(0) - particula2%posicion(0))**2) + ((particula1%posicion(1) - particula2%posicion(1))**2) + ((particula1%posicion(2) - particula2%posicion(2))**2) )
end function

!-----------------------------------------------
!Fuerza de gravedad entre 2 particulas
!Radio en metros, masa en kilogramos
!-----------------------------------------------
real*8 function fuerzaGravedad(particula1, particula2)
	use Constantes
	use Tipos

	implicit none

	type(Particula) particula1, particula2

	fuerzaGravedad = (G * particula1%masa * particula2%masa) / distanciaParticulas(particula1, particula2)

end function
!-----------------------------------------------
function centroMasaN(masas, coordenadas_x, coordenadas_y, coordenadas_z, N) result(cm)
	use Constantes
	use Tipos
	implicit none

	integer i, N
	real*8 masas(0:N-1), coordenadas_x(0:N-1), coordenadas_y(0:N-1), coordenadas_z(0:N-1)
	
	real*8 masaTotal
	
	real*8 cm(0:2)
	
	masaTotal = 0.0D+0
	
	cm(0) = 0.0D+0
	cm(1) = 0.0D+0
	cm(2) = 0.0D+0
	
	do i = 0, N - 1, 1
		masaTotal = masaTotal + masas(i)
	enddo
	
	do i = 0, N - 1, 1
		cm(0) = cm(0) + ((coordenadas_x(i) * masas(i)) / masaTotal)
		cm(1) = cm(1) + ((coordenadas_y(i) * masas(i)) / masaTotal)
		cm(2) = cm(2) + ((coordenadas_z(i) * masas(i)) / masaTotal)
	enddo
		
	
	return
	
end function

!-----------------------------------------------

!Construye una particula
!Masa la recibe en masas solares, la almacena en kg
!Posiciones vienen en parsec
!Densidad en gr/cm3
!-----------------------------------------------

function construirParticula(id, masa, x, y, z, vx, vy, vz, densidad) result(particulaX)
	use Constantes
	use Tipos
	implicit none

	type(Particula) particulaX

	integer :: id
	real*8 :: masa, x, y, z, vx, vy, vz
	real*8 :: densidad

	particulaX%posicion(0) = x
	particulaX%posicion(1) = y
	particulaX%posicion(2) = z
	particulaX%velocidad(0) = vx
	particulaX%velocidad(1) = vy
	particulaX%velocidad(2) = vz
	particulaX%id = id
	particulaX%masa = masa * SOLAR_MASS_KG	
	particulaX%densidad = densidad

	return

end function

!-----------------------------------------------
!Funcion auxiliar para producir un numero aleatorio entre low y high a partir de otro numero entre 0 y 1 (rand)
!-----------------------------------------------

real*8 function frand(low, high, rand)

	implicit none

	real*8 :: low, high, rand

	frand = low + (rand * (high - low))
	
	return
end function


!-----------------------------------------------
!Funcion para asignar la masa a una particula, a partir de una Masa total y un N de particulas pendientes|
!la masa de cada particula será M/N +- un % variacion
!-----------------------------------------------

real*8 function asignarMasa(M, N, Var)

	implicit none	

	!real*8, external :: frand

	real*8 :: M, Var, mi, rnd, rnd2, rndrg

	integer :: N, signo

	if(N == 1) then
		asignarMasa = M
		return
	end if

	rnd = 0.0D+0
	call random_number(rnd)
	call random_number(rnd2)
	
	if(rnd2 <= 0.5) then
		signo = 1
	else
		signo = -1
	endif

	rndrg = frand(0.0D+0, Var, rnd)

	mi = M/N	
	
	!write (*,*) 'mi uniforme: ', mi, ' - Random: ', rnd, 'Random en rango: ', rndrg

	mi = mi + (mi * rndrg * signo)

	if(mi < 0.0) then
		mi = M/N
	endif

	M = M - mi
	N = N - 1	

	asignarMasa = mi
	
	!write (*,*) 'Masa asignada: ', mi, '- Masa restante: ', M

	return

end function

SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
END SUBROUTINE


function resolverCuadratica(a, b, c) result(x)

	implicit none
	
	real*8 a, b, c, discriminante
	real*8 x(0:1)
	
	discriminante = b**2 - (4.0D+0 * a * c)
	
	write(*,*) "En resolverCuadratica. a, b, c = ", a, b, c, " discriminante: ", discriminante
	
	x(0) = ( (-1.0D+0 * b) + (SQRT(discriminante)) ) / (2.0D+0 * a)
	x(1) = ( (-1.0D+0 * b) - (SQRT(discriminante)) ) / (2.0D+0 * a)

	return

end function

END MODULE Funciones
