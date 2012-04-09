program principal
	use Constantes
	use Vectores
	use Funciones
	use Fisica
	use Tipos 
	use BarnesHut
	use Fuerzas	

	implicit none

	real*8 masa, radio, volumen, densidad, tff, acc, dt, temperatura, cs
	
	real*8 p1(0:2), p2(0:2), p3(0:2)
	
	real*8 p_vect(0:1)
	
	type(Particula) p, q
		
	volumen = 4.0D+0 * PI * (radio*PARSEC_CMS)**3
	masa =   9.97511514969359145E-6  * SOLAR_MASS_GR	
	densidad = 5.0E-18

	radio = calcularRadio(masa, densidad) / PARSEC_CMS	

	write(*,*) "Radio de particula: ", radio
end	
