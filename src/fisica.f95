MODULE Fisica
	CONTAINS
 
!Posiciones r1 y r2 en metros, velocidades v1 y v2 en metros por segundo
real*8 function Viscosidad_miu(r1, r2, v1, v2, soft_len)
	
	use Constantes
	use Vectores
	
	implicit none
	
	real*8 r1(0:2), r2(0:2), v1(0:2), v2(0:2), vij(0:2), rij(0:2)
	real*8 soft_len, magrij
	
	rij = diferenciaVectores3D(r1, r2)
	vij = diferenciaVectores3D(v1, v2)
	
	!Uso un soft len similar al gravitacional
	
	magrij = magnitudVector3D(rij)
	
	if(magrij >= soft_len) then
		Viscosidad_miu = productoEscalar3D(vij, rij) / magrij
	else
		Viscosidad_miu = productoEscalar3D(vij, rij) / soft_len
	endif
		
	!05/01/2012 Lo cambie quitando el cuadrado del denominador porque en las ecuaciones, h (de SPH) viene en numerador
	!y aqui no hay h, que viene siendo proporcional a distancia entre particulas vecinas
		

end function

function ArtificialViscosityAcc(idx_particula, temperatura, soft_len, n_vecinos, lista_vecinos, N, pos_x, pos_y, pos_z, v_x, v_y, v_z, densidades, densidades_locales, masas) result(aceleracion)

	use Constantes
	use Vectores
	
	implicit none
	
	integer idx_particula, n_vecinos, N, j
	integer lista_vecinos(0:n_vecinos-1)		
	real*8 pos_x(0:N-1), pos_y(0:N-1), pos_z(0:N-1), v_x(0:N-1), v_y(0:N-1), v_z(0:N-1), densidades_locales(0:N-1), masas(0:N-1), densidades(0:N-1)
	real*8 r1(0:2), r2(0:2), v1(0:2), v2(0:2), aceleracion(0:2), aceleracion_tmp(0:2), r12(0:2)
	real*8 resultado, total, densidad_local1, densidad_local2, densidad_local_promedio, temperatura, soft_len, radio1, radio2, mag_r12, factor

	resultado = 0.0D+0
	
	r1(0) = pos_x(idx_particula) * PARSEC_MTS
	r1(1) = pos_y(idx_particula) * PARSEC_MTS
	r1(2) = pos_z(idx_particula) * PARSEC_MTS
	
	v1(0) = v_x(idx_particula)
	v1(1) = v_y(idx_particula)
	v1(2) = v_z(idx_particula)
	
	aceleracion = 0.0D+0
	total = 0.0D+0
	densidad_local1 = densidades_locales(idx_particula)
				
	do j = 0, n_vecinos - 1, 1
						
		r2(0) = pos_x(lista_vecinos(j)) * PARSEC_MTS
		r2(1) = pos_y(lista_vecinos(j)) * PARSEC_MTS
		r2(2) = pos_z(lista_vecinos(j)) * PARSEC_MTS

		v2(0) = v_x(lista_vecinos(j))
		v2(1) = v_y(lista_vecinos(j))
		v2(2) = v_z(lista_vecinos(j))
		
		mag_r12 = magnitudVector3D(diferenciaVectores3D(r1, r2))

		!Vector unitario
		r12 = diferenciaVectores3D(r1, r2) / mag_r12
				
		densidad_local2 = densidades_locales(lista_vecinos(j))

		!write(*,*) "En ArtificialViscosityAcc: ", idx_particula, " densidades locales: ", densidad_local1, densidad_local2
		
		densidad_local_promedio = (densidad_local1 + densidad_local2) / 2.0D+0
		
		!Convierte a kg/m^3
		densidad_local_promedio = densidad_local_promedio / 1000.0D+0
		
		resultado = ArtificialViscosityPoints(r1, r2, v1, v2, temperatura, densidad_local1, densidad_local2, soft_len) 
			
		!write(*,*) idx_particula, lista_vecinos(j), " Viscosidad ij -> II: ", resultado, " mag_r12 ", mag_r12, (mag_r12**4)
		

		if(mag_r12 >= soft_len) then
			resultado = resultado * densidad_local_promedio / mag_r12
		else
			resultado = resultado * densidad_local_promedio / soft_len
		endif
		
		aceleracion_tmp(0) = resultado * r12(0)
		aceleracion_tmp(1) = resultado * r12(1)
		aceleracion_tmp(2) = resultado * r12(2)
		
		aceleracion(0) = aceleracion(0) + aceleracion_tmp(0)
		aceleracion(1) = aceleracion(1) + aceleracion_tmp(1)
		aceleracion(2) = aceleracion(2) + aceleracion_tmp(2)
	enddo
	
!	if( (densidad_local1 / densidades(idx_particula)) > 1.0D+0) then
!		write(*,*) "UY, llego al factor!: ", densidad_local1, densidades(idx_particula)
!		factor = (densidad_local1 / densidades(idx_particula))
!		if(factor >= n_vecinos) then
!			factor = 1.0D+0 * n_vecinos
!		endif
!	else
!		factor = 1.0D+0
!	endif

	aceleracion(0) = aceleracion(0) / n_vecinos
	aceleracion(1) = aceleracion(1) / n_vecinos
	aceleracion(2) = aceleracion(2) / n_vecinos
	
!	aceleracion(0) = aceleracion(0) / (n_vecinos / factor)
!	aceleracion(1) = aceleracion(1) / (n_vecinos / factor)
!	aceleracion(2) = aceleracion(2) / (n_vecinos / factor)
		
	return
	
end function

!Posiciones r1 y r2 en metros, velocidades v1 y v2 en metros por segundo, temperatura en Kelvins
!Densidades locales en gramos/cm^3
real*8 function ArtificialViscosityPoints(r1, r2, v1, v2, temperatura, densidad_local1, densidad_local2, soft_len)
	
	use Constantes
	use Vectores
	
	implicit none
	
	real*8 r1(0:2), r2(0:2), v1(0:2), v2(0:2), vij(0:2), rij(0:2)
	real*8 temperatura, cs, direccion, densidad_local1, densidad_local2, densidad_local_promedio, soft_len, miu
	
	cs = SoundSpeed(temperatura)	!velocidad de sonido en metros/segundo

	direccion = productoEscalar3D(diferenciaVectores3D(v1, v2), diferenciaVectores3D(r1, r2))
	
	!write(*,*) "En ArtificialViscosityPoints. r1 y r2 - v1 y v2: ", r1, r2, v1, v2, " densidades: ", densidad_local1, densidad_local2
		
	if(direccion < 0) then
		densidad_local_promedio = (densidad_local1 + densidad_local2) / 2.0D+0
		densidad_local_promedio = densidad_local_promedio / 1000.0D+0	!convierte densidad a kg/m^3
		miu = Viscosidad_miu(r1, r2, v1, v2, soft_len)
		ArtificialViscosityPoints = ( (-1.0D+0 * VISC_ALFA * cs * miu) + (VISC_BETA * (miu**2)) ) / densidad_local_promedio
	else
		ArtificialViscosityPoints = 0.0D+0
	endif
	
	return

end function

function GradientePresion(idx_particula, n_vecinos, lista_vecinos, presiones, N, pos_x, pos_y, pos_z, masas, densidades, soft_len) result(vector_gradiente)

	use Vectores
	use Constantes
	
	implicit none
	
	integer idx_particula, j, n_vecinos, N	
	integer lista_vecinos(0:n_vecinos-1)
	real*8 presiones(0:N-1), pos_x(0:N-1), pos_y(0:N-1), pos_z(0:N-1), masas(0:N-1), densidades(0:N-1)
	real*8 vector_gradiente(0:2), posicion1(0:2), posicion2(0:2), diferencia(0:2), unitario(0:2)
	
	real*8 magnitud_pos, soft_len, factor, diff_presion
	
	vector_gradiente = 0.0D+0

	posicion1(0) = pos_x(idx_particula) * PARSEC_MTS
	posicion1(1) = pos_y(idx_particula) * PARSEC_MTS
	posicion1(2) = pos_z(idx_particula) * PARSEC_MTS		
			
	do j = 0, n_vecinos - 1, 1		
								
		posicion2(0) = pos_x(lista_vecinos(j)) * PARSEC_MTS
		posicion2(1) = pos_y(lista_vecinos(j)) * PARSEC_MTS
		posicion2(2) = pos_z(lista_vecinos(j)) * PARSEC_MTS
		
		diferencia(0) = posicion1(0) - posicion2(0)
		diferencia(1) = posicion1(1) - posicion2(1)
		diferencia(2) = posicion1(2) - posicion2(2)

		magnitud_pos = magnitudVector3D(diferencia)

		unitario(0) = diferencia(0) / magnitud_pos
		unitario(1) = diferencia(1) / magnitud_pos
		unitario(2) = diferencia(2) / magnitud_pos
				
		!write(*,*) "En GradientePresion. Magnitud dif: ", magnitud_pos, " soft_len: ", soft_len, " - mayor: ", magnitud_pos >= soft_len

!		if(idx_particula == 100) then
!			write(*,*) "P1: ", presiones(idx_particula), " P2: ", presiones(lista_vecinos(j)), "P1 - P2: ", diff_presion, ". Unitario y posicion: ", unitario, diferencia, "mag diferencia: ", magnitud_pos, " Vector: ", vector_gradiente
!		endif
		
		diff_presion = presiones(idx_particula) - presiones(lista_vecinos(j))

		vector_gradiente(0) = vector_gradiente(0) + ( (diff_presion * unitario(0)) / magnitud_pos )
		vector_gradiente(1) = vector_gradiente(1) + ( (diff_presion * unitario(1)) / magnitud_pos )
		vector_gradiente(2) = vector_gradiente(2) + ( (diff_presion * unitario(2)) / magnitud_pos )

		if(magnitud_pos >= soft_len) then
			vector_gradiente(0) = vector_gradiente(0) + ( (diff_presion * unitario(0)) / magnitud_pos )
			vector_gradiente(1) = vector_gradiente(1) + ( (diff_presion * unitario(1)) / magnitud_pos )
			vector_gradiente(2) = vector_gradiente(2) + ( (diff_presion * unitario(2)) / magnitud_pos )
		else
			vector_gradiente(0) = vector_gradiente(0) + ( (diff_presion * unitario(0)) / soft_len )
			vector_gradiente(1) = vector_gradiente(1) + ( (diff_presion * unitario(1)) / soft_len )
			vector_gradiente(2) = vector_gradiente(2) + ( (diff_presion * unitario(2)) / soft_len )
		endif
		
	enddo
	
!	factor = 0.05

	vector_gradiente(0) = vector_gradiente(0) / (n_vecinos)
	vector_gradiente(1) = vector_gradiente(1) / (n_vecinos)
	vector_gradiente(2) = vector_gradiente(2) / (n_vecinos)
	
	return
end function

!Calcula la presion segun la ley de gas ideal
!Recibe: Densidad en g/cm^3, temperatura en Kelvin
!Resultado en pascales (Newton/m^2)
real*8 function PresionGasIdeal(densidad, temperatura)

	use Constantes
	
	implicit none
	
	real*8 densidad, temperatura, exponente

	!Debo convertir la densidad a kg/m^3

	!write(*,*) "En Fisica.PresionGasIdeal. densidad: ", densidad

	!Exponente politropico usado en <pendiente>
	
	exponente = 1.0D+0
	
	if( (densidad / 1.0E-13) <= 0.25D+0) then
		exponente = 1.0D+0
	else
	

		if( (densidad / 1.0E-13) > 0.25D+0 .and. (densidad / 1.0E-13) <= 5.0D+0) then
			exponente = 1.1D+0
		else
			if( (densidad / 1.0E-13) .gt. 5.0D+0) then
				exponente = 2.0D+0!4.0D+0 / 3.0D+0
			endif
		endif
	
	endif
	
	PresionGasIdeal = ((K_BOLTZMANN / (PESO_MOLECULAR_H * ATOMIC_MASS_CONSTANT)) * (temperatura**exponente) * ((densidad / 1000.0)))

	return

end function

real*8 function SoundSpeed(temperatura)

	use Constantes

	implicit none

	real*8 temperatura

	SoundSpeed = SQRT((K_BOLTZMANN / (PESO_MOLECULAR_H * ATOMIC_MASS_CONSTANT)) * temperatura)

end function

!Calcula el momento de inercia
!Recibe masa en kg y radio en metros
!
real*8 function momentoInercia(Masa, Radio)

	use Constantes

	implicit none
		
	real*8 Masa, Radio
	
	momentoInercia = Masa * Radio**2 * (2.0D+0 / 5.0D+0)

	return

end function

!Calcula la energia gravitacional
!Recibe masa en kg y radio en metros
!
real*8 function energiaGravitacional(Masa, Radio)

	use Constantes

	implicit none
	
	real*8 Masa, Radio
	
	energiaGravitacional = (G * Masa**2) / Radio

	return

end function

!----
real*8 function calculaVelocidadAngular(Masa, Radio, beta)

	use Constantes

	implicit none
	
	real*8 Masa, Radio, beta
	
	calculaVelocidadAngular = SQRT( (energiaGravitacional(Masa, Radio)  * 2.0D+0 * beta) / momentoInercia(Masa, Radio) )

	return

end function

!----
!Recibe aceleración en m/s^2, posiciones en metros
!Devuelve w en rad/seg

real*8 function velocidad_angular(acc_grav3D, vector_posicion3D, beta)

	use Constantes
	use Vectores
	
	implicit none
	
	real*8 beta
	real*8 acc_grav3D(0:2), vector_posicion3D(0:2)
	real*8 acc_grav2D(0:1), vector_posicion2D(0:1)
	
	acc_grav2D(0) = acc_grav3D(0)
	acc_grav2D(1) = acc_grav3D(1)
	
	vector_posicion2D(0) = vector_posicion3D(0)
	vector_posicion2D(1) = vector_posicion3D(1)
	
	velocidad_angular = SQRT( (beta * magnitudVector2D(acc_grav2D)) / magnitudVector2D(vector_posicion2D) )

	return

end function


!-----------------------------------------------
!Vector velocidad angular metros x segundo
!Recibo w en radianes/seg
!Posicion en parsecs
!
function vector_velocidad_angular(vector_posicion3D, w) result(vector_velocidad)
	use Constantes
	use Tipos
	use Vectores
	
	implicit none

	integer i, N
	
	real*8 radio_xy, circunferencia, seccion, w
	real*8 vector_posicion2D(0:1), vector_posicion3D(0:2), vector_w(0:2), vector_velocidad(0:2)
	
	vector_posicion2D(0) = vector_posicion3D(0)
	vector_posicion2D(1) = vector_posicion3D(1)
	
	radio_xy = magnitudVector2D(vector_posicion2D)
	circunferencia = 2.0D+0 * PI * radio_xy
	seccion = (w * circunferencia) / (2.0D+0 * PI)
	
!	write(*,*) "Seccion: ", seccion, " circunferencia: ", circunferencia, " radio x,y: ", radio_xy
	
	vector_w(0) = 0.0D+0
	vector_w(1) = 0.0D+0
	vector_w(2) = seccion
	
	!El producto cruz debe ser con el vector unitario del radio!
	
	vector_velocidad = productoCruz3D(vector_w,vector_posicion3D/magnitudVector3D(vector_posicion3D)) * PARSEC_MTS
	
	return
	
end function
!----------

function momento_angular_total(masas, pos_x, pos_y, pos_z, v_x, v_y, v_z, N) result(Lvect)
	use Constantes
	use Tipos
	use Vectores
	
	implicit none

	integer i, N
	
	real*8 vector_posicion3D(0:2), vector_velocidad3D(0:2) 
	real*8 Lvect(0:2), Li(0:2)
	real*8 masas(0:N-1), pos_x(0:N-1), pos_y(0:N-1), pos_z(0:N-1), v_x(0:N-1), v_y(0:N-1), v_z(0:N-1)
	
	Lvect(0) = 0.0D+0
	Lvect(1) = 0.0D+0
	Lvect(2) = 0.0D+0
	
	do i = 0, N - 1, 1
		vector_posicion3D(0) = pos_x(i)
		vector_posicion3D(1) = pos_y(i)
		vector_posicion3D(2) = pos_z(i)
		!
		vector_velocidad3D(0) = v_x(i)
		vector_velocidad3D(1) = v_y(i)
		vector_velocidad3D(2) = v_z(i)
		!
		Li = vector_momento_angular(vector_posicion3D, vector_velocidad3D, masas(i))
		
		Lvect(0) = Lvect(0) + Li(0)
		Lvect(1) = Lvect(1) + Li(1)
		Lvect(2) = Lvect(2) + Li(2)
		
	enddo		
	
	return
	
end function

!----------
!Masa en Kg, posicion en metros, velocidad en m/s
function vector_momento_angular(vector_posicion3D, vector_velocidad3D, masa) result(L)
	use Constantes
	use Tipos
	use Vectores
	
	implicit none

	integer i, N
	
	real*8 vector_posicion3D(0:2), vector_velocidad3D(0:2), L(0:2)
	real*8 masa	
		
	L = productoCruz3D(vector_posicion3D, vector_velocidad3D * masa) !!!* PARSEC_MTS
	
	return
	
end function

!Calcula el free fall time para una nube con la densidad dada en kg/m^3
real*8 function freefalltimeKgm(densidad)

	use Constantes

	implicit none
	
	real*8 densidad
	
	freefalltimeKgm = SQRT(1.0D+0 / (G * densidad) )

	return

end function

real*8 function freefalltimegcm(densidad)

	use Constantes

	implicit none
	
	real*8 densidad
	
	freefalltimegcm = SQRT(1.0D+0 / (G * 1000 * densidad) )

	return

end function
!-----------------------------------------------
!Calcula radio en funcion de masa y densidad
!Recibe masa en masas solares, y densidad en g/cm3
!Devuelve el resultado en cm
real*8 function calcularRadio(masa, densidad)
	use Constantes
	implicit none

	real*8 :: masa, densidad

	calcularRadio = ((3.0 * (masa * SOLAR_MASS_GR)) / (4.0 * PI * densidad)) ** (1.0/3.0)	
	
	return
end function
!-----------------------------------------------
!Calcula radio en funcion de masa y densidad
!Recibe masa en masas solares, y densidad en g/cm3, altura en cm
!Devuelve el resultado en cm
real*8 function calcularRadioDisco(masa, densidad, fraccion_altura)
	use Constantes
	implicit none

	real*8 :: masa, densidad, fraccion_altura

	calcularRadioDisco = ( (masa * SOLAR_MASS_GR * fraccion_altura) / (densidad * 2.0D+0 * PI) )**(1.0D+0/3.0D+0)
	
	return
end function
!-----------------------------------------------
real*8 function calcularVolumen(radio)
	use Constantes
	implicit none

	real*8 :: radio

	calcularVolumen = (4.0D+0/3.0D+0) * PI * (radio ** 3)

	return
end function
!-----------------------------------------------
real*8 function calcularVolumenDisco(radio, fraccion_altura)
	use Constantes
	implicit none

	real*8 :: radio, fraccion_altura

	calcularVolumenDisco = (2.0D+0 * PI * (radio**3)) / fraccion_altura

	return
end function
!-----------------------------------------------
!La densidad local a una particula se estima a partir de la lista de vecinos. Se suma la masa total contenida por la particula
!y sus vecinos, se calcula el volumen del paralelep�pedo que las contiene. El paralelep�pedo se define con los valores m�ximos de x, y, z de los vecinos
!Devuelve la densidad de g/cm^3
!-----------------------------------------------
real*8 function calcularDensidadLocal(N, id_particula, lista_vecinos, n_vecinos, pos_x, pos_y, pos_z, masas)
	use Tipos
	use Constantes
	implicit none

	integer N, id_particula, n_vecinos, i, j
	
	integer lista_vecinos(0:n_vecinos-1)
	
	real*8 radio, volumen, masa_total, volumen_total, max_x, max_y, max_z, min_x, min_y, min_z, resultado, lado1, lado2, lado3
	
	real*8 pos_x(0:N-1), pos_y(0:N-1), pos_z(0:N-1), masas(0:N-1)
	
	real*8 pos_x_v(0:n_vecinos-1), pos_y_v(0:n_vecinos-1), pos_z_v(0:n_vecinos-1), masas_v(0:n_vecinos-1)
	
	if(masas(id_particula) == 0.0D+0) then
		!Particulas eliminadas no tienen densidad local
		calcularDensidadLocal = 0.0D+0
		return	
	endif
	
	masa_total = masas(id_particula)
		
	volumen_total = 0.0D+0
	
	do i = 0, n_vecinos - 1, 1
		masa_total = masa_total + masas(lista_vecinos(i))
		pos_x_v(i) = pos_x(lista_vecinos(i))
		pos_y_v(i) = pos_y(lista_vecinos(i))
		pos_z_v(i) = pos_z(lista_vecinos(i))
	enddo
	
	max_x = maxval(pos_x_v)
	max_y = maxval(pos_y_v)
	max_z = maxval(pos_z_v)

	min_x = minval(pos_x_v)
	min_y = minval(pos_y_v)
	min_z = minval(pos_z_v)	

	lado1 = abs(max_x - min_x)
	lado2 = abs(max_y - min_y)
	lado3 = abs(max_z - min_z)

	volumen_total = (lado1* PARSEC_CMS) * (lado2 * PARSEC_CMS) * (lado3 * PARSEC_CMS)
	
	!write(*,*) "En fisica.calcularDensidadLocal. Maximos: ", max_x, max_y, max_z, "- Minimos: ", min_x, min_y, min_z, " - Lados: ", lado1, lado2, lado3, " - Volumen: ", volumen_total, " - Masa: ", masa_total
	
	resultado = ((masa_total * SOLAR_MASS_GR) / volumen_total) 
			
	calcularDensidadLocal = resultado
	
	return
end function
!-----------------------------------------------
!Calcula el centro de masa entre 2 particulas
function centroMasa(m1, x1, y1, z1, m2, x2, y2, z2) result(cm)
	use Constantes
	use Tipos
	implicit none

	integer i, N
	
	real*8 m1, x1, y1, z1, m2, x2, y2, z2
	
	real*8 cm(0:2)
	
	cm(0) = ((x1 * m1) + (x2 * m2)) / (m1 + m2)
	cm(1) = ((y1 * m1) + (y2 * m2)) / (m1 + m2)
	cm(2) = ((z1 * m1) + (z2 * m2)) / (m1 + m2)
	
	return
	
end function

!-----------------------------------------------
!Funcion para calcular la aceleracion de la gravedad entre 2 particulas
!Radio en metros, masa en kg, posiciones en metros
!Resultado en metros/segundo^2
!-----------------------------------------------

!Utiliza un softening de esfera uniforme para evitar aceleraciones irreales en encuentros cercanos debdos a timesteps muy grandes
!A una distancia mayor al soft_len, el potencial es Newtoniano

function aceleracion_g_vect(posicion1, posicion2, masa_total, soft_len) result(valor)
	use Constantes
	use Vectores
	implicit none

	real*8 posicion1(0:2), posicion2(0:2), valor(0:2), diferencia(0:2)

	real*8 masa_total, magnitud_pos, soft_len
	
	diferencia = diferenciaVectores3D(posicion1, posicion2)
	magnitud_pos = magnitudVector3D(diferencia)

	if(magnitud_pos >= soft_len) then
		!Se comporta como Newtoniano	
		valor(0) = (-1.0D+0 * G * masa_total * diferencia(0)) / (magnitud_pos**3)
		valor(1) = (-1.0D+0 * G * masa_total * diferencia(1)) / (magnitud_pos**3)
		valor(2) = (-1.0D+0 * G * masa_total * diferencia(2)) / (magnitud_pos**3)
	else
		!write(*,*) "Tuve que usar el softening. R: ", magnitud_pos/PARSEC_MTS, " Soft: ", soft_len/PARSEC_MTS, " ratio: ", soft_len/magnitud_pos
		valor(0) = (-1.0D+0 * G * masa_total * diferencia(0)) / (soft_len**3)
		valor(1) = (-1.0D+0 * G * masa_total * diferencia(1)) / (soft_len**3)
		valor(2) = (-1.0D+0 * G * masa_total * diferencia(2)) / (soft_len**3)	
	endif
	
	!write(*,*) "Aceleracion (x,y,z,m): ", valor, magnitudVector3D(valor)
	
	return

end function


!-----------------------------------------------
!Funcion para calcular la aceleracion de la gravedad de una distribucion de masa esferica con radio r
!Radio en metros, masa en kilogramos
!Resultado en metros/segundo^2
!-----------------------------------------------

real*8 function aceleracion_g(radio, masa)
	use Constantes
	implicit none

	real*8 radio, masa
	
!	print *, "Masa: ", masa, "distancia: ", radio
	
	aceleracion_g = -1.0D+0 * (G*masa)/(radio**2)

	return

end function

!Calcula el vector de fuerza gravitacional sobre una partícula, ejercido por N particulas en un espacio 3D

function calcularAceleracionGFB(p1, coordenadas_x, coordenadas_y, coordenadas_z, masas, N) result(fuerza_vect)
	use Constantes
	use Auxiliar
	use Tipos
	use Vectores
	use Funciones
	
	implicit none

	integer i, N
	real*8 coordenadas_x(0:N-1), coordenadas_y(0:N-1), coordenadas_z(0:N-1), masas(0:N-1)
	real*8 distancia, fuerza_esc, fuerza12_esc, macum, fuerza_ant_esc, conv_prom, acum_log
	real*8 fuerza_vect(0:2), fuerza12_vect(0:2)
	type(Particula) p1, p2

	fuerza_vect(0) = 0.0D+0
	fuerza_vect(1) = 0.0D+0
	fuerza_vect(2) = 0.0D+0
	
	fuerza_esc = 0.0D+0
	
	macum = 0.0D+0
	
	i = 0
		
	do i = 0, N - 1, 1
		if(.not. i == p1%id) then
			p2 = construirParticula(i, masas(i), coordenadas_x(i), coordenadas_y(i), coordenadas_z(i), 0.0D+0, 0.0D+0, 0.0D+0, 0.0D+0)
			fuerza12_vect = aceleracion_g_vect(p1%posicion*PARSEC_MTS, p2%posicion*PARSEC_MTS, p2%masa, 0.0D+0)
			fuerza_vect(0) = fuerza_vect(0) + fuerza12_vect(0)
			fuerza_vect(1) = fuerza_vect(1) + fuerza12_vect(1)
			fuerza_vect(2) = fuerza_vect(2) + fuerza12_vect(2)		
		endif
	enddo

end function

!Utilizando ecuación de http://www.wolframalpha.com/entities/calculators/jeans_mass
!Densidad en gramos/cm^3, temperatura en Kelvin
!Resultado en Kg
real*8 function masaJeansH2(densidad, temperatura)

	use Constantes
	
	implicit none

	real*8 densidad, temperatura
	
	masaJeansH2 = 1000.0 * (1 / SQRT(densidad)) * temperatura**(3.0/2.0) * ( ( (5*K_BOLTZMANN) / (G * MASA_H2_GR) )**(3.0/2.0) ) * SQRT(3.0 / (4.0 * PI))

	return
	
end function

function separarComponentesMomentum(r, L3, p_magnitud) result(p)
	
	use Funciones

	implicit none
	
	real*8 L3, p_magnitud, x, y, c, p1, p2
	
	real*8 r(0:2), p(0:1)
	
	x = L3 / r(0)
	
	y = r(1) / r(0)
	
	c = 1.0D+0 * (p_magnitud**2 - x**2 - (2.0D+0 * x * y))

	write(*,*) "En separarComponentesMomentum: ", r, L3, p_magnitud, x, y, c

	p = resolverCuadratica(1.0D+0, y**2, c)
	
!!!	write(*,*) "Cuadratica resuelta: ", p
	
	p1 = p(0)
	
	p2 = (L3 + (r(1)*p1)) / r(0)
	
	p(0) = p1
	p(1) = p2
	
	return

end function

END MODULE Fisica

