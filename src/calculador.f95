program principal
	use Constantes
	use Vectores
	use Funciones
	use Fisica
	use Tipos 
	
	implicit none

	real*8 masa, radio, volumen, densidad, tff, acc, dt, temperatura, cs
	
	real*8 p1(0:2), p2(0:2), p3(0:2)
	
	real*8 p_vect(0:1)
	
	type(Particula) p, q
	
	radio = 0.04323276852341028D+0
	
	volumen = 4.0D+0 * PI * (radio*PARSEC_CMS)**3
	masa = 5.0D+0 * SOLAR_MASS_GR
	
	densidad = 1.27098333127372E-15
	
	temperatura = 10.0D+0
		
	write(*,*) "Masa de Jeans: ", masaJeansH2(densidad, temperatura) / SOLAR_MASS_KG
	
	write(*,*) "Sound speed: ", SoundSpeed(temperatura)
	
	return
	
	cs = 166.0D+0	!velocidad del sonido en m/s
	
	write(*,*) "Presión con cs: ", cs**2 * (densidad / 1000.0D+0)
	write(*,*) "Presion con formula gas ideal: ", PresionGasIdeal(densidad, temperatura)
	write(*,*) "Miu despejado: ", (CONSTANTE_GAS * temperatura) / (cs**2)
	
	return
	
!!!	write(*,*) "Densidad: ", densidad, " gr/cm^3" 
		
	tff = SQRT(1.0D+0 / (G * 1000.0 * densidad) )
	
!!!	write(*,*) "Tff: ", tff, tff/SEGS_YR, " densidad: ", densidad
	
	masa = 5.0D+0 * SOLAR_MASS_KG
	radio = 0.04323276852341028 * PARSEC_MTS
	
	acc = (G * masa) / radio**2
	
!	write(*,*) "Aceleracion simple en R: ", acc

	p1(0) = 0.0125! -0.004113631208071325
	p1(1) = 0.0D+0!-0.0015364663865785467
	p1(2) = 0.0D+0!0.00680956379129703
	
	p2(0) = -0.013670176364150187!-0.004120780457042294
	p2(1) = -0.01493865657866434!-0.0015384446848856045
	p2(2) = 0.017203382326520585! 0.006823801362471346
	
!!!	write(*,*) "Distancia (43): " , distanciaPuntos(p1, p2)

	radio = distanciaPuntos(p1, p2) * PARSEC_MTS
	masa = (5.0D+0 + 0.0010123766676284724) * SOLAR_MASS_KG
	
	dt = SEGS_YR
	
	acc = (G * masa) / radio**2
!!!	write(*,*) "Aceleracion simple en P: ", acc, "Velocidad en medio yr: ", acc * dt, " m/s"
	
	radio = calcularRadio(0.0010123766676284724D+0, densidad) / PARSEC_CMS
	
!!!	write(*,*) "Radio de particula peq:  ", radio
	
	radio = calcularRadio(5.0D+0, 1.0D-18) / PARSEC_CMS
	
!!!	write(*,*) "Radio de nube gde:  ", radio

	radio = calcularRadio(5.0D+0, 1.0D-10) / PARSEC_CMS
	
!!	write(*,*) "Radio de particula pesada:  ", radio
	
	!p = construirParticula(0, 5.0D+0, p1(0), p1(1), p1(2), 0.0D+0, 0.0D+0, 0.0D+0, densidad)
	!q = construirParticula(0, 0.0010123766676284724D+0, p2(0), p2(1), p2(2), 0.0D+0, 0.0D+0, 0.0D+0, densidad)
	
	!write(*,*) "Colision? : ", colision(p, q)
	
	!r_vector
	p1(0) = 2.944207432196885E+14
	p1(1) = 1.427990562574884E+14
	p1(2) = 0

	!p_vector
	p2(0) = -1.1323E+30
	p2(1) = 2.3543E+30
	p2(2) = 0

	!L_vector	
	p3 = productoCruz3D(p1, p2)
	
	write(*,*) "Resultado de producto cruz: ", p3, " - ", magnitudVector3D(p3), " ---- ", magnitudVector3D(p1) * magnitudVector3D(p2)
	
	write(*,*) "Magnitudes: ", magnitudVector3D(p3) / magnitudVector3D(p1), " --- ", magnitudVector3D(p2)
	
	write(*,*) "Vector recuperado: ", (magnitudVector3D(p3) / magnitudVector3D(p1)), " --- ", magnitudVector3D(p2)
	
	p_vect = separarComponentesMomentum(p1, p3(2), magnitudVector3D(p2))
	
	write(*,*) "Solucion de ecuacion: ", p_vect, " -- ", magnitudVector2D(p_vect)
end	