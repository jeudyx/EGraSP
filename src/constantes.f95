MODULE Constantes
implicit none
	
	real*8, PARAMETER :: PI = 3.14159265358979323846
	real*8, PARAMETER :: SOLAR_MASS_GR = 1.98892E33			!en gramos
	real*8, PARAMETER :: SOLAR_MASS_KG = 1.98892E30			!en kilos	
	real*8, PARAMETER :: PARSEC_CMS = 3.0857E18			!en centimetros	
	real*8, PARAMETER :: PARSEC_MTS = 3.0857E16			!en metros
	real*8, PARAMETER :: UA = 149597870.7				!kilometros
	real*8, PARAMETER :: G = 6.67384E-11				!Newtons
	real*8, PARAMETER :: G_CGS = 6.67384E-8				!cgs
	real*8, PARAMETER :: SEGS_YR = 31536000.0D+0			!Cantidad de segundos en un año de 365 días
	real*8, PARAMETER :: MASA_H_GR = 1.67E-24			!Masa del hidrogeno atómico en gramos
	real*8, PARAMETER :: MASA_H2_GR = 3.9E-24			!Masa del hidrogeno molecular en gramos
	real*8, PARAMETER :: PESO_ATOMICO_H = 1.00794			!Peso atómico de hidrogeno
	real*8, PARAMETER :: PESO_MOLECULAR_H = 2.371 			!2.01594 !Peso molecular de H2
	real*8, PARAMETER :: K_BOLTZMANN = 1.3806504E-23		!Constante de Bolztmann en J/K
	real*8, PARAMETER :: ATOMIC_MASS_CONSTANT = 1.660538782E-27 	!en Kg
	real*8, PARAMETER :: CONSTANTE_GAS = 8.3144621			!Constante de gas R, valor en SI
	real*8, PARAMETER :: BIG_REAL = 99999.9999D+0
	real*8, PARAMETER :: VISC_ALFA = 0.5D+0
	real*8, PARAMETER :: VISC_BETA = 1.0D+0
	real*8, PARAMETER :: FACTOR_DISTANCIA_MAX = 3.0D+0
		
	integer, PARAMETER :: UNIT_ARBOL = 11
	integer, PARAMETER :: UNIT_NUBE = 10
	integer, PARAMETER :: NIVEL_MAX_VECINOS = 2
	integer, PARAMETER :: CANTIDAD_VECINOS = 50
	

END MODULE Constantes
