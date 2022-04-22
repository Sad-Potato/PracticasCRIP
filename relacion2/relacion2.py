"""
    CRIPTOGRAFIA Y COMPUTACIÓN
    Relación 2: Secuencias Pseudo-aleatorias

    @author: Sergio Fernández Vela

"""

import numpy as np
import random

from sympy import sec

######### Ejercicio 1 #########
"""

Enunciado: Ejercicio 1. Escribe una función que determine si 
una secuencia de bits cumple los postulados de Golomb

"""
def postuladosGolomb(secuencia):

    # Comprobamos cada uno de los 3 postulados con la 
    # secuencia proporcionada

    ##############################################

    # Primer postulado: que haya el mismo número de 0s y 1s, 
    # y en caso de que tenga un número impar de elementos
    # habrá a lo mucho 1 elemento de diferencia entre ambas
    size=np.size(secuencia)
    unos=np.count_nonzero(secuencia)
    ceros=size-unos

    # Depende de si tenemos un número par o impar de elementos
    # en la cadena, tambien compruebo que no tengamos cero 1s o cero 0s
    # para evitar problemas en los siguientes postulados
    if (size%2==0 and ceros!=unos) or ceros==0 or unos==0:
        return False
    elif ceros!=unos-1 and ceros!=unos+1:
        return False

    ##############################################
    
    # Segundo postulado: teniendo en cuenta que es una secuencia
    # cíclica el número de secuencias de 0s y 1s debe ser igual, 
    # además el número de secuencias de un tamaño debe de ser la mitad
    # cuando incrementamos la longitud de la secuencia en 1 elemento

    # Antes de ver si se cumple el segundo postulado, al tratarse de 
    # de secuencias cíclicas trasladamos la secuencia hasta que el primer
    # y último elemento sean distintos
    while secuencia[0]==secuencia[-1]:
        secuencia=np.roll(secuencia,1)


    # Implementación eficiente para encontrar y contar secuencias 
    # sacada de 
    # https://stackoverflow.com/questions/58221268/count-number-of-repeated-elements-in-a-row-in-a-numpy-array

    # En primer lugar lo que hace en esta implementación es para cada elemento ver si es igual
    # a su siguiente, obteniendo un array con [ True, ..., ¿n!=n+1?, ..., True], con esto pretendemos
    # tener los indices en los que cambian los elementos
    lista = np.r_[True,secuencia[:-1]!=secuencia[1:],True]

    # Aqui se calcula el número de elementos de cada secuencia, mediante la 
    # diferencia de la lista de los elementos que son True(elementos que son distintos de su siguiente)
    numElementos = np.diff(np.flatnonzero(lista))

    # Obtenemos a que elemento corresponde cada secuencia, tomando
    # los elementos de la secuencia para los que la lista calculada
    # es True
    elemento = secuencia[lista[:-1]]

    # Contamos rachas 

    # Array en el para la columna i tenemos las rachas de 
    # tamaño i+1 tanto para 0 y 1 en los indices 0 y 1
    rachas=np.zeros((2,np.max(numElementos)),dtype=np.int64)

    # Rellenamos el array a partir de la secuencia
    for i in np.arange(np.size(numElementos)):
        rachas[elemento[i]][numElementos[i]-1]+=1

    # Ahora comprobamos que se cumple el segundo postulado
    # a partir de lo calculado
    
    suma=rachas[:][0]+rachas[:][1]
    
    # No puede haber más de 2 tamaños de rachas con solo
    # 1 racha, en cuyo caso no se cumple el postulado
    if np.size(suma[np.where(suma==1)])>2:
        return False

    # No puede haber mas rachas de tamaño n que de n-1, 
    # La suma de los elementos de la diferencia del vector suma
    # que son mayores que 0 debe ser 0
    if np.sum(np.diff(suma)>0)!=0:
        return False

    # Compruebo que haya el mismo número de cadenas longitud n
    # para 0s y 1s (cuando no hay 1 sola cadena), que se cumpla
    # la mitad de cadenas etc
    for i in np.arange(np.size(suma)):
        if i==0:
            if suma[i]!=1 and rachas[0][0]!=rachas[1][0]:
                return False
        else:
            if suma[i]!=1 and ( rachas[0][0]!=rachas[1][0] or suma[i]!=suma[i-1]/2):
                return False
            if suma[i]==1 and suma[i]!=suma[i-1]/2 and i!=np.size(suma)-1:
                return False

    ##############################################

    # Tercer postulado: Se cumple si para todos los desplazamientos
    # de la secuencia el número de elementos con respecto a la secuencia
    # original se mantiene constante

    secuenciaDesp=np.roll(secuencia,1)
    desp=np.sum(secuenciaDesp!=secuencia)
    while(not np.array_equal(np.roll(secuenciaDesp,1), secuencia)):
        secuenciaDesp=np.roll(secuenciaDesp,1)
        if desp!=np.sum(secuenciaDesp!=secuencia):
            return False
        desp=np.sum(secuenciaDesp!=secuencia)

    # Si se cumplen los 3 postulados devolvemos true
    return True


###########################################################
# Main para elegir el ejercicio que queremos mostrar
###########################################################

def main():
    tecla=""
    random.seed(42)

    #################################################################################
    # PARAMETROS PARA LOS EJERCICIOS
    #################################################################################



    #################################################################################
    # Bucle para la eleccion del ejercicio
    #################################################################################
    
    while(tecla!="e"):
        print("Elija un ejercicio (1-5)")
        tecla=input()
        #clear() # Limpiamos la terminal

        #################################################################################

        if tecla=='1':
            # Tomo como secuencia binaria la usada como
            # ejemplo en la referencia del ejercicio
            secuencia=np.array([0,0,0,1,0,0,1,1,0,1,0,1,1,1,1])

            secuencias=np.array([
                [0,0,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1],
                [0,0,1,0,0,0,1,1,1,0,1],
                [0,0,0,0,0,1,0,0,0,1,1,0,1],
                [0,1,0,0,1,1,0,1]
            ])
            for sec in secuencias:
                print(postuladosGolomb(np.array(sec)))

            res=postuladosGolomb(secuencia)
            sol=""
            if not res:
                sol="no"
            print("La cadena",secuencia,sol,"cumple los postulados")


        #################################################################################

        tecla=input("\nCualquier tecla para volver al menu, \"e\" para salir\n")

if __name__ == "__main__":
    main()