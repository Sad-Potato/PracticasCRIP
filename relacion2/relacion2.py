"""
    CRIPTOGRAFIA Y COMPUTACIÓN
    Relación 2: Secuencias Pseudo-aleatorias

    @author: Sergio Fernández Vela

"""

from cv2 import threshold
import numpy as np
import random

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


######### Ejercicio 2 #########

# Función auxiliar para contar los bits de un 
# entero
def sumaBitsEntero(n):
    suma=0
    while(n):
        suma+=n&1
        n=n>>1
    return suma%2


"""
    Implementar LFSR.
    La funcion recibe los coeficientes del polinomio de
    conexión, la semilla o estado inicial del LFSR y la 
    longitud de la secuencia de salida.

    ##### FORMATO DE LA ENTRADA #####
        nSalida : Entero
        coeficientes : String de bits tal que "1001"
            de mayor a menor exponente, en este caso x4 + x2 + 1
            no tengo en cuenta el término independiente
        semilla : String de bits "0110"
"""
def lfsr(nSalida, coeficientes, semilla):

    # Al tener un tamaño de salida requerido, nos da 
    # igual la periodicidad de la salida, además el número
    # de stages o L va a ser el número de elementos de la
    # semilla
    lenght=len(semilla) # Guardo el tamaño de la semilla para tener
            # constancia de los 0's a la izquierda
    actual=int(semilla,2)
    coeficientes=int(coeficientes,2)
    salida=np.empty(nSalida)

    # Para el tamaño de la salida realizamos lo siguiente
    for i in np.arange(nSalida):
        # Añadimos a la salida el último elemento del 
        # estado actual( si el valor actual no es par, acaba en 1
        # , en otro caso 0)
        salida[i]=int(actual%2!=0)

        # Antes de hacer el shift, calculamos el 
        # último valor
        calculo=coeficientes & actual

        # Hacemos un shift de 1 bit
        actual=actual>>1

        # Añadimos un 0 o un 1 segundo corresponda a partir del 
        # polinomio de conexión, para ello hago un AND entre el polinomio
        # de conexión y el estado actual.
        actual=(sumaBitsEntero(calculo)<<lenght-1) | actual

    return salida


######### Ejercicio 3 #########

"""
    La implementación es equivalente a la de LFSR pero
    en este caso tenemos una función no lineal

    ##### FORMATO DE LA ENTRADA #####

"""

def nlfsr(funcionP, semilla, k):
    
    print("TEST")

    lenght=len(semilla)
    actual=int(semilla,2)
    coeficientes=int(coeficientes,2)
    salida=np.empty(k)

    for i in np.arange(k):
        salida[i]=int(actual%2!=0)
        calculo=coeficientes & actual
        actual=actual>>1
        
        # Función no lineal
        actual=(sumaBitsEntero(calculo)<<lenght-1) | actual

    return salida

######### Ejercicio 4 #########

"""
    Generador de Geffe

    ##### FORMATO DE LA ENTRADA #####

"""

def geffe(lfsr1, lfsr2, lfsr3, salida):

    # Transformo las salidas de los lfsr para poder operar
    # con ellas con la operaciones de bits inherentes al 
    # tipo int
    # print(np.array2string(c3.astype(np.int64), max_line_width=10000).translate({32: None, 91: None, 93: None, 44: None}))

    x1=int(np.array2string(lfsr1.astype(np.int64), max_line_width=500000, threshold=500000)
        .translate({32: None, 91: None, 93: None, 44: None}),2)
    x1ex=int("".join(map(str, lfsr1.astype(int))),2)
    x2=int(np.array2string(lfsr2.astype(np.int64), max_line_width=500000, threshold=500000)
        .translate({32: None, 91: None, 93: None, 44: None}),2)
    x3=int(np.array2string(lfsr3.astype(np.int64), max_line_width=500000, threshold=500000)
        .translate({32: None, 91: None, 93: None, 44: None}),2)

    # Combinamos las salidas de los LFSR's mediante
    # la función no lineal de la definicion y las funciones de
    # bits que tiene el tipo int

    # Al usar valores de L que son primos entre ellos 
    # tenemos que el periodo del resultado será de 
    # (2^L1-1)(2^L2-1)(2^L3-1)
    
    resultado=""
    for i in np.arange(salida):
        lbit="0"
        if i%50000==0:
            print(i)

        x1i, x2i, x3i= x1 & 1, x2 & 1, x3 & 1
        
        # Toma el elemento más a la derecha y opero
        x1, x2, x3=x1>>1, x2>1, x3>>1
        nextValor=x1i & x2i ^ x2i & x3i ^ x3i

        if nextValor:
            lbit="1"
        
        resultado=resultado+lbit

    # print("2", len(resultado), resultado[-20:-2])




###########################################################
# Main para elegir el ejercicio que queremos mostrar
###########################################################

def main():
    tecla=""
    random.seed(42)

    #################################################################################
    # PARAMETROS PARA LOS EJERCICIOS
    #################################################################################

    # Secuencias de ejemplo, solo la primera de ellas cumple los 
    # postulados
    
    # Ejercicio 1
    secuencias=np.array([
        [0,0,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1],
        [0,0,1,0,0,0,1,1,1,0,1],
        [0,0,0,0,0,1,0,0,0,1,1,0,1],
        [0,1,0,0,1,1,0,1]
    ], dtype=object)

    # Ejercicio 2
    cadR1=lfsr(15, "111", "100")
    cadR2=lfsr(15, "111", "101")

    cadI1=lfsr(15, "1100", "0110")
    cadI2=lfsr(15, "1100", "1001")

    cadP1=lfsr(15, "1001", "0110")
    cadP2=lfsr(7, "101", "010")

    #################################################################################
    # Bucle para la eleccion del ejercicio
    #################################################################################
    
    while(tecla!="e"):
        print("Elija un ejercicio (1-5)")
        tecla=input()
        #clear() # Limpiamos la terminal

        #################################################################################

        if tecla=='1':
            for sec in secuencias:
                res=postuladosGolomb(np.array(sec))
                sol=""
                if not res:
                    sol="no"
                print("La cadena",np.asarray(sec),sol,"cumple los postulados")

        if tecla=='2':
            """
            ##### FORMATO DE LA ENTRADA #####
                nSalida : Entero
                coeficientes : String de bits tal que "1001"
                    de mayor a menor exponente, en este caso x4 + x2 + 1
                    no tengo en cuenta el término independiente
                semilla : String de bits "0110"

            """
            # Supongo que el grados de los polinomios es igual a 
            # L (tamaño de la semilla)
            ####### Ejemplo con polinomios reducibles en Z2 #######

            print("Polinomios reducibles")
            print("\tSecuencias para x3 + x2 + x1 + 1")
            print("\t",cadR1,"\n\t",cadR2)

            # Como vemos, para la primera cadena "100" tenemos un periodo
            # de 4, mientras que para "101" tenemos periodo 2, demostrando
            # que para polinomios reducibles el periodo depende de la semilla

            ####### Ejemplo con polinomio irreducible en Z2 #######
            # Tomo 1 polinomio irreducible, pero no primitivo

            print("Polinomios irreducibles")
            print("\tSecuencias para x4 + x3 + 1",)
            print("\t",cadI1,"\n\t",cadI2)

            # Como vemos tienen el mismo periodo, ambas tienen periodo 3, en la primera
            # comienza a repetirse desde el primer elemento mientras que en la segunda el 
            # periodo comienza desde el tercer elemento (0 1 1), por tanto demostrando la 
            # independencia de la semilla con el periodo

            ####### Ejemplos con polinomios primitivos #######
            # Tomo como ejemplo 2 polinomios primitivos sacados de 
            # la tabla, los cuales deben tener como periodo 15 y 7
            # respectivamente para unas semillas cualesquiera

            print("Polinomios primitivos")
            print("\tSecuencia para x4 + x1 + 1",cadP1)
            print("\tSecuencia para x3 + x1 + 1",cadP2)

            # Y como vemos no se repiten secuencias (si duplico el tamaño
            # de la cadena de salida obtengo 2 cadenas identicas), es decir, 
            # tienen periodo maximo

            # Ahora compruebo que las cadenas de los polinomios primitivos
            # cumplen los postulados Golomb
            if postuladosGolomb(cadP1.astype(int)) and postuladosGolomb(cadP2.astype(int)):
                print("Ambas cadenas cumplen los postulados")
            else:
                print("Alguna de las 2 cadenas o ambas no cumplen los postulados")

        if tecla=='3':
            # Función polinomica
            f=np.array([[]])

            # Semilla
            semilla="1011"

            nlfsr()

        if tecla=='4':

            a="10000000"
            n=len(a)
            # print(a, n)
            periods = [i for i in range(2,n//2+1) if a[:i]*(n//i)==a[:n - n % i]]
            if np.size(periods)==0:
                print("Tiene periodo", n, "es decir del tamaño del vector")
            else:
                print("Tiene periodo minimo de", periods[0])

            # Generador de Geffe

            # Tomamos 3 valores de L que sean primos relativos entre
            # ellos, tomo como valores 4, 5, 9 y estos serán los grados
            # de los polinomios primitivos de los 3 LFSR, los listo
            # a continuación (no se representa el termino independiente)
            # además de las semillas de tamaño L y aleatorias
            p1, p2, p3="1001", "10010", "100000001"
            s1, s2, s3="1101", "10101", "010111001"

            # Guardamos el L para cada LFSR
            l1, l2, l3=len(p1), len(p2), len(p3)

            # Periodo máximo
            periodo=((2**l1)-1)*((2**l2)-1)*((2**l3)-1)

            # Calculamos los 3 LFSR de periodo máximo
            lfsr1=lfsr(periodo,p1, s1)
            lfsr2=lfsr(periodo,p2, s2)
            lfsr3=lfsr(periodo,p3, s3)
            print("test")

            # print("1", periodo)
            geffe(lfsr1, lfsr2, lfsr3, periodo)

        #################################################################################

        tecla=input("\nCualquier tecla para volver al menu, \"e\" para salir\n")

if __name__ == "__main__":
    main()