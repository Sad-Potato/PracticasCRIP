"""
    CRIPTOGRAFIA Y COMPUTACIÓN
    Relación 3: Funciones de un solo sentido

    @author: Sergio Fernández Vela

"""

import numpy as np
import random
import os
import math
import random
import hashlib
import relacion1 as rel1 # Funciones de la primera relación

######### Ejercicio 1 #########

"""
    Función para obtener la llave pública a partir
    de una secuencia super-creciente, un n y un u
"""

def llavPublica(secuencia, n, u):
    return [modPotencia(u*secuencia[i],1,n) for i in np.arange(np.size(secuencia))]



"""

    La función de cifrado de la mochila toma el mensaje y hace la 
    sumatoria del mensaje con la llave

"""

def mochilaCipher(mensaje, llavePublica):
    return np.sum([mensaje[i]*llavePublica[i] for i in np.arange(np.size(mensaje))])


"""

    La función de descifrado toma el mensaje cifrado y la llave 
    privada (la cadena super-creciente, el u y el n)

"""

def mochilaDecipher(mensajeCifrado, cadenaA, u, n):
    # Obtengo el inverso de u módulo n, para deshacer
    # la multiplicacion de u con la cadena supercreciente
    uInv=modInverso(u, n)

    # Ahora obtengo el b de la mochila, al multiplicar por
    # el inverso de u obtenemos de forma efectiva la cadena 
    # super-creciente ya que el mensaje es una cadena binaria 
    b=modPotencia(mensajeCifrado*uInv, 1, n)

    mensajeDescifrado=np.zeros((np.size(cadenaA)), dtype=int)

    # Al quedarnos la secuencia super-creciente basta con aplicar un 
    # greedy de la mochila al valor b
    suma=0
    for i in np.arange(np.size(cadenaA))[::-1]:
        if suma+cadenaA[i]<=b:
            mensajeDescifrado[i]=1
            suma+=cadenaA[i]
    
    return mensajeDescifrado

######### Ejercicio 2 #########

""" Desarrollo en el main """

######### Ejercicio 3 #########

""" Desarrollo en el main """

######### Ejercicio 4 #########

# Mismo valor para n que en el
# ejercicio 3
n=48478872564493742276963

# Funcion h de compresión definida
# por Goldwasser, Micali y Rivest, 
# con a0 y a1 cuadrados arbitrarios
def funcionCompresion(b, x):
    a0, a1 = 144, 169

    return modPotencia(
                modPotencia(x, 2, n)*
                modPotencia(a0, b ,n)*
                modPotencia(a1, 1-b, n),
                1, n)
    
# Uso la construccion de Merkle-Damgàrd para 
# implementar una función hash sin colisiones.
# - f es una función de compresión resistente a colisiones,
# en nuestro caso la dada por el enunciado
# - v es el vector inicial con los mismos bits
# que n
# - M es el mensaje 
def hashMerkleDamgard(f, v, M):

    b=1 # Tamaño de bloque
    Msize=len(M) # Tamaño del mensaje
    M=int(M, 2)

    for i in np.arange(Msize):
        bloque=M & 0x1
        M=M >> b
        v=f(bloque, v)

    return v


######### Ejercicio 5 #########

def funcionRSA(n, e, x):
    return modPotencia(x, e, n)

######### Ejercicio 6 #########

# Factores primos p y q a partir del 
# método explicado en Notes on Cryptography,
# asumimos con este método que n es el 
# producto de 2 primos, y sabemos e y d inversos
# uno del otro
    
def primeFactors(N, d, e):
    
    # Descomponemos en potencias de 
    # 2 el resultado de d * e - 1
    a, b=rel1.descomposicionUyS(d*e-1)

    # Escogemos un valor de x para el cual
    # se cumpla que mcd(N, x)=1
    x=random.randint(1,n-1)
    while rel1.mcd(x, N)[0]!=1:
        x=random.randint(1,n-1)

    # Calculamos valor de y a partir de 
    # x y b 
    y=rel1.modPotencia(x, b, N)

    # Si se cumple esto en una primera
    # instancia el algoritmo falla
    if rel1.modPotencia(y+1, 1, N)==0 or\
        rel1.modPotencia(y-1, 1, N)==0:
        return None 


    # Buscamos "y" tal que se cumpla lo anterior
    z=y
    y=rel1.modPotencia(y, 2, N)
    while not rel1.modPotencia(y-1, 1, N)==0 and\
        not rel1.modPotencia(y+1, 1, N)==0:
        z=y
        y=rel1.modPotencia(y, 2, N)

    # Si se ha cumplido que y es congruente con -1
    # entonces el algoritmo falla, mientras que si 
    # lo ha sido con 1 podemos calcular los factores primos 
    # de N
    if rel1.modPotencia(y+1, 1, N)==0:
        return None
    
    # Hemos encontrado un z con el que calcular
    # los factores primos
    return rel1.mcd(N, z+1)[0], rel1.mcd(N, z-1)[0]

######### Ejercicio 7 #########

# Función para generar las claves, uso las funciones
# de los apartados anteriores
def keygen():
    # Escojo 2 número primos aleatorios de gran tamaño y con 
    # una lontitud de bits parecida, en concreto escojo 2 valores entre
    # 90 y 100 bits de longitud [90, 100)

    p=random.randint(618970019642690137449562112, 633825300114114700748351602688-1)
    while not rel1.esPrimo(p, 10):
        p+=1
    
    q=random.randint(618970019642690137449562112, 633825300114114700748351602688-1)
    while not rel1.esPrimo(q, 10):
        q+=1

    n=p*q

    # Calculo un valor e para la llave pública 
    # menor que el valor de phi de n
    phi_n=(p-1)*(q-1)
    e=2
    while e<phi_n and rel1.mcd(e, phi_n)[0]!=1:
        e+=1
    
    if e>=phi_n:
        print("Error")
        return None
    
    # Calculo un valor de d, siendo d 
    # el inverso de e módulo phi(n)
    d=rel1.mcd(e, phi_n)[1]

    # Devuelvo la llave pública y privada
    return [n, e], [n, d]

def firmagen(mensaje, clavePrivada):
    # Leo el mensaje y la clave privada


    hash_object = hashlib.sha1(b)
    print(hash_object.digest())

###########################################################
# Main para elegir el ejercicio que queremos mostrar
###########################################################

# TODO Tamaño de bloque del ejercicio 4
# TODO Arreglar formato ejercicio 2
# TODO Funciones relación 1 
# TODO Explicación ejercicio 6

def main():
    tecla=""
    np.random.seed(42)
    random.seed(42)

    while(tecla!="e"):
        print("Elija un ejercicio (1-7), e para salir")
        tecla=input()
        #clear() # Limpiamos la terminal

        #################################################################################
        # PARAMETROS PARA LOS EJERCICIOS
        #################################################################################
        
        n1=48478872564493742276963
        dni=54312680
        fechaNacimiento=20000501


        #################################################################################

        if tecla=='1':

            print("Introduce un valor de k")
            k=int(input())

            # Mensaje de tamaño k
            mensaje=np.array([np.random.randint(0,2) for i in np.arange(k)])
            # mensaje=np.array([0,1,1,0,0,1,0,1]) # Letra e en binario
            print("Mensaje original: ",mensaje)

            # Elijo secuencia super-creciente de tamaño
            # k, por ejemplo las k primeras potencias de 2 
            secuenciaA=[pow(2,k1) for k1 in np.arange(k)]

            # Eligo un n acorde a la secuencia super-creciente
            n=np.sum(secuenciaA)+1

            # Escojo un u tal que u y v sean coprimos
            u=5039
            while mcd(u,n)[0]!=1:
                u+=1
            
            ##########################################
            # CALCULO LA LLAVE PÚBLICA 
            ##########################################
            llavePublica=llavPublica(secuenciaA, n, u)

            ##########################################
            # FUNCION DE ENCRIPTACIÓN
            ##########################################
            mensajeCifrado=mochilaCipher(mensaje, llavePublica)
            print("Mensaje cifrado: ", mensajeCifrado)

            ##########################################
            # FUNCION DE LA MOCHILA O DESENCRIPTACION
            ##########################################
            mensajeDescifrado=mochilaDecipher(mensajeCifrado, secuenciaA, u, n)
            print("Mensaje descifrado: ", mensajeDescifrado,"\n\n")
            
        if tecla=='2':
            

            # Primo p mayor o igual que el número de 
            # identidad y que (p-1)/2 tambien lo sea para 
            # facilitar el proceso
            p=dni
            while not esPrimo(p, 10) or not esPrimo((p-1)/2, 10):
                p+=1

            # Encuentro el valor de la función phi de euler
            # para p y al ser p primo esto es p-1 (2.100)
            phi=p-1
            
            # Usamos el algoritmo de paso enano-paso gigante
            # para determinar si para un α cualquier el menor t
            # tal que a^t≡1 mod n es igual al valor de phi de p
            """ for i in np.arange(21):
                print(i, ":",pasoEnanoGigante(i, 1, 21)) """
            
            alpha=np.random.randint(0,p)
            while pasoEnanoGigante(alpha, 1, p)!=phi:
                alpha=np.random.randint(0,p)

            print(alpha)

            """ i=3
            primos=[2]
            while i<phi//2:
                if esPrimo(i, 10) and phi%i==0:
                    primos.append(i)
                i+=1

            for x in primos:
                if modPotencia(alpha, phi/x, p)-1!=0:
                    print(x, phi, modPotencia(alpha, phi/x, p))
                    print("No es elemento primitivo") """

            

            print("El inverso de",fechaNacimiento,"es",modInverso(fechaNacimiento, p))

        if tecla=='3':
            # A partir del n dado y del lemma 2.43 encontrar 2 
            # divisores no triviales de n se resume en lo siguiente

            # A partir de la información dada se 
            # deduce facilmente lo siguiente, es decir que el valor
            # de x es 12 mientras que el "y" calculado en el lemma a partir del 
            # 144 que resulta de hacer el cuadrado de 12 modulo n y n es 
            # 37659670402359614687722
            x, y=12, 37659670402359614687722

            # Y aplicando lo que se demuestra 
            # el lemma 2.43 encontramos unos valores para 
            # p y q

            if not modPotencia(x-y,1,n1)==0 and not modPotencia(x+y,1,n1)==0:
                div1=mcd(x-y, n1)[0]
                print("Un divisor no trivial de",n1 ,"es", div1)

            # Calculo el otro divisor de n a partir
            # calculado
            div2=n1//div1
            print("siendo el otro divisor", div2)

            # Compruebo con rabin la primalidad de ambos
            if not esPrimo(div1, 10) or not esPrimo(div2, 10):
                print("Error, uno o ambos divisores no es primo")
        
        if tecla=='4':
            # print(len(bin(n)[2:]))

            # Tomo un valor inicial arbitrario y con 
            # como mucho los mismos bits que n
            x=123456
            # print(len(bin(x)[2:]))

            # Mensaje
            M="00000001001" # Considero que el mensaje se 
            # lee de izquierda a derecha, en caso contrario
            # bastaría con invertir el mensaje

            Ht=hashMerkleDamgard(funcionCompresion, x, M)

            print("El valor hash para M es", Ht)
        
        if tecla=='5':
            # Encuentro p menor primo mayor o igual que mi dni
            # y q menor primo mayor o igual que mi fecha de 
            # nacimiento
            p=dni
            q=fechaNacimiento
            while not esPrimo(p, 10):
                p+=1
            while not esPrimo(q, 10):
                q+=1

            # Selecciono "e" tal que se cumple que
            # gcd(e, (p−1)(q−1)) = 1
            e=12379
            assert mcd(e, (p-1)*(q-1))[0]==1

            # n como ya se ha mencionado si no se 
            # dice lo contrario es igual a p*q
            n=p*q

            # Funcion RSA para un m igual a 10
            m=10
            print("Valor de la función RSA para m =", 10, ":", funcionRSA(n, e, 10))

            # Inverso para 1234567890
            print("El inverso para 1234567890 es", rel1.modInverso(1234567890, n))
            # print(modPotencia(1234567890*636007744202166, 1, n))

        if tecla=='6':
            n2=50000000385000000551

            # Sabemos la llave pública (n, e)
            # módulo y exponente de cifrado
            e=5

            # Sabemos la llave privada (n, d)
            # módulo y exponente de descifrado
            d=10000000074000000101
            
            # Calculamos p y q
            p, q = primeFactors(n2, d, e)
            
            if p*q==n2:
                print("Los valores de p y q para n son",p,"y",q)

            # TODO
            # El procedimiento de este ejercicio y el 3 son similares en que
            # 
        
        if tecla=='7':
            # Implemento sistema de firma RSA

            ##########################################
            # GENERACIÓN DE CLAVES
            ##########################################

            # Genero la clave publica y privada, las cuales 
            # son de la forma [n, e], [n, d]
            pubKey, privKey=keygen()

            # Guardo la llave privada en un fichero 
            f=open("./privKey.txt","w+")
            f.write(str(privKey))
            f.close()

            ##########################################
            # GENERACIÓN DE FIRMA
            ##########################################

            fichero="./mensaje.txt"
            ficheroLlavePrivada="./privKey.txt"

            firmaFichero=firmagen(fichero, ficheroLlavePrivada)

            f=open("./firmaFichero.txt","w+")
            f.write(str(firmaFichero))
            f.close()

            # Hago uso de SHA1 para generar un resumen del mensaje

            ##########################################
            # VERIFICACIÓN DE FIRMA
            ##########################################
            print()


if __name__ == "__main__":
    main()