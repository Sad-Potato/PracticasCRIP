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
import relacion1

######### Funciones Extras #########

def mcd_r(m):
    if m[0][0]%m[1][0]!=0:
        aux=(m[0][0]//m[1][0])
        c=[m[0][1]-m[1][1]*aux,m[0][2]-m[1][2]*aux]
        m[0][1:]=m[1][1:]
        m[1][1:]=c
        aux=m[1][0]
        m[1][0]=m[0][0]%m[1][0]
        m[0][0]=aux
        m=mcd_r(m)
    return m

def mcd(a,b):
    m=[[a,1,0],[b,0,1]]
    m=mcd_r(m)
    return m[1][:]

def modPotencia(a,b,n):
    p=1
    while b>0:
        r=b%2
        if r==1:
            p=p*a%n
        a=(a*a)%n
        b=(b-r)//2
    return p

def modInverso(a,b):
    # Primero comprobamos que a y b sean primos relativos
    # es decir que el mcd de a y b sea 1

    # Calculamos con la función del apartado anterior
    # el mcd de a y b y sus coeficientes de bezout
    m=mcd(a,b)

    # Si el mcd es 1 continuamos si no, devolvemos
    # que no tiene inverso
    if m[0]==1:
        return m[1]%b
    else:
        return "Error: a no tiene inverso"

# Miller-Rabin

def descomposicionUyS(n):
    u=0
    while(modPotencia(n,1,2)==0):
        n=n//2
        u+=1
    return u,n

def listaL_apuntes(a,u,s,p):
    # Para cada a_i comprobamos unas condiciones
    # para sus potencias

    # Partimos de a^s
    aS=modPotencia(a,s,p)
    if aS==1 or aS==p-1:
        return True

    # Usamos u-1 ya que ya hemos comprobado 
    # para u=0
    for k in np.arange(u-1):
        aS=modPotencia(aS,2,p)
        # Si hay algun "a" que valga -1 devolvemos primo
        if aS==p-1:
            return True
        # Si hay algun "a" que valga 1 devolvemos que no es primo
        if aS==1:
            return False
    return False

def esPrimo(p, k):

    # Primer caso, si p es par y p ≠ 2
    if p!=2 and p%2==0:
        return False

    # Escogemos k números aleatorios que 
    # estén en Zₚ y realizamos los siguientes
    # calculos para cada uno de ellos
    u,s=descomposicionUyS(p-1)

    for i in np.arange(k):
        # Escogemos el número aleatorio
        a_i=random.randint(2,p-1)

        # Llamamos a la función que se 
        # encarga de hacer las comprobaciones
        # pertinentes para el número en Zₚ
        # Si para alguno de los números aleatorios
        # devuelve que no es primo paramos y terminamos
        # ya que tenemos la seguridad de que no lo es
        if not listaL_apuntes(a_i,u,s,p):
            return False

    # Si aguantan las condiciones devolvemos
    # que es primo
    return True

def newtonRaphson(n, l):
    # Primero hacemos una suposición con
    # respecto al valor de la raiz cuadrada de n
    raiz=n

    # Condicion, iteramos mientras el error 
    # entre iteraciones sea mayor a l
    cond=True

    while(cond):

        raiz1=0.5*(raiz+(n/raiz))

        # Si hemos llegado a una raiz con un error menor a 
        # l paramos
        if(abs(raiz1-raiz)<l):
            cond=False
        
        # Actualizamos el valor de la raiz
        # con el calculado a partir de este
        raiz=raiz1

    return raiz

def pasoEnanoGigante(a,c,p):
    
    # El primer paso es encontrar un "s" que 
    # sea mayor o igual que √p, lo calculamos 
    # con una precisión suficiente
    raizP=newtonRaphson(p,0.00000001)

    # El valor de s sera el siguiente entero
    # que exista a partir de la raiz de p
    # (s es un entero mayor o igual que raiz de p)
    s=math.ceil(raizP)
    
    # Calculamos la primera de las listas
    listaS={}
    elemento=c
    
    for r in np.arange(s): # Es decir [0,s-1]
        # Añadimos el valor y el indice a la tabla hash
        listaS[elemento]=r

        # Calculamos el valor de la 
        # siguiente iteración
        elemento=modPotencia(elemento*a,1,p)
        

    # A partir de la primera lista vamos calculando 
    # valores de la segunda hasta que encontremos 
    # una coincidencia o lleguemos al final

    # No necesitamos guardar los valores de la segunda lista
    # ya que estamos comprobandolos 1 a 1

    encontrado=False
    elemento=inicial=modPotencia(a,s,p)
    for t in np.arange(s)+1:

        # Paramos si el elemento calculado en la iteración actual
        # coincide con alguno de la primera lista
        if elemento in listaS:
            encontrado=True
            logaritmo=t*s-listaS[elemento] # r
            break

        # Calculamos el valor de la 
        # siguiente iteración
        elemento=modPotencia(elemento*inicial,1,p)

    # Si se ha encontrado alguna coincidencia entre las 2 listas
    # se devuelve el algoritmo calculado, en otro caso se vuelve
    # que no existe
    if encontrado:
        return logaritmo
    
    return None

def modPotencia(a,b,n):
    p=1
    while b>0:
        r=b%2
        if r==1:
            p=p*a%n
        a=(a*a)%n
        b=(b-r)//2
    return p

def modInverso(a,b):
    # Primero comprobamos que a y b sean primos relativos
    # es decir que el mcd de a y b sea 1

    # Calculamos con la función del apartado anterior
    # el mcd de a y b y sus coeficientes de bezout
    m=mcd(a,b)

    # Si el mcd es 1 continuamos si no, devolvemos
    # que no tiene inverso
    if m[0]==1:
        return m[1]%b
    else:
        return "Error: a no tiene inverso"

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
    # A partir de la función f creo una 
    # función h hash 
    
    # Valor de n, numéro de bits de la salida
    bits=len(bin(n)[2:])
    # print(bits)

    # Tomo como valor de b 128 bits para el tamaño 
    # del bloque
    b=128

    while(M):
        bloque=M & 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        # print(bin(bloque))
        M=M >> b

        v=funcionCompresion(bloque, v)
        # print(bloque, len(bin(bloque)[2:]))

    return v

"""     x=int("0000000100000001",2)
    c = (x >> 8) & 0xff
    f = x & 0xff
    print(c, f, "0xf") """

    # Voy tomando bloques de 
    
    

###########################################################
# Main para elegir el ejercicio que queremos mostrar
###########################################################

def main():
    tecla=""
    np.random.seed(42)
    random.seed(42)

    while(tecla!="e"):
        print("Elija un ejercicio (1-5), e para salir")
        tecla=input()
        #clear() # Limpiamos la terminal

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
            dni=54312680

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

            fechaNacimiento=20000501

            print("El inverso de",fechaNacimiento,"es",modInverso(fechaNacimiento, p))

        if tecla=='3':
            # A partir del n dado y del lemma 2.43 encontrar 2 
            # divisores no triviales de n se resume en lo siguiente
            n=48478872564493742276963

            # A partir de la información dada se 
            # deduce facilmente lo siguiente, es decir que el valor
            # de x es 12 mientras que el "y" calculado en el lemma a partir del 
            # 144 que resulta de hacer el cuadrado de 12 modulo n y n es 
            # 37659670402359614687722
            x, y=12, 37659670402359614687722

            # Y aplicando lo que se demuestra 
            # el lemma 2.43 encontramos unos valores para 
            # p y q

            if not modPotencia(x-y,1,n)==0 and not modPotencia(x+y,1,n)==0:
                div1=mcd(x-y, n)[0]
                print("Un divisor no trivial de",n ,"es", div1)

            # Calculo el otro divisor de n a partir
            # calculado
            div2=n//div1
            print("siendo el otro divisor", div2)

            # Compruebo con rabin la primalidad de ambos
            if not esPrimo(div1, 10) or not esPrimo(div2, 10):
                print("Error, uno o ambos divisores no es primo")
        
        if tecla=='4':
            # Mismo n que en el apartado anterior
            n=48478872564493742276963
            # print(len(bin(n)[2:]))

            # Tomo un valor inicial arbitrario y con 
            # como mucho los mismos bits que n
            x=n-1000000000
            # print(len(bin(x)[2:]))

            # Mensaje
            M=2342974023840238424820420834024082204272932793972
            print(len(bin(M)[2:])) # 1 bloque de 128 y otro de 33

            Ht=hashMerkleDamgard(funcionCompresion, x, M)

            print("El valor hash para M es", Ht)
            




if __name__ == "__main__":
    main()