"""
    CRIPTOGRAFIA Y COMPUTACIÓN
    Relación 1: Aritmética modular

"""

import numpy as np
import os
import random
import math

######### Ejercicio 1 #########

def mcd_r(m):
    if m[0,0]%m[1,0]!=0:
        c=m[0,1:]-(m[0,0]//m[1,0])*m[1,1:]
        m[0,1:]=m[1,1:]
        m[1,1:]=c
        aux=m[1,0]
        m[1,0]=m[0,0]%m[1,0]
        m[0,0]=aux
        m=mcd_r(m)
    return m

def mcd(a,b):
    m=np.array([[a,1,0],[b,0,1]])
    m=mcd_r(m)
    return m[1,:]


######### Ejercicio 2 #########

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

######### Ejercicio 3 #########

# Función para calcular a elevado a b modulo n
def modPotencia(a,b,n):
    p=1
    while b>0:
        r=b%2
        if r==1:
            p=p*a%n
        a=(a*a)%n
        b=(b-r)//2
    return p

######### Ejercicio 4 #########

# Función para calcular los valores de 
# u y s para un número n-1 que es par,
# esto consiste en calcular la descomposición en 
# números primos de un número par

def descomposicionUyS(n):
    u=0
    while(modPotencia(n,1,2)==0):
        n=n//2
        u+=1
    return u,n

# Funciónes para hacer las comprobaciones pertinentes a 
# cada numero aleatorio que escogemos de Zₚ y sus potencias 
# de 2ᵘ*s, esta función se llama k veces desde "esPrimo" con 
# k dependiendo del grado de confianza que necesitemos

def listaL_apuntes(a,u,s,p):
    # Para cada a_i comprobamos unas condiciones
    # para sus potencias

    # Partimos de a^s
    a_s=modPotencia(a,s,p)
    if a_s==1 or a_s==p-1:
        return True

    # Usamos u-1 ya que ya hemos comprobado 
    # para u=0
    for k in np.arange(u-1):
        a_s=modPotencia(a_s,2,p)
        if a_s==p-1:
            return True
        if a_s==1:
            return False
    return False


# Función para determinar si dado un número p es (probablemente)
# primo usando el método de Miller-Rabin, si devolvemos que no es primo
# lo hacemos con total seguridad y si decimos que lo es lo hacemos 
# con un error de 1/4¹⁰ 

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


######### Ejercicio 5 #########

# Para el cálculo de un s que sea mayor o igual que 
# √p usamos el método de Newton-Raphson para el calculo
# de raices cuadradas aunque sabemos la función 
# math.sqrt de python tambien funciona para valores grandes

# Valor n para el que queremos la raiz cuadrada y l
# para el nivel de precisión
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


# Funcion para el algoritmo de paso enano-paso gigante
# El algoritm tiene como entrada un a y c no nulos, un p 
# primo y buscamos un b tal que a^b = c en Zₚ

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
    
    return "no existe el logaritmo"
        




    
    



###########################################################
# Main para elegir el ejercicio que queremos mostrar
###########################################################

def main():
    clear = lambda: os.system('clear')
    tecla=""
    random.seed(42)

    # Lista de primos de gran tamaño para hacer 
    # pruebas
    primos=[46381, 768479, 9476407, 36780481, 562390847, 1894083629,
    65398261921, 364879542899, 8590365927553, 28564333765949, 123456789101119,
    623084000430982607975364879776457600049]

    while(tecla!="e"):
        print("Elija un ejercicio (1-8)")
        tecla=input()
        #clear() # Limpiamos la terminal

        

        #################################################################################
        a,b=28,13
        c,d=252336560693540533935881068298825202079,38942750026936412998460304986028600003

        if tecla=='1':
            m=mcd(a,b)
            print("MCD positivo de",a,"y "+str(b)+"="+str(m[0])+" ;u:",m[1],"v:",m[2])
        if tecla=='2':
            print(modInverso(a,b))
        if tecla=='3':
            # Probamos nuestra función con números muy grandes
            print(modPotencia(c,d,primos[-1]))
        if tecla=='4':
            # Probamos miller rabin para primos muy grandes
            print("Lista de números primos:")
            for p in primos:
                res=": no es primo"
                if esPrimo(p,10):
                    res=": es primo"
                print("\t",p,res)
        
        if tecla=='5':
            a,c,p=13,2,19
            print("El logaritmo en base",a,"de",c,"modulo",p,"es: ",
                pasoEnanoGigante(a,c,p))

        if tecla=='6':
            print("Work in progress")




        #################################################################################

        tecla=input("\nCualquier tecla para volver al menu, \"e\" para salir\n")

if __name__ == "__main__":
    main()
