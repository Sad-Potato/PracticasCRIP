"""
    CRIPTOGRAFIA Y COMPUTACIÓN
    Relación 1: Aritmética modular

"""
import numpy as np
import os

from sqlalchemy import false

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
        a=pow(a,2)%n
        b=(b-r)/2
    return p

######### Ejercicio 4 #########

# Función para calcular los valores de 
# u y s para un número n-1 que es par,
# esto consiste en calcular la descomposición en 
# números primos de un número par

def descomposicionUyS(n):
    u=0
    while(n%2==0):
        n/=2
        u+=1
    return u,int(n)

# Funciónes para hacer las comprobaciones pertinentes a 
# cada numero aleatorio que escogemos de Zₚ y sus potencias 
# de 2ᵘ*s, esta función se llama k veces desde "esPrimo" con 
# k dependiendo del grado de confianza que necesitemos

# Función de los apuntes de clase, aqui 
# la probabilidad esta en que sea primo

def listaL_wiki(a,u,s,p):
    # Para cada a_i comprobamos unas condiciones
    # para sus potencias

    aPot=modPotencia(a,s,p)
    """ print(aPot) """
    if aPot==1 or aPot==p-1:
        return True

    for k in np.arange(u-1):
        aPot=modPotencia(aPot,2,p)
        """ print(u,aPot) """
        if aPot==p-1:
            return True
    return False




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
        if a_s!=1 and a_s!=p-1 and modPotencia(a_s,2,p)==1:
            return False
        if a_s==1:
            return False
    return False



# Función de "Lecture Notes on Cryptography", la
# probabilidad esta en que no sea primo

# ESTA MAL?¿?¿?

def listaL_lnc(a,u,s,p):
    # Para cada a_i comprobamos unas condiciones
    # para sus potencias

    if modPotencia(a,p-1,p)!=1:
        return False
    for k in np.arange(u)+1:
        # Calculamos a^2^(k)*s mod p con 
        # 0<=k<=u
        aModP=modPotencia(a,pow(2,k)*s,p)
        k_1=modPotencia(a,pow(2,k-1)*s,p)
        if aModP==1 and k_1!=1 and k_1!=-1:
            return False
    return True

# Función para determinar si dado un número p es (probablemente)
# primo usando el método Miller-Rabin, si devolvemos que no es primo
# lo hacemos con total seguridad y si decimos que lo es lo hacemos 
# con un error de 1/4¹⁰ 

def esPrimo(p, k):

    # Primer caso, si p es par y ≠ a 2
    if p!=2 and p%2==0:
        return False

    # Escogemos k números aleatorios que 
    # estén en Zₚ y realizamos los siguientes
    # calculos para cada uno de ellos
    u,s=descomposicionUyS(p-1)

    for i in np.arange(k):
        # Escogemos el número aleatorio
        a_i=np.random.randint(2,p-1)

        # Llamamos a la función que se 
        # encarga de hacer las comprobaciones
        # pertinentes para el número en Zₚ
        # Si para alguno de los números aleatorios
        # devuelve que no es primo paramos y terminamos
        # ya que tenemos la seguridad de que no lo es
        """ print("*****") """
        if not listaL_apuntes(a_i,u,s,p):
            return False

    # Devolvemos si es o no primo
    return True


######### Ejercicio 5 #########





###########################################################
# Main para elegir el ejercicio que queremos mostrar
###########################################################

def main():
    clear = lambda: os.system('clear')
    tecla=""
    np.random.seed(42)

    # Lista de primos de gran tamaño para hacer 
    # pruebas
    primos=[46381, 768479, 9476407, 36780481, 562390847, 1894083629,
    65398261921, 364879542899, 8590365927553, 28564333765949, 123456789101119]

    while(tecla!="e"):
        print("Elija un ejercicio (1-8)")
        tecla=input()
        clear() # Limpiamos la terminal

        a,b=28,13

        # Serie de "ifs" para cada uno de los ejercicios
        if tecla=='1':
            m=mcd(a,b)
            print("MCD positivo de",a,"y "+str(b)+"="+str(m[0])+" ;u:",m[1],"v:",m[2])
        if tecla=='2':
            print(modInverso(a,b))
        if tecla=='3':
            print(modPotencia(2,13,5))
        if tecla=='4':
            for p in primos:
                print(esPrimo(1569,10))


        tecla=input("\nCualquier tecla para volver al menu, \"e\" para salir\n")

if __name__ == "__main__":
    main()
