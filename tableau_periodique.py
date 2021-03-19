# Ceci est le fichier .py permettant de commander la base de données SQL
# nommée mendeleiv.db que vous devriez normalement trouver dans le dossier
# du projet. Seules les fonctions ask_element() et search_and_extract sont utiliser
# au niveau du programme principal. Les autres n'ayant servi qu'à remplir la base de
# données et verifier son contenu.




###############################################################
#                           IMPORTS                           #
###############################################################

import sqlite3


###############################################################
#                          FONCTIONS                          #
###############################################################

#Connection à la base de données ( permet des tests des focntions ),
#donc inutiels pour le reste

def ask_element():#Permet de demander les infos necessaires au rajout d'un atome dans le tabelau periodique
    element=[]
    name=str(input('Nom ='))
    if test_presence(name)==True:
        print('Cet élément est déjà dans le Tableau')
        return '?'
    else :
        symbol=str(input('Symbole ='))
        num=int(input('Numéro Atomique ='))
        element.append(name)
        element.append(symbol)
        element.append(num)
        return element

#Commande typique SQL    
content = """CREATE TABLE tableau ( name VARCHAR(30),
symbol VARCHAR(4), atomic_number INTEGER);"""

def adding_elements():# Nous permet de rajouter plusieurs élements à la suite dans la base de données après avoir demander les infos nécessaires à l'utilisateur
    New_elements=[]
    continuer=input('Voulez vous rajoutez un élément au tableau ? (O/N)')
    while continuer=='o' or continuer=='O':
        element=ask_element()
        if element=='?'or element=='!':
            print('Erreur')
        else:
            New_elements.append(element)   
        continuer=input('Voulez vous rajoutez un autre element au tableau ? (O/N)')
        
    n=len(New_elements)
    add_element(New_elements)
    print(n,' elements ont été rajoutés.')

def add_element(elements): # Nous permet de rajouter les listes d'élements dans la base de données une fois la liste d'infos donnée
    connection = sqlite3.connect("mendeleiv.db")
    cursor = connection.cursor()
    n=len(elements)
    for i in range(n):
        element=elements[i]
        format_str = """INSERT INTO tableau (name, symbol, atomic_number)
        VALUES ('{name}', '{symbol}', '{atomic_number}');"""
        sql_command = format_str.format(name=element[0], symbol=element[1], atomic_number=element[2])
        print(sql_command)
        cursor.execute(sql_command)
    connection.commit()
    connection.close()

def extract_element(): # Nous permet d'extraire toutes les données du tableau périodique
    connection = sqlite3.connect("mendeleiv.db")
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM tableau") 
    result = cursor.fetchall() 
    connection.commit()
    connection.close()
    return result

def test_presence(element): # Nous permet de savoir si un élement est déjà présent dans le tableau ou non
    Catalogue=extract_element()
    Names=list_name(Catalogue)
    if element in Names :
        return True
    else :
        return False


def list_name(Catalogue): # Nous permet d'accéder à la liste des noms et symboles des élements du tableau
    n=len(Catalogue)
    Names=[]
    for i in range(n):
        name=Catalogue[i][0]
        symbol=Catalogue[i][1]
        Names.append(name)
        Names.append(symbol)
    return Names

def search_and_extract(element): # Nous permet de chercher un élement dans le tableau et d'en extraire ces données.
    Catalogue=extract_element()
    n=len(Catalogue)
    k='?'
    for i in range(n):
        name=Catalogue[i][0]
        symbol=Catalogue[i][1]
        if name==element or symbol==element:
            return Catalogue[i]
    return k
