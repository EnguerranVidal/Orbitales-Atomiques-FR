# PROJET DE CALCUL SCIENTIFIQUE / MODULE DE S4 EN LICENCE L2 PARCOURS SPECIAL
# PROGRAMME D'AFFICHAGE D'ORBITALES ATOMIQUES PAR ISOSURFACE
# Par Guillaume SIBRA et Enguerran VIDAL




###############################################################
#                           IMPORTS                           #
###############################################################

from tableau_periodique import*
import numpy as np
from skimage import measure
from plotly import tools
from plotly.offline import download_plotlyjs, plot



###############################################################
#                           CLASSES                           #
###############################################################

class Element():
    '''
    La classe Element permet de créer un Atome et donne sa configuration éléctronique.

    Paramètres :
    
    - name : nom ou symbole de l'élement
    - charge : charge de l'Atome

    Attributs :
    
    - name = nom de l'élement
    - symbol = symbole de l'élement
    - atomic_number = numéro atomique de l'atome ( Z )
    - number_electrons = nombre d'éléctrons que possède l'Atome
    - orbitals = liste des orbitales éléctroniques de l'Atome
    - valence = liste des orbitales éléctroniques de l'Atome étant de Valence
    - Slater_Zeff = liste des charges effectives dans l'approximation de Slater, respectivement aux orbitales de Valence

    Fonctions :

    - Slater_list() : permet de donner la liste des Z effectifs par la méthode de Slater
    - display_valence() : permet d'afficher les isosurfaces des orbitales de Valence ( à utiliser une fois Slater_list() utilisée sous peine d'erreurs )
    

    '''
    def __init__(self,name,charge):
        element = search_and_extract(name)# Invoque une recherche de l'Element dans la base de données à partir du nom ou symbole ( tableau périodique ) et extrait :
        self.name = element[0]            # Le nom de l'element qui est mis dans l'attribut name
        self.symbol = element[1]          # Le symbole de l'élement qui est mis dans l'attribut symbol
        self.atomic_number = element[2]   # La charge du noyau ou numéro atomique qui est mis dans l'attribut atomic_number
        self.charge = charge              # On definit la charge de l'Atome
        self.number_electrons = self.atomic_number-self.charge # on calcule le nombre d'éléctrons en soustrayant la charge de l'Atome à la charge du noyau ( égale à Z en absence de charge )
        self.orbitals = []
        self.valence = []
        self.Slater_Zeff = []
        liste_OAs = orbital_occupation(self.number_electrons) # Invocation de la fonction utilisant la methode de Klechkowsky
        n = len(liste_OAs)
        for i in range(n): # Création des OAs et remplissage de celles ci avant de les placer dans self.orbitals
            triplet = liste_OAs[i] # triplet et une liste sous la forme [n,l,nombre d'éléctrons dans l'orbitale]
            orbital = OA(triplet[0],triplet[1],triplet[2])
            self.orbitals.append(orbital)
        prim_nums = []
        for i in self.orbitals: # On prend les orbitales non remplies et de n maximal, elles deviennent de Valence
            prim_nums.append(i.prim_num)
        n = max(prim_nums)
        for i in self.orbitals:
            if i.prim_num == n or i.remplie == 0:
                i.De_Valence()
                self.valence.append(i) # On les rajoute dans self.valence


    def Slater_list(self):
        ''' Fonction utilisant la méthode de Slater afin de donner une liste de Z effectifs relatifs aux orbitales de Valence '''
        for i in range(len(self.valence)): #Traitement de toutes les orbitales de valence une à une
          Zeff = self.atomic_number #On donne une valeur de départ pour Zeff, correspondant à la charge du noyau
          ov=self.valence[i] #Traitement de la sous-couche de valence n°i, en définissant les paramètres dont nous avons besoin (n et l)
          na = ov.prim_num
          la = ov.sec_num
          if la==0 or la==1: #Traitement du cas où l'orbitale de valence est une s ou p
            for j in range(len(self.orbitals)): #On parcourt toutes les orbitales occupées de l'atome de sorte à couvrir tous les électrons
              oa=self.orbitals[j] #Comme pour l'OV, les paramètres qui nous intéressent sont n et l, et on regarde l'écran des électrons de la sous-couche n°j
              nb = oa.prim_num
              lb = oa.sec_num
              if na == nb == 1 : #série de conditions pour couvrir tous les effets d'écran possibles et appliquer le bon coefficient d'écrantage
                Zeff = Zeff -0.30*(oa.electrons - 1) #le "-1" permet d'éviter de prendre en compte l'électron dont on cherche l'effet d'écran
              if na == nb != 1 and la == lb:
                Zeff = Zeff - 0.35*(oa.electrons - 1)
              if na == nb !=1 and la != lb :
                Zeff = Zeff - 0.35*(oa.electrons)
              if na >= nb+2 :
                Zeff = Zeff-oa.electrons
              if na == nb+1:
                Zeff = Zeff-0.85*(oa.electrons)
            self.Slater_Zeff.append(Zeff) #Ajout de la charge effective de l'électron à une liste dans le même ordre que celui des OV pour l'implanter dans la fonction d'onde
          if la==2 or la == 3 : #Même chose que tout à l'heure, sauf qu'on traite le cas où l'OV est une d ou f
            assert na >= 3
            for j in range(len(self.orbitals)):
              oa=self.orbitals[j]
              nb = oa.prim_num
              lb = oa.sec_num
              if na == nb and lb == la:
                Zeff = Zeff - 0.35*(oa.electrons - 1)
              if na == nb and lb != na:
                Zeff = Zeff - 0.35*(oa.electrons)
              if na>=nb+1:
                Zeff = Zeff-oa.electrons
            self.Slater_Zeff.append(Zeff)

    def display_valence(self,valeur=20000,multiple=6):
        '''
        La fonction display_valence permet d'afficher les orbitales de Valence de l'Atome.

        Paramètres :
    
        - valeur : valeur de la Fonction d'Onde selon laquelle sera tracée les isosurfaces
        ( valeur par défaut de 20000 )
        - multiple : demi épaisseur de la largeur en rayons de Bohr ( 1 a0 = 0.529*10**(-10) m ) du repère cartésien cubiques dans lequel sera calculer la Focntion d'Onde
        ( valeur par défaut de 6 )

        '''
        n_rows = len(self.valence) # On choisit autant de lignes qu'il y a d'orbitales à afficher
        n_columns = 0
        for i in range(n_rows): # Determination du nombre de colonnes ( en mesurant le maximum de nombres de sous orbitales parmis les orbitales de Valence )
            oa = self.valence[i]
            print(oa)
            m_liste = oa.magnetic_number
            if n_columns < len(m_liste):
                n_columns = len(m_liste)
        types = []
        for i in range(n_rows): # Création de la matrice de configuration des subplots
            line = []
            for j in range(n_columns):
                line.append({'is_3d': True}) # On doit definir tout les subplots dans un tableau avec la mention {'is_3d': True}
            types.append(line)               # si jamais le graphique sera un plot 3D ( c'est notre cas )
        names = []
        for i in range(n_rows): # Création de la liste de configuration des noms des subplots ( afin de pouvoir afficher le nom des orbitales affichées )
            row=[]
            oa = self.valence[i]
            m_liste = oa.magnetic_number
            for j in range(n_columns):
                if j < len(m_liste):
                    row.append(orbitale_name(oa.prim_num,oa.sec_num,m_liste[j])) # On rajoute dans une liste 'row' les noms des sous orbitales une par une,
                else:                                                            # choisi par la fonction orbitale_name par leur n, l et m
                    row.append('')  # Si il n'y a pas d'orbitale affichée sur ce subplot, on décide de ne pas afficher de nom, mais il faut tout de même remplir la liste, d'où le ''                                    
            names = names+row # on ajoute les listes ' row' une par une dans une liste générale de noms d'orbitales qui servira peut après
        names = tuple(names) # On transforme la liste des noms d'orbitales en tuple
        print(names)
        fig = tools.make_subplots(rows = n_rows, cols = n_columns,specs = types) # Création des subplots avec la configuration créée plus haut
        annotations = tools.make_subplots(rows = n_rows, cols = n_columns,subplot_titles = names)['layout']['annotations'] # Création des noms des subplots
        fig['layout'].update(annotations=annotations)# On met à jour les noms des subplots avec la configuration 'annotations'
        for i in range(n_rows): # On parcourt les lignes du tableau de subplots
            oa = self.valence[i] 
            Zeff = self.Slater_Zeff[i]
            borne = multiple*a0()# Creation d'une variable 'borne' qui nous servira à donner les limites du mgrid que l'on créer à la ligne suivante
            X,Y,Z = np.mgrid[-borne:borne:100j, -borne:borne:100j, -borne:borne:100j] # Création d'un mgrid pour X,Y,Z ( un tableau de coordonnées de 220 points dans chaque direction
            m_liste = oa.magnetic_number #Extraction de la liste des m de l'orbitale étudiée
            for j in range(len(m_liste)): # On parcourt les colonnes du tableau pour la ligne n°i
                Fonction = Fonction_Onde(X,Y,Z,oa.prim_num,oa.sec_num,m_liste[j],Zeff)  # Calcule de la Focntion d'Onde dans l'espace ( le mgrid )
                verts,faces = measure.marching_cubes_lewiner(abs(Fonction),valeur,spacing = (X[1,0, 0]-X[0,0,0],Y[0,1, 0]-Y[0,0,0],Z[0,0, 1]-Z[0,0,0]))[:2] #Création des faces de
                # l'isosurface ainsi que des lignes entre ces faces, (l'isosurface est defini pour la valeur absolue de la fonction d'onde lorsque celle ci est égale à la variable
                # donnée ' valeur '
                verts = verts-borne # on recentre les surfaces sur l'Origine du repère cartésien de départ ( mgrid )
                data = plotly_triangular_mesh_oa(verts,faces,n=oa.prim_num,l=oa.sec_num,index=m_liste[j],Zeff=Zeff,intensities = Fonction_Onde,colorscale=color(),showscale=False)
                # On invoque alors la fonction 'plotly_triangular_mesh' qui nous crée un dictionnaire contenant toutes les informations sur la figure
                fig.add_traces(data, rows = [i+1]*len(data), cols = [j+1]*len(data))# On place la figure 'data' à sa place dans le tableau de subplots
        for i in range(len(fig['data'])): # Modification de l'opacité de toutes les surfaces ( 0: transparent , 1: opaque )
            fig['data'][i].update(opacity=0.7)
        plot(fig) # On affiche les subplots

                
                

class OA():
    '''
    La classe OA permet de créer une Orbitale Atomique

    Paramètres :
    
    - prim_num : nombre quantique principal / premier
    - sec_num : nombre quantique orbital / second / azimutale
    - electrons : nombre d'éléctrons que contient l'Orbitale

    Attributs :
    
    - prim_num : nombre quantique principal / primaire
    - sec_num : nombre quantique orbital / secondaire / azimutal
    - electrons : nombre d'éléctrons que contient l'Orbitale
    - magnetic_number : liste des nombres quantiques magnétiques pour cette Orbitale
    - remplie : permet de savoir si l'Orbitale est remplie ( 0: non remplie, 1: remplie )
    - valence : permet de savoir si l'Orbitale est de Valence ( 0 Orbitale de Coeur, 1: Orbitale de Valence )

    Fonctions :

    - De_Valence() : permet de modifier self.valence en 1
    - De_Coeur() : permet de modifier self.valence en 0
    

    '''
    def __init__(self,prim_num,sec_num,electrons):
        self.prim_num = prim_num
        self.sec_num = sec_num
        self.electrons = electrons
        self.magnetic_number = []
        
        if electrons == 2+self.sec_num*4: #On vérifie que l'Orbiatle est remplie ( 4xl + 1 éléctrons si remplie ) 
            self.remplie = 1
        else:
            self.remplie = 0
        self.valence = 0 # Le paramètre nous permettant de savoir si l'Orbitale est de Valence est mis à 0 par défaut

        for i in range(self.sec_num+1): # Création de la liste des nombres quantiques magnétiques ( liste de tout les entiers relatifs appartenant à [-l ; l]
            if i == 0:
                self.magnetic_number.append(i)
            else:
                self.magnetic_number.append(i)
                self.magnetic_number.append(-i)
        self.magnetic_number.sort() # triage de la liste des nombres quantiques magnétiques
            

    def De_Valence(self):
        self.valence = 1
        
    def De_Coeur(self):
        self.valence = 0


###############################################################
#                          FONCTIONS                          #
###############################################################

def orbital_occupation(number_electrons):
    '''
    La fonction orbital_occupation permet de donner la confiquration éléctronique d'un Atome sous forme d'une liste de triplets
    sous la forme ( nombre quantique primaire, nombres quantiques azimutal, nombres d'électrons de la couche ) .

    Paramètres :
    
    - number_electrons : nombres d'électrons que possède l'Atome étudié

    '''
    electrons_left = number_electrons
    orbitals = []
    occupation_list = [(1,0),(2,0),(2,1),(3,0),(3,1),(4,0),(3,2),(4,1),(5,0),(4,2),(5,1),(6,0),(4,3),(5,2),(6,1),(7,0),(5,3),(6,2),(7,1)]
    i = 0
    while electrons_left > 0:
        n = occupation_list[i][0]
        l = occupation_list[i][1]
        OA = [n,l]
        capacite = 2+l*4
        if electrons_left >= capacite:
            electrons_left = electrons_left-capacite
            OA.append(capacite)
        else:
            OA.append(electrons_left)
            electrons_left = 0
        orbitals.append(OA)
        i = i+1
    return orbitals

def a0():
    '''
    La fonction a0 nous permet de facilement accéder à la constante ' Rayon de Bohr'

    '''
    return 0.529*10**(-10)

def factorielle (n):
    ''' Nous permet d'accéder à la factorielle de n  '''
    s = 1
    for i in range(1,n+1):
        s = s*i
    return s   

def Laguerre(X,i,j):
    ''' Nous permet d'accéder au polynômes de Laguerre généralisés associé aux valeurs i et j puis de calculer les ordonnées pour une valeur X '''
    laguerre = [lambda x,a : 1,
              lambda x,a : -x+a+1,
              lambda x,a : (x**2)/2 -(a+2)*x+(a+2)*(a+1)/2,
              lambda x,a : -(x**3)/(6)+(a+3)*(x**2)/(2)-(a+2)*(a+3)*x/2 +(a+1)*(a+2)*(a+3)/6,
              lambda x,a : (1+a)*(2+a)*(3+a)*(4+a)/(24)-(2+a)*(3+a)*(4+a)*x/(6)+(3+a)*(4+a)*(x**2)/(4) -(4+a)*(x**3)/(6)+(x**4)/(24),
              lambda x,a : -(x**5)/(120)+(5+a)*(x**4)/(24)-(4+a)*(5+a)*(x**3)/(12)+(3+a)*(4+a)*(5+a)*(x**2)/(12)-(2+a)*(3+a)*(4+a)*(5+a)*x/(24)+(1+a)*(2+a)*(3+a)*(4+a)*(5+a)/120,
              lambda x,a : (x**6)/(720)-(6+a)*(x**5)/(120)+(6+a)*(5+a)*(x**4)/(48)-(6+a)*(5+a)*(4+a)*(x**3)/(36)+(6+a)*(5+a)*(4+a)*(3+a)*(x**2)/(48)-(6+a)*(5+a)*(4+a)*(3+a)*(2+a)*x/(120)+(6+a)*(5+a)*(4+a)*(3+a)*(2+a)*(1+a)/720 ]
    return laguerre[i](X,j)

def Fonction_Onde(X,Y,Z,n,l,index,Zeff):
    '''
    Nous rend la valeur de la Fonction d'Onde pour une Orbitale.

    Paramètres :
    
    - X : valeurs/liste/tableau des X
    - Y : valeurs/liste/tableau des Y
    - Z : valeurs/liste/tableau des Z
    - n : nombre quantique primaire
    - l : nombre quantique azimutal
    - index : valeur de l'index du m ( nombre quantique magnétique ) étudié dans la liste des m
    - Zeff : charge effective du noyau calculé par méthode de Slater

        '''
    R = np.sqrt(X**2+Y**2+Z**2)# Calcul de la norme à l'Origine
    Rad = partie_Radiale(R,Zeff,n,l) # Calcul de la partie Radiale de la Focntion d'Onde
    Ang = partie_Angulaire(X,Y,Z,R,l,index) # Calcul de la partie Radiale de la Focntion d'Onde
    return Rad*Ang   # On rend le produit des deux parties      
            
def partie_Radiale(R,Zeff,n,l):
    # Permet de calculer la partie Radiale de la Fonction d'Onde
    Rho = R*Zeff/a0()
    return np.sqrt((2*Zeff*factorielle(n-l-1))/(n*a0()*2*n*factorielle(n+l)))*(Rho**l)*Laguerre(Rho,n-l-1,2*l+1)*np.exp(-Rho/n)

def partie_Angulaire(X,Y,Z,R,l,index): # Permet de calculer la partie Angulaire de la Focntion d'Onde
    spherical = [[lambda x,y,z,r : 0.5*np.sqrt(1/np.pi)], # On cherche l'harmonique sphérique correspondant à l'orbitale dans ce tableau
                 [lambda x,y,z,r : np.sqrt(3/(4*np.pi))*(y/r),
                  lambda x,y,z,r : np.sqrt(3/(4*np.pi))*(z/r),
                  lambda x,y,z,r : np.sqrt(3/(4*np.pi))*(x/r)],
                 [lambda x,y,z,r : 0.5*np.sqrt(15/np.pi)*(x*y/(r**2)),
                  lambda x,y,z,r : 0.5*np.sqrt(15/np.pi)*(z*y/(r**2)),
                  lambda x,y,z,r : 0.25*np.sqrt(5/np.pi)*(-x**2-y**2+2*z**2)/(r**2),
                  lambda x,y,z,r : 0.5*np.sqrt(15/np.pi)*(x*z/(r**2)),
                  lambda x,y,z,r : 0.5*np.sqrt(15/np.pi)*((x**2-y**2)/(r**2))],
                 [lambda x,y,z,r : 0.25*np.sqrt(35/(2*np.pi))*((3*x**2-y**2)*y)/(r**3),
                  lambda x,y,z,r : 0.5*np.sqrt(105/np.pi)*(x*y*z)/(r**3),
                  lambda x,y,z,r : 0.25*np.sqrt(21/(2*np.pi))*(y*(4*z**2-x**2-y**2))/(r**3),
                  lambda x,y,z,r : 0.25*np.sqrt(7/np.pi)*(z*(2*z**2-3*x**2-3*y**2))/(r**3),
                  lambda x,y,z,r : 0.25*np.sqrt(21/(2*np.pi))*(x*(4*z**2-x**2-y**2))/(r**3),
                  lambda x,y,z,r : 0.25*np.sqrt(105/np.pi)*(z*(x**2-y**2))/(r**3),
                  lambda x,y,z,r : 0.25*np.sqrt(35/(2*np.pi))*(x*(x**2-3*y**2))/(r**3)]]
    return spherical[l][index](X,Y,Z,R) # On donne à la fin la valeur pour un point X,Y,Z avec une norme à l'origine de R
    
def orbitale_name(n,l,index): # Nous permet d'acceder au nom d'une orbitale atomique en la prennant dans ce tableau
    orbitale = [['s'],
                ['p {y}',
                 'p {z}',
                 'p {x}'],
                ['d {xy}',
                 'd {yz}',
                 'd {z'+chr(0x00B2)+'}',
                 'd {xz}',
                 'd {x'+chr(0x00B2)+'-y'+chr(0x00B2)+'}'],
                ['f {y(3x'+chr(0x00B2)+'-y'+chr(0x00B2)+')}',
                 'f {xyz}',
                 'f {yz'+chr(0x00B2)+'}',
                 'f {z'+chr(0x00B3)+'}',
                 'f {xz'+chr(0x00B2)+'}',
                 'f {z(x'+chr(0x00B2)+'-y'+chr(0x00B2)+')}',
                 'f {x(x'+chr(0x00B2)+'-3y'+chr(0x00B2)+')}']]
    return str(n)+orbitale[l][index]

       
def color(): # Nous permet d'acceder à un gradient de couleur pour les valeurs des isosurfaces
    pl_BtoR=[[0.0, 'rgb(165,0,38)'],
             [0.1111111111111111, 'rgb(215,48,39)'],
             [0.2222222222222222, 'rgb(244,109,67)'],
        [0.3333333333333333, 'rgb(253,174,97)'],
             [0.4444444444444444, 'rgb(254,224,144)'],
             [0.5555555555555556, 'rgb(224,243,248)'],
        [0.6666666666666666, 'rgb(171,217,233)'],
             [0.7777777777777778, 'rgb(116,173,209)'],
             [0.8888888888888888, 'rgb(69,117,180)']]
    return pl_BtoR

def intensity_func(x,y,z): # Fonction inutile servant à defaut pour la fonction 'plotly_triangular_mesh_oa' 
    return -x * np.exp(-(x**2 + y**2 + z**2))

def plotly_triangular_mesh_oa(vertices, faces,n,l,index,Zeff, intensities = intensity_func, colorscale = "Viridis", 
                           showscale = False, reversescale = False, plot_edges = False):
    '''
    La fonction plotly_triangular_mesh_oa permet de créer un dictionnaire acceptable par Plotly sous forme d'un 'mesh3d' .

    Paramètres :
    
    - vertices : un array de numpy de shape (n_vertices, 3)
    - faces = un array de numpy de shape (n_faces, 3)
    - n : nombre quantique premier
    - l : nombre quantique second
    - index : index du m ( nombre quantique magnétique ) étudié dans la liste des nombres magnétique de la classe OA()
    - Zeff : charge effective du noyau calculer par méthode de Slater
    - colorscale : gradient de couleur des surfaces
    - showscale : permet si mis sur True d'afficher l'echelle de couleur ( valeur par défaut : False )
    - reverscale : permet si mis sur True d'inverser l'echelle de couleur ( valeur par défaut : False )
    - plot_edges : permet si mis sur True d'afficher les lignes séparant les éléments de surface triangulaires ( valeur par défaut : False )
    - intensities : peut être soit une fonction f(x,y,z) ou une liste, permet de calculer une valeur afin d'associer une couleur à l'élement de surface
    ( valeur par défaut : intensity_func )
    
    '''
    x, y, z = vertices.T # Extraction des listes X,Y,Z
    I, J, K = faces.T # Extraction des coordonnées des faces 
    if hasattr(intensities, '__call__'):
        intensity = intensities(x,y,z,n,l,index,Zeff) # Si jamais on donne une fonction d'intensité
          
    elif  isinstance(intensities, (list, np.ndarray)):
        intensity = intensities # Si jamais les intensités sont données dans une liste
    else:
        raise ValueError("intensities can be either a function or a list, np.array") # erreur si aucune valeur n'a été donnée pour 'intensities' 
        
    mesh = dict(type='mesh3d',x = x,y = y,z = z,colorscale = colorscale,reversescale = reversescale,intensity = intensity,i = I,j = J,k = K,name = '',showscale = showscale)
    # Création du dictionnaire 'mesh' donnant les informations de la figure de facon à pouvoir être lu par le module Plotly
    
        
    if  showscale is True:  # Si showscale=True, on update le dictionnaire de l'Isosurface afin de montrer l'echelle de couleur
            mesh.update(colorbar = dict(thickness = 20, ticklen = 4, len = 0.75))
    
    if plot_edges is False: # Si on ne veut pas montrer les lignes separant les elements de surface 
        return  [mesh]
    else: # Si on veut montrer les lignes separant les elements de surface 
        tri_vertices = vertices[faces]
        Xe = []
        Ye = []
        Ze = []
        for T in tri_vertices:
            Xe += [T[k%3][0] for k in range(4)]+[ None]
            Ye += [T[k%3][1] for k in range(4)]+[ None]
            Ze += [T[k%3][2] for k in range(4)]+[ None]
       
        #define the lines to be plotted
        lines = dict(type = 'scatter3d',x = Xe,y = Ye,z = Ze,mode = 'lines',name = '',line = dict(color = 'rgb(70,70,70)', width = 1))
        
            
        return [mesh, lines] # On rend le mesh3d et les lignes si plot_edges=True

   



