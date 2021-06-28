# Orbitales Atomiques ( 2018 )
 Ce projet fait en collaboration avec Guillaume Sibra fut realisé dans le cadre d'un projet du module de deuxième année de Licence de Calcul Scientifique. Il propose un affichage d'isosurfaces des orbitales atomiques de valence d'un élément donné. On utilsie pour cela les expressions analytiques des orbitales atomiques par séries d'harmoniques sphériques.

Ce programme utilise la version de Python 3.6

**CONTENU DU DOSSIER :**

Le dossier qui vous est présenté contient plusieurs éléments du projet :
- un fichier   mendeleiv.db          	 : base de données contenant notre tableau périodique
- un fichier   main_program.py    	 : fichier Python à lancer sur le Shell afin d'accéder à l'affichage d'Orbitales Atomiques
- un fichier   tableau_periodique.py 	 : fichier Python qui contient les fonctions permettant d'utiliser la base de données

**MODULES A INSTALLER :**
Voici la liste des modules Python à installer par différentes méthodes afin de pouvoir faire fonctionner le programme :
- Skymage (utiliser plutôt scikit-image)
- Numpy
- Plotly

**MODE D'EMPLOI :**

Pour afficher les Orbitales Atomiques  de Valence d'un atome, rien de plus simple, il suffit de rentrer les lignes suivantes dans un Shell Python une fois 
que les modules mentionnés plus haut on été installés.

- Tout d'abord rentrer la ligne  " X=Element('Fluor',0)" pour initialiser un atome de Fluor sans charge par exemple.

- Il faut ensuite rentrer la ligne  " X.Slater_list() " afin que l'ordinateur puisse calculer les charges effectives du noyau atomique par méthode de Slater. Ces charges seront utiliser plus tard pour le calcul de la Fonction d'Onde .

-  Il faut ensuite rentrer la ligne  " Xe.display_valence(valeur=7000,mutliple=6) " pour afficher des isosurfaces prises pour une valeur de 700 dans une grille de coordonnées cubique de demi-largeur de 6 Rayon de Bohr ( Rayon de Bohr = 0.529*10**(-10) m )

Le Shell risque de répondre avec une erreur du type ' Memory Error ', cela vient de la taille de la matrice étant alors trop volumineuse en données pour le module numpy, il vaut alors mieux jouer avec les valeurs des variables 'valeur' et 'multiple' en entrée de la fonction 'display_valence()'. C'est aussi tout à fait normal de finir avec un résultat graphique peu convaincant ( portions d'isosurfaces coupées, isosurfaces de couleurs blanches ou mal séparées ou même présentant des patterns grossiers ...) et il faut de nouveau pouvoir jouer avec les variables d'entrée
de la fonction d'affichage.
