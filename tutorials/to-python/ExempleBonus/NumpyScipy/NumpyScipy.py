# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # MEDCoupling, NumPy et SciPy
#
# [NumPy](https://numpy.org/doc/stable/user/quickstart.html) est un package additionnel de python qui permet de manipuler des tableaux de manière optimisée. Il s’agit d’un prérequis optionnel de `MEDCoupling`.
#
# NumPy est une passerelle vers le HPC Python (multiprocessing, pyCUDA, SciPy…) qui offre de puissantes fonctions de calcul vectoriel. `MEDCoupling` est capable d'interagir avec `NumPy`.
#
# `SciPy` est aussi un package de python nécessitant `NumPy`. Il s’agit également d’un prérequis optionnel de MEDCoupling. `SciPy` offre , entre autres, des services d’algèbre linéaire et de transformé de Fourrier.
#
# Ce tutoriel montre des manières d'interagir entre les librairies `MEDCoupling`, `NumPy` et `SciPy`.

#
# Pour commencer l’exercice importer le module Python `medcoupling`:

import medcoupling as mc
from MEDCouplingRemapper import MEDCouplingRemapper
import numpy as np
import gc

# Comme `NumPy` est un prérequis optionnel, on vérifie qu'il est disponible.

assert mc.MEDCouplingHasNumPyBindings()

# ## Conversion entre DataArray et NumPy array
#
# On crée une instance de `DataArrayDouble` ayant une composante et 12 tuples. On initialise toutes les valeurs à 4.

arr = mc.DataArrayDouble(12)
arr[:] = 4.0

# On crée un tableau NumPy reposant sur les mêmes données que `arr`.

nparr = arr.toNumPyArray()
print(nparr.__repr__())
print(nparr.tolist())

# Une question intéressante est de savoir si `arr` et `nparr` partagent le même bloc mémoire.
#
# Pour le vérifier, on assigne la valeur 7.0 à un tuple sur 2 avec `nparr` et on vérifie que `arr` et `nparr` sont simultanément modifiés.

nparr[::2] = 7.0
print(nparr.__repr__())
print(arr.__repr__())

# Si on détruit `arr` (le premier à avoir alloué la mémoire), est-ce que `nparr` est détruit aussi ? Cela cause-t-il une erreur (type SIGSEGV) ?

# +
del arr

gc.collect()  # Make sure the object has been deleted
print(nparr.__repr__())
# -

# Inversement, puis-je faire une instance de `DataArrayDouble` avec `nparr` ? Oui,
# en utilisant le constructeur qui prend un `nparray` en entrée :
# contenu.:

arr2 = mc.DataArrayDouble(nparr)
print(arr2.__repr__())

# On modifie `nparr` en assignant la valeur 5.0 à tous les tuples et on vérifie que les 2 représentations ont bien été modifiées simultanément.

nparr[:] = 5.0
print(nparr.__repr__())
print(arr2.__repr__())

# Nous en profitons pour montrer un petit service pratique avec NumPy, à savoir, l’écriture optimisée. Ecrivons le contenu binaire de `nparr` dans un fichier.

f = open("toto.data", "w+b")
a = np.memmap(f, dtype="float64", mode="w+", offset=0, shape=nparr.shape)
a[:] = nparr[:]
f.flush()

# Relisons “toto.data”.

f2 = open("toto.data", "r+b")
b = np.memmap(f2, dtype="float64", mode="r", offset=0, shape=(12,))

# Pour rigoler, assignons 3.14 à `a`, flushons et relisons.

a[:] = 3.14
f.flush()
b = np.memmap(f2, dtype="float64", mode="r", offset=0, shape=(12,))
print(b.__repr__())

# On voit donc que le passage de MEDCoupling à NumPy se fait directement et de manière optimisée. Donc ca peut valoir le coup ! Tout ce qui vient d’être montré marche aussi avec des `DataArrayInt`. Regardons la suite.
#
# ## Jouons avec SciPy
#
# Nous allons créer un maillage non conforme. Le but sera de trouver la peau de ce maillage sans les surfaces non conformes.
#
# Nous allons faire cela en jouant avec les matrices creuses de SciPy (sparse matrix). Nous interpolons ce maillage non conforme sur lui même, ce qui devrait donner une matrice diagonale si le maillage était conforme.
#
# Avant nous vérifions que l’on peut jouer avec SciPy !

assert mc.MEDCouplingHasSciPyBindings()

# Pour le moment créons un maillage non conforme. Nous collons simplement deux maillages structurés avec des discrétisations spatiales différentes.:

c1 = mc.MEDCouplingCMesh()
arr1 = mc.DataArrayDouble(7)
arr1.iota()
c1.setCoords(arr1, arr1, arr1)
c2 = mc.MEDCouplingCMesh()
arr2 = mc.DataArrayDouble(9)
arr2.iota()
arr2 *= 6.0 / 8.0
c2.setCoords(arr2, arr2, arr2)

# Dégénérons `c1` et `c2` en non-structuré, une translation de `[6.,0.,0.]` de `c2`, et en faisant l’agrégation des deux, c’est dans la poche.

c1 = c1.buildUnstructured()
c2 = c2.buildUnstructured()
c2.translate([6.0, 0.0, 0.0])
c = mc.MEDCouplingUMesh.MergeUMeshes([c1, c2])

# Attention des nœuds sont dupliqués, il faut invoquer `mergeNodes()`.

c.mergeNodes(1e-12)

# Récupérons la peau et les faces non conformes. Ca nous savons faire, car nous avons fait les exercices avant :-)

skinAndNCFaces = c.computeSkin()

# Retirons les nœuds non utilisés. Cette étape n’est pas obligatoire.

skinAndNCFaces.zipCoords()

# Voici à quoi cela ressemble:
#
# <img src="skinandcells.png" style="width:500px;">
#
# OK maintenant on va séparer les cellules de bord des cellules non conformes grâce au `MEDCouplingRemapper`. Interpolons `skinAndNCFaces` sur lui-même. On acceptera un écart entre face de `1.0e-12` et un warping maximal de `0.01`.

# +
rem = MEDCouplingRemapper()
rem.setMaxDistance3DSurfIntersect(1e-12)
rem.setMinDotBtwPlane3DSurfIntersect(0.99)
rem.prepare(skinAndNCFaces, skinAndNCFaces, "P0P0")
# -

# Récupérer la matrice creuse au format CSR du remapper.

mat = rem.getCrudeCSRMatrix()

# <div class="alert alert-block alert-success">
# <b>Note:</b> Le format CSR est un format de stockage efficace des matrices creuses : Sparse matrix CSR</div>
#
# Comme nous avons bien suivi les exos sur NumPy, grâce au NumPy array `mat.indptr` on peut récupérer l’ensemble des lignes de la matrice `mat` ayant exactement un élément non nul.

indptr = mc.DataArrayInt(mat.indptr.astype(np.int64))
indptr2 = indptr.deltaShiftIndex()
cellIdsOfSkin = indptr2.findIdsEqual(1)

# C’est presque fini. Créer le sous maillage contenant uniquement la peau et l’écrire dans un fichier VTK ou MED pour le visualiser avec ParaView.

skin = skinAndNCFaces[cellIdsOfSkin]
skin.writeVTK("skin.vtu")

# <div class="alert alert-block alert-success">
# <b>Note:</b> skin contient des nœuds orphelins, on peut les retirer avec skin.zipCoords().</div>
#
# Et voilà ce que cela donne :
#
# <img src="skinonly.png" style="width:500px;">
