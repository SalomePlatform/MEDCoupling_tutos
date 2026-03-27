---
title: Présentation générale de MEDCoupling
author: Aymeric SONOLET
date: 08-12-2025
---

# 0. Architecture C++ / Interface Python et ouverture du code

## **Base C++**

- Performances élevées pour les opérations lourdes sur les maillages.
- Intégration naturelle avec les codes de calcul existants.
- Gestion fine de la mémoire, adaptée au HPC.

## **Interface Python**

- Utilisation quotidienne simple pour les études, le pré/post-traitement et le prototypage.
- Accès direct aux structures C++ sans perte de performance sur les calculs.
- Intégration facile dans les workflows Python déjà largement utilisés.

## **Ouverture du code**

- Code open source, facile à intégrer dans des projets variés du CEA ou en
  collaboration externe (EdF).
- Transparence sur les algorithmes, possibilité d’amélioration collaborative.
- Favorise la pérennité et la mutualisation entre équipes.

## **Intérêt dans les projets**

- À destination de développeurs de pré-pro complexes construits en python
  (comethe pour Flica5, Apollo)
- À destination de développeurs de codes de calcul pour des besoins avancés
  (TRUST, Apollo pour le couplage)
- il ne s'agit pas d'un outil disponible avec une interface graphique pour un
  utilisateur qui réalise des études.

# 1. Génération procédurale et manipulation de maillages

MEDCoupling fournit un ensemble d’outils permettant de créer, transformer et
analyser des maillages. Ces fonctionnalités couvrent les besoins classiques des
outils de pré/post-traitement pour plusieurs physiques (thermo-hydraulique,
mécanique, neutronique).

## Principales capacités

### **Construction automatique de maillages**

- maillages cartésiens et extrudés (1D, 2D, 3D) ;
- maillages non structurés (triangles, quadrangles, tétraèdres, hexaèdres, polyèdres)
- génération programmatique à partir de lignes, contours, surfaces ou volumes.

### **Opérations de modification**

- concaténation, aggrégation, conformisation
- transformations géométriques (translation, rotation, déformation).

### **Opérations géométriques avancées**

- calcul d’intersections (2D/2D, 3D/plan, etc) entre maillages ;
- localisation de points ;
- calcul de volumes, surfaces, centres de gravité.

### **Gestion des champs**

- définition de champs nodaux (P1) ou cellulaires (P0) ;
- support des champs intensifs/extensifs ;
- prise en charge du temps et de l’évolution transitoire.

### **Export direct en VTK (format standard de maillage)**

- pour la visualisation ou le post-traitement.

## Intérêt pour les projets

- Permet de remplacer des scripts hétérogènes par une API python unifiée.
- Propose des algorithmes efficaces en C++.
- Offre une base commune entre équipes simulation, CFD, mécanique, neutronique, etc.

# 2. Interface avec la librairie MED (EDF) pour la lecture/écriture

MEDLoader est un module qui repose sur la bibliothèque MED développée par EDF.
Son rôle est de fournir une couche facilitant la **conversion des données MED
vers les objets internes de MEDCoupling** et réciproquement. Il s’agit d’un
outil de traduction entre deux représentations différentes : la représentation
du format MED et la représentation interne simplifiée utilisée par MEDCoupling.

## Fonctionnalités principales

MEDLoader expose une interface permettant :

- la lecture et l’écriture de maillages, champs, groupes, familles et profils,
- la gestion du multi-pas de temps tel que défini dans le format MED,
- l’écriture incrémentale (append),
- la prise en charge des maillages multi-niveaux (maillage 3D, maillage de
  faces, maillage d’arêtes partageant les mêmes coordonnées),

Ces opérations reposent sur la conversion des données issues de MED vers les
structures MEDCoupling plus simples (DataArray, Mesh, Field), ce qui facilite
ensuite les manipulations géométriques ou les projections.

## Concepts manipulés

- **Groupes** : ensembles nommés de cellules ou entités, souvent utilisés pour
  les conditions limites.
- **Familles** : mécanisme interne du format MED permettant d’encoder la
  composition des groupes.
- **Profils** : sous-ensembles d’indices, permettant de restreindre un champ à
  une partie du maillage.

## Intérêt pour les projets

- Permet de réutiliser les structures MEDCoupling dans des workflows où le
  stockage se fait en format MED.
- Facilite la mise en place de chaînes de traitement impliquant MEDCoupling
  (manipulations géométriques, projections, transformations).
- Garantit la compatibilité avec les outils de visualisation habituels,
  notamment via l’export VTK.
- Sert d’intermédiaire lorsqu'un code ou script doit convertir des données MED
  vers un format directement exploitable dans MEDCoupling.

# 3. Couplage de codes (séquentiel et parallèle)

MEDCoupling offre des outils dédiés au **transfert de champs entre codes**, que
ce soit en mode séquentiel ou via MPI en parallèle.

Ces outils sont essentiels dans les couplages multiphysiques (neutronique -
thermique, fluide - thermique, fluide - structure, etc.), en particulier
lorsque les codes utilisent des maillages différents.

## 3.1 Couplage séquentiel (Remapper)

Le **Remapper** assure la projection de champs entre maillages, même de nature
très différente.

### Fonctionnement général

1. **Préparation** : construction de la correspondance géométrique entre
   maillage source et maillage cible.
2. **Transfert** : application d’une projection conservatrice ou non, selon la
   nature du champ.

### Choix de la nature du champ

- **Intensive** : indépendantes de la maille (température, vitesse).
- **Extensive** : proportionnelles au volume/aire (énergie, masse).
- Modes : _conservation_, _maximum_, _moyenne_ selon les besoins physiques.

### Compatibilité

- Maillages structurés, non structurés ou extrudés.
- Dimensions 1D, 2D, 3D et surfaces 2D dans un espace 3D.
- Champs P0, P1 ou points de Gauss.

## 3.2 Couplage parallèle (DEC / ParaMEDMEM)

Pour les applications HPC, MEDCoupling gère le transfert de champs entre codes
parallèles via **Data Exchange Channels (DEC)**.

### Principes de fonctionnement

- Échange des boîtes englobantes entre processus MPI.
- Détermination des sous-domaines qui se recoupent.
- Calcul des projections localement, sans reconstituer le maillage global.
- Construction distribuée de la matrice de transfert.
- Communication optimisée pour minimiser les échanges inter-processus.

### Intérêt pour les projets

- Permet le couplage performant de codes massivement parallèles.
- Évite la duplication des maillages complets sur chaque processus.
- Maintient la cohérence du transfert même sur des domaines fortement découpés.

# Conclusion

- Outil central pour plusieurs projets de simulation au CEA.
- Noyau C++ performant et interface Python favorable au développement rapide.
- Approche open source facilitant l’intégration et la collaboration.
- Large éventail de fonctionnalités : génération, manipulation et analyse de maillages.
- Socle commun pour les chaînes de pré/post-traitement.
- Outils indispensable pour le transfert de champs et le couplage
  multiphysique, en séquentiel comme en parallèle.
