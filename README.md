# CalibGEOXIM
Code R : calibration d'un simulateur d'écoulement 3D (GEOXIM) à l'aide de données synthétiques.

# Codes R et descriptions

- **AnalysePostCalib** : Analyse et évaluation des résultats de calibration.
- **AnalyseSimulation** : Analyse des simulations du code de calcul, avec représentation graphique des moyennes et des variations dans le temps et l'espace.
- **Dinit** : Sélection du plan initial à l'aide du critère **Maxmin** et des simulations correspondantes.
- **GPinitial** : Construction du modèle initial de processus gaussien.
- **GPseq** : Amélioration du modèle gaussien initial par planification séquentielle (EI+Var).
- **MCMCsampling** : Échantillonnage de la densité a posteriori par MCMC de type Metropolis dans Gibbs.
- **OptimTempsdeplacements** : Optimisation du temps de déplacement des capteurs mobiles.
- **OptimXobs** : Optimisation des positions des puits d'observations.
- **PlotMCMC** : Visualisation des ACF, élagage des chaînes et tracé des résultats MCMC.
- **PlotMeanVarspatialettemporel** : Graphiques des moyennes et des variances en fonction de l'espace et du temps.
- **PlotQ2models** : Tracé des coefficients de prédictivité des modèles de processus gaussien.
- **PlotSVDFinit** : Tracé des fonctions propres et des valeurs propres issues de la SVD appliquée sur **Finit**.
- **PlotYtop** : Visualisation des observations physiques après optimisation du temps de déplacement des capteurs.
- **TraitementSimulationBrut** : Sélection et save des simulations brutes dans une matrice.

