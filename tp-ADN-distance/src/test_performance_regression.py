import numpy as np
from sklearn.linear_model import LinearRegression


"""
# Données d'entraînement

#itératif

X_temps = np.array([10000,20000,30000,40000]).reshape(-1, 1)  # Variables indépendantes (caractéristiques)
y_temps = np.array([2.81769, 11.3032, 25.2021, 44.3899])  # Variable dépendante (cible)
X_energie = np.array([10000,20000,30000,40000]).reshape(-1, 1)  # Variables indépendantes (caractéristiques)
y_energie = np.array([1.96746e-05,7.32657e-05,0.0001414,0.0002387])
# Créez un modèle de régression linéaire
model_temps = LinearRegression()
model_energie = LinearRegression()
# Ajustez le modèle aux données d'entraînement
model_temps.fit(X_temps, y_temps)
model_energie.fit(X_energie,y_energie)

# Faites des prédictions
X_new = np.array([ 20236404]).reshape(-1, 1)  # Nouvelles données d'entrée pour les prédictions
y_temps_pred = model_temps.predict(X_new)
y_energie_pred= model_energie.predict(X_new)
print("Prédictions temps :", y_temps_pred)
print("Prédictions energie : ",y_energie_pred)
# temps: ~ 7.5heures 
# énergie : 0.14669348 kWh
"""

"""
# Données d'entraînement

#cache aware

X_temps = np.array([10000,20000,30000,40000]).reshape(-1, 1)  # Variables indépendantes (caractéristiques)
y_temps = np.array([3.74311, 16.9161, 36.237, 61.6594])  # Variable dépendante (cible)
X_energie = np.array([10000,20000,30000,40000]).reshape(-1, 1)  # Variables indépendantes (caractéristiques)
Y_energie = np.array([2.04487e-05,9.26128e-05,0.0002049,0.0003539])

# Créez un modèle de régression linéaire
model_temps = LinearRegression()
model_energie = LinearRegression()
# Ajustez le modèle aux données d'entraînement
model_temps.fit(X_temps, y_temps)
model_energie.fit(X_energie,Y_energie)


# Faites des prédictions
X_new = np.array([ 20236404]).reshape(-1, 1)  # Nouvelles données d'entrée pour les prédictions
y_temps_pred = model_temps.predict(X_new)
y_energie_pred= model_energie.predict(X_new)
print("Prédictions temps :", y_temps_pred)
print("Prédictions energie : ",y_energie_pred)
# temps ~ 11heures
# energie : 0.22504835 kWh
"""


"""
#cache oblivious

# Données d'entraînement

X_temps = np.array([10000,20000,30000,40000]).reshape(-1, 1)  # Variables indépendantes (caractéristiques)
y_temps = np.array([8.96335, 35.3481, 79.7995, 142.422])  # Variable dépendante (cible)
X_energie = np.array([10000,20000,30000,40000]).reshape(-1, 1)  # Variables indépendantes (caractéristiques)
Y_energie= np.array([5.28187e-05,0.0001909,0.0003930,0.0007156])

# Créez un modèle de régression linéaire
model_temps = LinearRegression()
model_energie = LinearRegression()
# Ajustez le modèle aux données d'entraînement
model_temps.fit(X_temps, y_temps)
model_energie.fit(X_energie,Y_energie)


# Faites des prédictions
X_new = np.array([ 20236404]).reshape(-1, 1)  # Nouvelles données d'entrée pour les prédictions
y_temps_pred = model_temps.predict(X_new)
y_energie_pred= model_energie.predict(X_new)
print("Prédictions temps :", y_temps_pred)
print("Prédictions energie : ",y_energie_pred)


# temps : ~ 25heures
# energie : 0.44305755 kWh
"""


