from Fonction import*
import numpy as np
import random
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')


def RSA(func,Solution,nb_solu,taille_solu,itération):
    T = itération
    taille_solution = taille_solu
    nb_solution = nb_solu
    Matrice_Solution=np.copy(Solution)

    # ******** Trouver the best Solution and the worst *************
    idxMin = np.argmin(Matrice_Solution[:,taille_solution])
    solution_best_in_itra=np.copy(Matrice_Solution[idxMin,:])
    vecteur_Fonc_objectiv_meillure=np.full(shape=(T+1), fill_value=Matrice_Solution[idxMin,taille_solution],dtype=float)
    _,UB,_,LB=func(([1,2]))
    #******Inialiser les parametre utilise dans algorithme RSA******
    t = 1
    eps = 0.1
    B = 0.005
    a = 0.1
    #*******************The beginning of the RSA*******************
    while t <= T:
        ESt=2*random.uniform(0,1)*(1-(t/T))
        for i in range(nb_solution):
            r2 = random.randint(0,nb_solution-1)
            r1 = random.randint(0,nb_solution-1)
            while  r1== idxMin or r2== idxMin:
                r2 = random.randint(0,nb_solution-1)
                r1 = random.randint(0,nb_solution-1)
            rand = random.uniform(0,1)
            P = a +(Matrice_Solution[i,:-1] - (np.mean(Matrice_Solution[i,:-1])))/(solution_best_in_itra[:-1]*(UB - LB) + eps)
            N = solution_best_in_itra[:-1] * P
            R = (solution_best_in_itra[:-1] - Matrice_Solution[r2,:-1]) / (solution_best_in_itra[:-1] + eps)
            if t <= (T/4):
                Matrice_Solution[i,:-1] = solution_best_in_itra[:-1] - N * B - R * rand
            if t > (T/4) and t <= (T/2) :
                Matrice_Solution [i,:-1] = solution_best_in_itra[:-1] * Matrice_Solution[r1,:-1] * ESt * rand
            if t > (T/2) and t <= ((3*T)/4)  :
                Matrice_Solution[i,:-1] = solution_best_in_itra[:-1] * P * rand
            if t > ((3*T)/4):
                Matrice_Solution[i,:-1] = solution_best_in_itra[:-1] - N * eps - R * rand
            for j in range(taille_solution):
                if Matrice_Solution[i,j] < LB:
                    Matrice_Solution[i,j] = "{0:.3f}".format(random.uniform(LB,0))
                if Matrice_Solution[i,j] > UB:
                    Matrice_Solution[i,j] = "{0:.3f}".format(random.uniform(0,UB))
            # ***********************calcule la fonction objective***************
            Matrice_Solution[i,taille_solution],_,_,_=func(Matrice_Solution[i,:-1])
        # **********************************Fin de traitment RSA****************************************
        idxMin = np.argmin(Matrice_Solution[:,taille_solution])
        if (solution_best_in_itra[taille_solution]> Matrice_Solution[idxMin,taille_solution]):
            solution_best_in_itra=np.copy(Matrice_Solution[idxMin,:])
        vecteur_Fonc_objectiv_meillure[t] = solution_best_in_itra[taille_solution]

            
        t = t + 1

    return(vecteur_Fonc_objectiv_meillure)


T=100
n=10
fonction=Quartic
_,UB,taille_solution,LB=fonction(([1,2]))
nb_solution=n
# crer matrice du Solution
M1 = np.full(shape=(nb_solution, taille_solution+1), fill_value=0,dtype=np.float64)
# inialiser le Matrice du solution random
for j in range(nb_solution):
    for i in range(taille_solution):
        M1[j,i] = "{0:.3f}".format(random.uniform(LB, UB)) 
# calcule la fonction objective 
for j in range(nb_solution):
    M1[j,taille_solution],_,_,_=fonction(M1[j,:-1])
vecteur_best=RSA(fonction,M1,nb_solution,taille_solution,T)
ax1 =plt.axes()
ax1.plot(vecteur_best,label="RSA",c="#4B2991")
ax1.set_ylabel('Coût F')
ax1.set_xlabel("Nombre d'itération")
if fonction== Ackley or fonction == Michalewicz or fonction== Schwefel1 or fonction== Schwefel2_26 or fonction== Himmelblau:
    ax1.legend()
    ax1.set_title("Le Coût f en fonction de Nombre d'itération")
    plt.show()
else:
    ax1.set_yscale('log')
    ax1.legend()
    ax1.set_title("Le Coût f en fonction de Nombre d'itération")
    plt.show()