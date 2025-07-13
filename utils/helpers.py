import numpy as np

def newton_raphson_solver(func, jacobian, x0, tol=1e-6, max_iter=100):
    """
    Risolve un sistema di equazioni non lineari f(x) = 0 usando il metodo di Newton-Raphson.

    Args:
        func (callable): La funzione vettoriale f(x) da risolvere. Deve restituire un array 1D.
        jacobian (callable): La funzione che calcola la matrice Jacobiana J(x) di f(x).
                             Deve restituire un array 2D.
        x0 (np.ndarray): Il guess iniziale per la soluzione.
        tol (float): La tolleranza per la convergenza.
        max_iter (int): Il numero massimo di iterazioni.

    Returns:
        np.ndarray: La soluzione x che rende f(x) = 0.

    Raises:
        RuntimeError: Se il solutore non converge entro il numero massimo di iterazioni.
    """
    x = np.array(x0, dtype=float)

    for i in range(max_iter):
        f_val = func(x)
        J_val = jacobian(x)

        # Controlla la convergenza
        if np.linalg.norm(f_val) < tol:
            return x

        # Calcola il passo di Newton: J * delta_x = -f
        # np.linalg.solve è più stabile e veloce di np.linalg.inv(J_val) @ -f_val
        try:
            delta_x = np.linalg.solve(J_val, -f_val)
        except np.linalg.LinAlgError:
            raise RuntimeError("Errore di singolarità nella matrice Jacobiana.")

        x = x + delta_x

    raise RuntimeError(f"Newton-Raphson non è converso dopo {max_iter} iterazioni. Ultimo errore: {np.linalg.norm(func(x))}")

def numerical_jacobian(func, x, h=1e-6):
    """
    Calcola la matrice Jacobiana numericamente usando differenze finite centrali.

    Args:
        func (callable): La funzione vettoriale f(x) da derivare.
        x (np.ndarray): Il punto in cui calcolare la Jacobiana.
        h (float): La dimensione del passo per la differenza finita.

    Returns:
        np.ndarray: La matrice Jacobiana.
    """
    x = np.array(x, dtype=float)
    n = len(x)
    jacobian_matrix = np.zeros((n, n))
    f_x = func(x)

    for j in range(n):
        x_plus_h = np.copy(x)
        x_minus_h = np.copy(x)
        x_plus_h[j] += h
        x_minus_h[j] -= h
        jacobian_matrix[:, j] = (func(x_plus_h) - func(x_minus_h)) / (2 * h)
    return jacobian_matrix

# Esempio di utilizzo (solo per dimostrazione)
if __name__ == "__main__":
    # Esempio: risolvere f(x,y) = [x^2 + y - 1, x + y^2 - 1] = 0
    def my_func(xy):
        x, y = xy
        return np.array([x**2 + y - 1, x + y**2 - 1])

    def my_jacobian(xy):
        x, y = xy
        return np.array([[2*x, 1],
                         [1, 2*y]])

    # Guess iniziale
    x0_guess = np.array([0.5, 0.5])

    print("--- Test del solutore Newton-Raphson con Jacobiana analitica ---")
    try:
        solution_analytical = newton_raphson_solver(my_func, my_jacobian, x0_guess)
        print(f"Soluzione analitica: {solution_analytical}")
        print(f"Valore di f(soluzione): {my_func(solution_analytical)}")
    except RuntimeError as e:
        print(f"Errore: {e}")

    print("\n--- Test del solutore Newton-Raphson con Jacobiana numerica ---")
    try:
        # Usa numerical_jacobian come funzione jacobian per newton_raphson_solver
        solution_numerical = newton_raphson_solver(my_func, lambda x: numerical_jacobian(my_func, x), x0_guess)
        print(f"Soluzione numerica: {solution_numerical}")
        print(f"Valore di f(soluzione): {my_func(solution_numerical)}")
    except RuntimeError as e:
        print(f"Errore: {e}")
