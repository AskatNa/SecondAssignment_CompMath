def cramer_method(matrix, results):
    def determinant(mat):
        if len(mat) == 2:
            return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]
        return sum((-1) ** i * mat[0][i] * determinant([row[:i] + row[i+1:] for row in mat[1:]]) for i in range(len(mat)))

    detMain = determinant(matrix)

    return [determinant([row[:i] + [results[j]] + row[i+1:] for j, row in enumerate(matrix)]) / detMain for i in range(len(matrix))]
def gauss_method(matrix, results):
    matrix = [row[:] for row in matrix]
    results = results[:]
    n = len(matrix)
    for i in range(n):
        for j in range(i, n):
            matrix[i][j] /= matrix[i][i]
        results[i] /= matrix[i][i]
        for k in range(i + 1, n):
            factor = matrix[k][i]
            for j in range(i, n):
                matrix[k][j] -= factor * matrix[i][j]
            results[k] -= factor * results[i]
    return [results[i] - sum(matrix[i][j] * results[j] for j in range(i + 1, n)) for i in range(n-1, -1, -1)]

def jacobi_method(matrix, results, MaxIter=100, tol=1e-6):
    n = len(matrix)
    x = [0] * n
    for _ in range(MaxIter):
        new_x = [(results[i] - sum(matrix[i][j] * x[j] for j in range(n) if j != i)) / matrix[i][i] for i in range(n)]
        if all(abs(new_x[i] - x[i]) < tol for i in range(n)):
            return new_x
        x = new_x
    return x
def gauss_seidel_method(matrix, results, maxIter=100, tol=1e-6):
    n = len(matrix)
    x = [0] * n
    for _ in range(maxIter):
        x_new = x[:]
        for i in range(n):
            x_new[i] = (results[i] - sum(matrix[i][j] * x_new[j] for j in range(n) if j != i)) / matrix[i][i]
        if all(abs(x_new[i] - x[i]) < tol for i in range(n)):
            return x_new
        x = x_new
    return x

matrix = [
    [10, -2, 67, 23],
    [9, 13, 20, 15],
    [54, 25, 12, -10],
    [19, 46, -52, 8]
]
results = [18, 26, 34, 82]

print("Cramer:", cramer_method(matrix, results))
print("Gauss:", gauss_method(matrix, results))
print("Jacobi:", jacobi_method(matrix, results))
print("Gauss-Seidel:", gauss_seidel_method(matrix, results))
