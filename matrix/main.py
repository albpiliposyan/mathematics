from 




matrix = [[0, 1, 1], 
          [1, 0, 4],
          [9, 3, 12]]


l = [[0] * 3, [0] * 3, [0] * 3]
u = [[0] * 3, [0] * 3, [0] * 3]

n = 3
sum = 0 
for i in range(n):
    for j in range(i + 1):
        sum = 0
        for k in range(j):
            sum += l[i][k] * u[k][j]
        l[i][j] = matrix[i][j] - sum
    u[i][i] = 1
    j = i + 1
    while j < n:
        sum = 0
        for k in range(i):
            sum += l[i][k] * u[k][j]
        u[i][j] = (matrix[i][j] - sum) / l[i][i]

        j += 1
    
    print("L:")
    for r in l:
        print(r)

    print("U:")
    for r in u:
        print(r)

