sim = [[1.    ,0.921 ,0.    ,0.    ,0.    ,0.   ],
[0.921 ,1.    ,0.    ,0.    ,0.    ,0.   ],
[0.    ,0.    ,1.    ,0.003 ,0.003 ,0.003],
[0.    ,0.    ,0.003 ,1.    ,0.768 ,0.001],
[0.    ,0.    ,0.003 ,0.768 ,1.    ,0.001],
[0.    ,0.    ,0.003 ,0.001 ,0.001 ,1.   ]]

result = []
for row in sim:
    new_row = []
    for j in row:
        new_row.append(1.0-(2.0*float(j)/(1.0+float(j))) ** (1.0 / 21.0))
    result.append(new_row)

print(result)

result = []
for row in sim:
    new_row = []
    for j in row:
        new_row.append(1.0-(2.0*float(j)/(1.0+float(j))) ** (1.0 / 21.0))
    result.append(new_row)

print(result)