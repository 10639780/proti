x = [i for i in range(5)]
y = [5] * 5

ziplist = list(zip(x, y))
print(ziplist)

if (1, 5) in ziplist:
    print("found coordinates")