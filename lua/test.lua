n = 100

a_table = {}
b_table = {}
c_table = {}
d_table = {}

index = {}

function to_1D_func(i, j, k, n)
    return  3*n*(i-1) + 3*(j-1) + k
end

for i = 1,n do
    a_table[i] = {}
    for j = 1,n do
        a_table[i][j] = {x=0, y=1, z=2} 
    end
end

for i = 1,n do
    b_table[i] = {}
    for j = 1,n do
        b_table[i][j] = {0,1,2} 
    end
end

for i = 1,n do
    for j = 1,n do
        for k = 1,3 do
            c_table[to_1D_func(i, j, k, n)] = k
        end
    end
end

for i = 1,n do
    d_table[i] = {}
    for j = 1,n do
        d_table[i][3 * (j-1) + 1] = 0.0 
        d_table[i][3 * (j-1) + 2] = 0.0 
        d_table[i][3 * (j-1) + 3] = 0.0 
    end
end

c_table[to_1D_func(1, 2, 3, n)] = 47

for i = 1,n do
    index[i] = math.random(1,n)
end

local start = os.clock()
for i = 1,n do
    for j = 1,n do
        local x = a_table[index[i]][index[j]].x
        local y = a_table[index[i]][index[j]].y
        local z = a_table[index[i]][index[j]].z
    end
end
print("Tabela A", -start + os.clock())

start = os.clock()
for i = 1,n do
    for j = 1,n do
        local x = b_table[index[i]][index[j]][1]
        local y = b_table[index[i]][index[j]][2]
        local z = b_table[index[i]][index[j]][3]
    end
end
print("Tabela B",- start + os.clock())

#define("to_1D(x, y, z, dimen)", "(3*(dimen*(x-1) + (y-1)) + z)")
--to_1D(42, 2) -- 42

start = os.clock()
for i = 1,n do
    for j = 1,n do
        local major = 3*(n*(index[i]-1) + (index[j]-1))
        local x = c_table[to_1D(index[i], index[j], 1, n)]
        local x = c_table[major + 2]
        local y = c_table[major + 3]
    end
end
print("Tabela C",- start + os.clock())


print(c_table[to_1D(1, 2, 3, n)])

start = os.clock()
for i = 1,n do
    for j = 1,n do
        local x = d_table[index[i]][3 * (index[j]-1) + 1]
        local y = d_table[index[i]][3 * (index[j]-1) + 2]
        local z = d_table[index[i]][3 * (index[j]-1) + 3]
    end
end

print("Tabela D",- start + os.clock())

