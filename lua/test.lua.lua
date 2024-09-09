n = 100 -- test.lua:1
a_table = {} -- test.lua:3
b_table = {} -- test.lua:4
c_table = {} -- test.lua:5
d_table = {} -- test.lua:6
b_maj_table = {} -- test.lua:7
index = {} -- test.lua:9
to_1D_func = function(i, j, k, n) -- test.lua:11
return 3 * n * (i - 1) + 3 * (j - 1) + k -- test.lua:12
end -- test.lua:12
for i = 1, n do -- test.lua:15
a_table[i] = {} -- test.lua:16
for j = 1, n do -- test.lua:17
a_table[i][j] = { -- test.lua:18
["x"] = 0, -- test.lua:18
["y"] = 1, -- test.lua:18
["z"] = 2 -- test.lua:18
} -- test.lua:18
end -- test.lua:18
end -- test.lua:18
for i = 1, n do -- test.lua:23
b_table[i] = {} -- test.lua:24
for j = 1, n do -- test.lua:25
b_table[i][j] = { -- test.lua:26
0, -- test.lua:26
1, -- test.lua:26
2 -- test.lua:26
} -- test.lua:26
end -- test.lua:26
end -- test.lua:26
for k = 1, 3 do -- test.lua:30
b_maj_table[k] = {} -- test.lua:31
for i = 1, n do -- test.lua:32
b_maj_table[k][i] = {} -- test.lua:33
for j = 1, n do -- test.lua:34
b_table[k][i][j] = 0 -- test.lua:35
end -- test.lua:35
end -- test.lua:35
end -- test.lua:35
for i = 1, n do -- test.lua:40
for j = 1, n do -- test.lua:41
for k = 1, 3 do -- test.lua:42
c_table[to_1D_func(i, j, k, n)] = k -- test.lua:43
end -- test.lua:43
end -- test.lua:43
end -- test.lua:43
for i = 1, n do -- test.lua:48
d_table[i] = {} -- test.lua:49
for j = 1, n do -- test.lua:50
d_table[i][3 * (j - 1) + 1] = 0.0 -- test.lua:51
d_table[i][3 * (j - 1) + 2] = 0.0 -- test.lua:52
d_table[i][3 * (j - 1) + 3] = 0.0 -- test.lua:53
end -- test.lua:53
end -- test.lua:53
c_table[to_1D_func(1, 2, 3, n)] = 47 -- test.lua:57
for i = 1, n do -- test.lua:59
index[i] = math["random"](1, n) -- test.lua:60
end -- test.lua:60
print("Testes Aceso aleatorio") -- test.lua:63
local start = os["clock"]() -- test.lua:65
for oo = 1, 10000 do -- test.lua:67
for i = 1, n do -- test.lua:69
for j = 1, n do -- test.lua:70
local x = a_table[index[i]][index[j]]["x"] -- test.lua:71
local y = a_table[index[i]][index[j]]["y"] -- test.lua:72
local z = a_table[index[i]][index[j]]["z"] -- test.lua:73
end -- test.lua:73
end -- test.lua:73
end -- test.lua:73
print("Tabela A", - start + os["clock"]()) -- test.lua:78
start = os["clock"]() -- test.lua:80
for oo = 1, 10000 do -- test.lua:82
for i = 1, n do -- test.lua:84
for j = 1, n do -- test.lua:85
local x = b_table[index[i]][index[j]][1] -- test.lua:86
local y = b_table[index[i]][index[j]][2] -- test.lua:87
local z = b_table[index[i]][index[j]][3] -- test.lua:88
end -- test.lua:88
end -- test.lua:88
end -- test.lua:88
print("Tabela B", - start + os["clock"]()) -- test.lua:93
start = os["clock"]() -- test.lua:98
for oo = 1, 10000 do -- test.lua:100
for i = 1, n do -- test.lua:102
for j = 1, n do -- test.lua:103
local major = 3 * (n * (index[i] - 1) + (index[j] - 1)) -- test.lua:104
local x = c_table[(3 * (n * (index[i] - 1) + (index[j] - 1)) + 1)] -- test.lua:105
local x = c_table[major + 2] -- test.lua:106
local y = c_table[major + 3] -- test.lua:107
end -- test.lua:107
end -- test.lua:107
end -- test.lua:107
print("Tabela C", - start + os["clock"]()) -- test.lua:112
print(c_table[(3 * (n * (1 - 1) + (2 - 1)) + 3)]) -- test.lua:115
start = os["clock"]() -- test.lua:117
for oo = 1, 10000 do -- test.lua:118
for i = 1, n do -- test.lua:120
for j = 1, n do -- test.lua:121
local x = d_table[index[i]][3 * (index[j] - 1) + 1] -- test.lua:122
local y = d_table[index[i]][3 * (index[j] - 1) + 2] -- test.lua:123
local z = d_table[index[i]][3 * (index[j] - 1) + 3] -- test.lua:124
end -- test.lua:124
end -- test.lua:124
end -- test.lua:124
print("Tabela D", - start + os["clock"]()) -- test.lua:130
print("Testes Acesso sequencial") -- test.lua:132
start = os["clock"]() -- test.lua:134
for oo = 1, 10000 do -- test.lua:136
for i = 1, n do -- test.lua:138
for j = 1, n do -- test.lua:139
local x = a_table[i][j]["x"] -- test.lua:140
local y = a_table[i][j]["y"] -- test.lua:141
local z = a_table[i][j]["z"] -- test.lua:142
end -- test.lua:142
end -- test.lua:142
end -- test.lua:142
print("Tabela A", - start + os["clock"]()) -- test.lua:147
start = os["clock"]() -- test.lua:149
for oo = 1, 10000 do -- test.lua:151
for i = 1, n do -- test.lua:153
for j = 1, n do -- test.lua:154
local x = b_table[i][j][1] -- test.lua:155
local y = b_table[i][j][2] -- test.lua:156
local z = b_table[i][j][3] -- test.lua:157
end -- test.lua:157
end -- test.lua:157
end -- test.lua:157
print("Tabela B", - start + os["clock"]()) -- test.lua:162
start = os["clock"]() -- test.lua:164
for oo = 1, 10000 do -- test.lua:165
for i = 1, n do -- test.lua:166
for j = 1, n do -- test.lua:167
local x = b_table[1][j][i] -- test.lua:168
local y = b_table[2][j][i] -- test.lua:169
local z = b_table[3][j][i] -- test.lua:170
end -- test.lua:170
end -- test.lua:170
end -- test.lua:170
print("Tabela B major", - start + os["clock"]()) -- test.lua:175
start = os["clock"]() -- test.lua:180
for oo = 1, 10000 do -- test.lua:182
for i = 1, n do -- test.lua:184
for j = 1, n do -- test.lua:185
local x = c_table[(3 * (n * (i - 1) + (j - 1)) + 1)] -- test.lua:186
local x = c_table[(3 * (n * (i - 1) + (j - 1)) + 2)] -- test.lua:187
local y = c_table[(3 * (n * (i - 1) + (j - 1)) + 3)] -- test.lua:188
end -- test.lua:188
end -- test.lua:188
end -- test.lua:188
print("Tabela C", - start + os["clock"]()) -- test.lua:193
print(c_table[(3 * (n * (1 - 1) + (2 - 1)) + 3)]) -- test.lua:196
start = os["clock"]() -- test.lua:198
for oo = 1, 10000 do -- test.lua:199
for i = 1, n do -- test.lua:201
for j = 1, n do -- test.lua:202
local x = d_table[i][3 * (j - 1) + 1] -- test.lua:203
local y = d_table[i][3 * (j - 1) + 2] -- test.lua:204
local z = d_table[i][3 * (j - 1) + 3] -- test.lua:205
end -- test.lua:205
end -- test.lua:205
end -- test.lua:205
print("Tabela D", - start + os["clock"]()) -- test.lua:211
