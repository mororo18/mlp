dofile("Data.lua") -- main.canc.lua:1
s_print = function(solut) -- main.canc.lua:3
for i = 1, # solut["s"] do -- main.canc.lua:4
io["write"](solut["s"][i], " ") -- main.canc.lua:5
end -- main.canc.lua:5
print() -- main.canc.lua:7
end -- main.canc.lua:7
table_print = function(tbl) -- main.canc.lua:10
for i = 1, # tbl do -- main.canc.lua:11
io["write"](tbl[i], " ") -- main.canc.lua:12
end -- main.canc.lua:12
print() -- main.canc.lua:14
end -- main.canc.lua:14
matrix_print = function(info) -- main.canc.lua:17
for i = 1, # info["c"] do -- main.canc.lua:18
for j = 1, # info["c"] do -- main.canc.lua:19
io["write"](info["c"][i][j], " ") -- main.canc.lua:20
end -- main.canc.lua:20
print() -- main.canc.lua:22
end -- main.canc.lua:22
end -- main.canc.lua:22
seq_print = function(solut, info) -- main.canc.lua:27
for i = 1, info["dimension"] + 1 do -- main.canc.lua:37
for j = i, info["dimension"] + 1 do -- main.canc.lua:38
for k = 1, 3 do -- main.canc.lua:39
io["write"](solut["seq"][(3 * ((info["dimension"] + 1) * (i - 1) + (j - 1)) + k)], " ") -- main.canc.lua:40
end -- main.canc.lua:40
io["write"]("| ") -- main.canc.lua:42
end -- main.canc.lua:42
print() -- main.canc.lua:44
end -- main.canc.lua:44
end -- main.canc.lua:44
table["clone"] = function(org) -- main.canc.lua:48
if type(jit) == "table" then -- main.canc.lua:49
return { unpack(org) } -- main.canc.lua:51
else -- main.canc.lua:51
return { table["unpack"](org) } -- main.canc.lua:53
end -- main.canc.lua:53
end -- main.canc.lua:53
solut_clone = function(solut) -- main.canc.lua:58
local cpy = {} -- main.canc.lua:59
cpy["seq"] = table["clone"](solut["seq"]) -- main.canc.lua:61
cpy["s"] = table["clone"](solut["s"]) -- main.canc.lua:62
cpy["cost"] = solut["cost"] -- main.canc.lua:63
return cpy -- main.canc.lua:65
end -- main.canc.lua:65
subseq_fill = function(seq, info) -- main.canc.lua:70
for i = 1, info["dimension"] + 1 do -- main.canc.lua:71
for j = i, info["dimension"] + 1 do -- main.canc.lua:72
for k = 1, 3 do -- main.canc.lua:73
seq[(3 * ((info["dimension"] + 1) * (i - 1) + (j - 1)) + k)] = 0.0 -- main.canc.lua:74
end -- main.canc.lua:74
end -- main.canc.lua:74
end -- main.canc.lua:74
end -- main.canc.lua:74
subseq_load = function(solut, info) -- main.canc.lua:80
local dimen = info["dimension"] -- main.canc.lua:81
for i = 1, dimen + 1 do -- main.canc.lua:82
local k = 1 - i - (i ~= 1 and 0 or 1) -- main.canc.lua:83
solut["seq"][(3 * ((dimen + 1) * (i - 1) + (i - 1)) + info["T"])] = 0.0 -- main.canc.lua:85
solut["seq"][(3 * ((dimen + 1) * (i - 1) + (i - 1)) + info["C"])] = 0.0 -- main.canc.lua:86
solut["seq"][(3 * ((dimen + 1) * (i - 1) + (i - 1)) + info["W"])] = (i ~= 1 and 1 or 0) -- main.canc.lua:87
for j = i + 1, dimen + 1 do -- main.canc.lua:88
local j_prev = j - 1 -- main.canc.lua:89
solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["T"])] = info["c"][solut["s"][j_prev]][solut["s"][j]] + solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j_prev - 1)) + info["T"])] -- main.canc.lua:91
solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["C"])] = solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["T"])] + solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j_prev - 1)) + info["C"])] -- main.canc.lua:92
solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["W"])] = j + k -- main.canc.lua:93
end -- main.canc.lua:93
end -- main.canc.lua:93
solut["cost"] = solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (dimen + 1 - 1)) + info["C"])] - info["EPSILON"] -- main.canc.lua:97
end -- main.canc.lua:97
sort = function(arr, r, info) -- main.canc.lua:101
for i = 1, # arr do -- main.canc.lua:102
for j = 1, # arr - i do -- main.canc.lua:103
if info["c"][r][arr[j]] > info["c"][r][arr[j + 1]] then -- main.canc.lua:104
local tmp = arr[j] -- main.canc.lua:105
arr[j] = arr[j + 1] -- main.canc.lua:106
arr[j + 1] = tmp -- main.canc.lua:107
end -- main.canc.lua:107
end -- main.canc.lua:107
end -- main.canc.lua:107
end -- main.canc.lua:107
construction = function(alpha, info) -- main.canc.lua:113
local s = { 1 } -- main.canc.lua:114
local cList = {} -- main.canc.lua:116
for i = 2, info["dimension"] do -- main.canc.lua:117
cList[# cList + 1] = i -- main.canc.lua:118
end -- main.canc.lua:118
local r = 1 -- main.canc.lua:123
while # cList > 0 do -- main.canc.lua:126
sort(cList, r, info) -- main.canc.lua:129
local range = math["floor"](# cList * alpha) + 1 -- main.canc.lua:132
local i = math["random"](1, range) -- main.canc.lua:133
i = info["rnd"][info["rnd_index"]] + 1 -- main.canc.lua:134
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:135
local c = table["remove"](cList, i) -- main.canc.lua:137
table["insert"](s, c) -- main.canc.lua:138
r = c -- main.canc.lua:140
end -- main.canc.lua:140
table["insert"](s, 1) -- main.canc.lua:142
return table["clone"](s) -- main.canc.lua:144
end -- main.canc.lua:144
swap = function(s, i, j) -- main.canc.lua:147
s[j], s[i] = s[i], s[j] -- main.canc.lua:148
end -- main.canc.lua:148
reverse = function(s, i, j) -- main.canc.lua:151
local l = j -- main.canc.lua:152
for k = i, math["floor"]((j + i) / 2) do -- main.canc.lua:153
swap(s, k, l) -- main.canc.lua:155
l = l - 1 -- main.canc.lua:156
end -- main.canc.lua:156
end -- main.canc.lua:156
reinsert = function(s, i, j, pos) -- main.canc.lua:160
if i < pos then -- main.canc.lua:161
for k = i, j do -- main.canc.lua:162
local e = s[i] -- main.canc.lua:163
table["insert"](s, pos, e) -- main.canc.lua:164
table["remove"](s, i) -- main.canc.lua:165
end -- main.canc.lua:165
else -- main.canc.lua:165
for k = i, j do -- main.canc.lua:168
local e = s[j] -- main.canc.lua:169
table["remove"](s, j) -- main.canc.lua:170
table["insert"](s, pos, e) -- main.canc.lua:171
end -- main.canc.lua:171
end -- main.canc.lua:171
end -- main.canc.lua:171
search_swap = function(solut, info) -- main.canc.lua:176
local cost_best = math["huge"] -- main.canc.lua:177
local I = - 1 -- main.canc.lua:178
local J = - 1 -- main.canc.lua:179
local dimen = info["dimension"] -- main.canc.lua:180
local cost_concat_1 = 0.0 -- main.canc.lua:182
local cost_concat_2 = 0.0 -- main.canc.lua:183
local cost_concat_3 = 0.0 -- main.canc.lua:184
local cost_concat_4 = 0.0 -- main.canc.lua:185
local cost_new = 0.0 -- main.canc.lua:186
for i = 2, dimen - 1 do -- main.canc.lua:188
local i_prev = i - 1 -- main.canc.lua:189
local i_next = i + 1 -- main.canc.lua:190
cost_concat_1 = solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (i_prev - 1)) + info["T"])] + info["c"][solut["s"][i_prev]][solut["s"][i_next]] -- main.canc.lua:192
cost_concat_2 = cost_concat_1 + solut["seq"][(3 * ((dimen + 1) * (i - 1) + (i_next - 1)) + info["T"])] + info["c"][solut["s"][i]][solut["s"][i_next + 1]] -- main.canc.lua:193
cost_new = solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (i_prev - 1)) + info["C"])] + solut["seq"][(3 * ((dimen + 1) * (i - 1) + (i_next - 1)) + info["W"])] * (cost_concat_1) + info["c"][solut["s"][i_next]][solut["s"][i]] + solut["seq"][(3 * ((dimen + 1) * (i_next + 1 - 1) + (dimen + 1 - 1)) + info["W"])] * (cost_concat_2) + solut["seq"][(3 * ((dimen + 1) * (i_next + 1 - 1) + (dimen + 1 - 1)) + info["C"])] -- main.canc.lua:197
if cost_new < cost_best then -- main.canc.lua:199
cost_best = cost_new - info["EPSILON"] -- main.canc.lua:200
I = i -- main.canc.lua:201
J = i_next -- main.canc.lua:202
end -- main.canc.lua:202
for j = i_next + 1, dimen do -- main.canc.lua:205
local j_prev = j - 1 -- main.canc.lua:206
local j_next = j + 1 -- main.canc.lua:207
cost_concat_1 = solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (i_prev - 1)) + info["T"])] + info["c"][solut["s"][i_prev]][solut["s"][j]] -- main.canc.lua:210
cost_concat_2 = cost_concat_1 + info["c"][solut["s"][j]][solut["s"][i_next]] -- main.canc.lua:211
cost_concat_3 = cost_concat_2 + solut["seq"][(3 * ((dimen + 1) * (i_next - 1) + (j_prev - 1)) + info["T"])] + info["c"][solut["s"][j_prev]][solut["s"][i]] -- main.canc.lua:212
cost_concat_4 = cost_concat_3 + info["c"][solut["s"][i]][solut["s"][j_next]] -- main.canc.lua:213
cost_new = solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (i_prev - 1)) + info["C"])] + cost_concat_1 + solut["seq"][(3 * ((dimen + 1) * (i_next - 1) + (j_prev - 1)) + info["W"])] * cost_concat_2 + solut["seq"][(3 * ((dimen + 1) * (i_next - 1) + (j_prev - 1)) + info["C"])] + cost_concat_3 + solut["seq"][(3 * ((dimen + 1) * (j_next - 1) + (dimen + 1 - 1)) + info["W"])] * cost_concat_4 + solut["seq"][(3 * ((dimen + 1) * (j_next - 1) + (dimen + 1 - 1)) + info["C"])] -- main.canc.lua:220
if (cost_new < cost_best) then -- main.canc.lua:222
cost_best = cost_new - info["EPSILON"] -- main.canc.lua:223
I = i -- main.canc.lua:224
J = j -- main.canc.lua:225
end -- main.canc.lua:225
end -- main.canc.lua:225
end -- main.canc.lua:225
if cost_best < solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (dimen + 1 - 1)) + info["C"])] - info["EPSILON"] then -- main.canc.lua:231
swap(solut["s"], I, J) -- main.canc.lua:234
subseq_load(solut, info) -- main.canc.lua:235
return true -- main.canc.lua:237
end -- main.canc.lua:237
return false -- main.canc.lua:240
end -- main.canc.lua:240
search_two_opt = function(solut, info) -- main.canc.lua:243
local cost_best = math["huge"] -- main.canc.lua:244
local I = - 1 -- main.canc.lua:245
local J = - 1 -- main.canc.lua:246
local dimen = info["dimension"] -- main.canc.lua:247
local cost_concat_1 = 0.0 -- main.canc.lua:249
local cost_concat_2 = 0.0 -- main.canc.lua:250
local cost_new = 0.0 -- main.canc.lua:251
for i = 2, dimen - 1 do -- main.canc.lua:253
local i_prev = i - 1 -- main.canc.lua:254
local rev_seq_cost = solut["seq"][(3 * ((dimen + 1) * (i - 1) + (i + 1 - 1)) + info["T"])] -- main.canc.lua:255
for j = i + 2, dimen do -- main.canc.lua:256
local j_next = j + 1 -- main.canc.lua:257
rev_seq_cost = rev_seq_cost + info["c"][solut["s"][j - 1]][solut["s"][j]] * (solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["W"])] - 1.0) -- main.canc.lua:259
cost_concat_1 = solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (i_prev - 1)) + info["T"])] + info["c"][solut["s"][j]][solut["s"][i_prev]] -- main.canc.lua:261
cost_concat_2 = cost_concat_1 + solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["T"])] + info["c"][solut["s"][j_next]][solut["s"][i]] -- main.canc.lua:262
cost_new = solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (i_prev - 1)) + info["C"])] + solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["W"])] * cost_concat_1 + rev_seq_cost + solut["seq"][(3 * ((dimen + 1) * (j_next - 1) + (dimen + 1 - 1)) + info["W"])] * cost_concat_2 + solut["seq"][(3 * ((dimen + 1) * (j_next - 1) + (dimen + 1 - 1)) + info["C"])] -- main.canc.lua:266
if cost_new < cost_best then -- main.canc.lua:268
cost_best = cost_new - info["EPSILON"] -- main.canc.lua:269
I = i -- main.canc.lua:270
J = j -- main.canc.lua:271
end -- main.canc.lua:271
end -- main.canc.lua:271
end -- main.canc.lua:271
if cost_best < solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (dimen + 1 - 1)) + info["C"])] - info["EPSILON"] then -- main.canc.lua:276
reverse(solut["s"], I, J) -- main.canc.lua:277
subseq_load(solut, info) -- main.canc.lua:278
return true -- main.canc.lua:279
end -- main.canc.lua:279
return false -- main.canc.lua:282
end -- main.canc.lua:282
search_reinsertion = function(solut, info, opt) -- main.canc.lua:285
local cost_best = math["huge"] -- main.canc.lua:286
local I = - 1 -- main.canc.lua:287
local J = - 1 -- main.canc.lua:288
local POS = - 1 -- main.canc.lua:289
local dimen = info["dimension"] -- main.canc.lua:290
local cost_concat_1 = 0.0 -- main.canc.lua:292
local cost_concat_2 = 0.0 -- main.canc.lua:293
local cost_concat_3 = 0.0 -- main.canc.lua:294
local cost_new = 0.0 -- main.canc.lua:295
for i = 2, dimen - opt + 1 do -- main.canc.lua:297
local j = opt + i - 1 -- main.canc.lua:298
local i_prev = i - 1 -- main.canc.lua:299
local j_next = j + 1 -- main.canc.lua:300
for k = 1, i_prev - 1 do -- main.canc.lua:302
local k_next = k + 1 -- main.canc.lua:303
cost_concat_1 = solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (k - 1)) + info["T"])] + info["c"][solut["s"][k]][solut["s"][i]] -- main.canc.lua:305
cost_concat_2 = cost_concat_1 + solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["T"])] + info["c"][solut["s"][j]][solut["s"][k_next]] -- main.canc.lua:306
cost_concat_3 = cost_concat_2 + solut["seq"][(3 * ((dimen + 1) * (k_next - 1) + (i_prev - 1)) + info["T"])] + info["c"][solut["s"][i_prev]][solut["s"][j_next]] -- main.canc.lua:307
cost_new = solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (k - 1)) + info["C"])] + solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["W"])] * cost_concat_1 + solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["C"])] + solut["seq"][(3 * ((dimen + 1) * (k_next - 1) + (i_prev - 1)) + info["W"])] * cost_concat_2 + solut["seq"][(3 * ((dimen + 1) * (k_next - 1) + (i_prev - 1)) + info["C"])] + solut["seq"][(3 * ((dimen + 1) * (j_next - 1) + (dimen + 1 - 1)) + info["W"])] * cost_concat_3 + solut["seq"][(3 * ((dimen + 1) * (j_next - 1) + (dimen + 1 - 1)) + info["C"])] -- main.canc.lua:312
if cost_new < cost_best then -- main.canc.lua:314
cost_best = cost_new - info["EPSILON"] -- main.canc.lua:315
I = i -- main.canc.lua:316
J = j -- main.canc.lua:317
POS = k -- main.canc.lua:318
end -- main.canc.lua:318
end -- main.canc.lua:318
for k = i + opt, dimen do -- main.canc.lua:323
local k_next = k + 1 -- main.canc.lua:324
cost_concat_1 = solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (i_prev - 1)) + info["T"])] + info["c"][solut["s"][i_prev]][solut["s"][j_next]] -- main.canc.lua:326
cost_concat_2 = cost_concat_1 + solut["seq"][(3 * ((dimen + 1) * (j_next - 1) + (k - 1)) + info["T"])] + info["c"][solut["s"][k]][solut["s"][i]] -- main.canc.lua:327
cost_concat_3 = cost_concat_2 + solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["T"])] + info["c"][solut["s"][j]][solut["s"][k_next]] -- main.canc.lua:328
cost_new = solut["seq"][(3 * ((dimen + 1) * (1 - 1) + (i_prev - 1)) + info["C"])] + solut["seq"][(3 * ((dimen + 1) * (j_next - 1) + (k - 1)) + info["W"])] * cost_concat_1 + solut["seq"][(3 * ((dimen + 1) * (j_next - 1) + (k - 1)) + info["C"])] + solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["W"])] * cost_concat_2 + solut["seq"][(3 * ((dimen + 1) * (i - 1) + (j - 1)) + info["C"])] + solut["seq"][(3 * ((dimen + 1) * (k_next - 1) + (dimen + 1 - 1)) + info["W"])] * cost_concat_3 + solut["seq"][(3 * ((dimen + 1) * (k_next - 1) + (dimen + 1 - 1)) + info["C"])] -- main.canc.lua:333
if cost_new < cost_best then -- main.canc.lua:335
cost_best = cost_new - info["EPSILON"] -- main.canc.lua:336
I = i -- main.canc.lua:337
J = j -- main.canc.lua:338
POS = k -- main.canc.lua:339
end -- main.canc.lua:339
end -- main.canc.lua:339
end -- main.canc.lua:339
if cost_best < solut["cost"] then -- main.canc.lua:346
reinsert(solut["s"], I, J, POS + 1) -- main.canc.lua:350
subseq_load(solut, info) -- main.canc.lua:352
if cost_best ~= solut["cost"] then -- main.canc.lua:355
print("ERROR") -- main.canc.lua:356
os["exit"](1) -- main.canc.lua:357
end -- main.canc.lua:357
return true -- main.canc.lua:359
end -- main.canc.lua:359
return false -- main.canc.lua:362
end -- main.canc.lua:362
RVND = function(solut, info) -- main.canc.lua:365
local SWAP = 0 -- main.canc.lua:366
local REINSERTION = 1 -- main.canc.lua:367
local OR_OPT_2 = 2 -- main.canc.lua:368
local OR_OPT_3 = 3 -- main.canc.lua:369
local TWO_OPT = 4 -- main.canc.lua:370
local neighbd_list = { -- main.canc.lua:372
SWAP, -- main.canc.lua:372
TWO_OPT, -- main.canc.lua:372
REINSERTION, -- main.canc.lua:372
OR_OPT_2, -- main.canc.lua:372
OR_OPT_3 -- main.canc.lua:372
} -- main.canc.lua:372
while # neighbd_list > 0 do -- main.canc.lua:380
local index = math["random"](1, # neighbd_list) -- main.canc.lua:381
index = info["rnd"][info["rnd_index"]] + 1 -- main.canc.lua:383
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:384
local neighbd = neighbd_list[index] -- main.canc.lua:386
local improve = false -- main.canc.lua:388
if neighbd == SWAP then -- main.canc.lua:390
improve = search_swap(solut, info) -- main.canc.lua:391
elseif neighbd == REINSERTION then -- main.canc.lua:393
improve = search_reinsertion(solut, info, REINSERTION) -- main.canc.lua:394
elseif neighbd == OR_OPT_2 then -- main.canc.lua:396
improve = search_reinsertion(solut, info, OR_OPT_2) -- main.canc.lua:397
elseif neighbd == OR_OPT_3 then -- main.canc.lua:399
improve = search_reinsertion(solut, info, OR_OPT_3) -- main.canc.lua:400
elseif neighbd == TWO_OPT then -- main.canc.lua:402
improve = search_two_opt(solut, info) -- main.canc.lua:403
end -- main.canc.lua:403
if improve == true then -- main.canc.lua:408
neighbd_list = { -- main.canc.lua:409
SWAP, -- main.canc.lua:409
TWO_OPT, -- main.canc.lua:409
REINSERTION, -- main.canc.lua:409
OR_OPT_2, -- main.canc.lua:409
OR_OPT_3 -- main.canc.lua:409
} -- main.canc.lua:409
else -- main.canc.lua:409
table["remove"](neighbd_list, index) -- main.canc.lua:412
end -- main.canc.lua:412
end -- main.canc.lua:412
end -- main.canc.lua:412
perturb = function(solut, info) -- main.canc.lua:419
local s = table["clone"](solut["s"]) -- main.canc.lua:420
local A_start = 1 -- main.canc.lua:422
local A_end = 1 -- main.canc.lua:423
local B_start = 1 -- main.canc.lua:424
local B_end = 1 -- main.canc.lua:425
local size_max = math["floor"](# s / 10) -- main.canc.lua:427
size_max = size_max >= 2 and size_max or 2 -- main.canc.lua:428
local size_min = 2 -- main.canc.lua:429
while (A_start <= B_start and B_start <= A_end) or (B_start <= A_start and A_start <= B_end) do -- main.canc.lua:431
A_start = math["random"](2, # s - 1 - size_max) -- main.canc.lua:432
A_end = A_start + math["random"](size_min, size_max) -- main.canc.lua:433
B_start = math["random"](2, # s - 1 - size_max) -- main.canc.lua:435
B_end = B_start + math["random"](size_min, size_max) -- main.canc.lua:436
A_start = info["rnd"][info["rnd_index"]] + 1 -- main.canc.lua:439
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:440
A_end = A_start + info["rnd"][info["rnd_index"]] -- main.canc.lua:441
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:442
B_start = info["rnd"][info["rnd_index"]] + 1 -- main.canc.lua:444
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:445
B_end = B_start + info["rnd"][info["rnd_index"]] -- main.canc.lua:446
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:447
end -- main.canc.lua:447
if A_start < B_start then -- main.canc.lua:450
reinsert(s, B_start, B_end - 1, A_end) -- main.canc.lua:451
reinsert(s, A_start, A_end - 1, B_end) -- main.canc.lua:452
else -- main.canc.lua:452
reinsert(s, A_start, A_end - 1, B_end) -- main.canc.lua:454
reinsert(s, B_start, B_end - 1, A_end) -- main.canc.lua:455
end -- main.canc.lua:455
return s -- main.canc.lua:458
end -- main.canc.lua:458
GILS_RVND = function(Imax, Iils, R, info) -- main.canc.lua:461
local solut_partial = { -- main.canc.lua:463
["s"] = {}, -- main.canc.lua:464
["seq"] = {} -- main.canc.lua:465
} -- main.canc.lua:465
subseq_fill(solut_partial["seq"], info) -- main.canc.lua:472
solut_partial["cost"] = 0 -- main.canc.lua:473
local solut_crnt = solut_clone(solut_partial) -- main.canc.lua:475
local solut_best = solut_clone(solut_crnt) -- main.canc.lua:476
solut_best["cost"] = math["huge"] -- main.canc.lua:478
for i = 1, Imax do -- main.canc.lua:480
local Rsz = # R -- main.canc.lua:481
local alpha = R[math["random"](1, Rsz)] -- main.canc.lua:482
alpha = R[info["rnd"][info["rnd_index"]] + 1] -- main.canc.lua:483
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:485
print("[+] Local Search", i) -- main.canc.lua:488
print("\9[+] Constructing Inital Solution..") -- main.canc.lua:489
solut_crnt["s"] = construction(alpha, info) -- main.canc.lua:490
subseq_load(solut_crnt, info) -- main.canc.lua:491
s_print(solut_crnt) -- main.canc.lua:492
print("\9Construction cost  ", solut_crnt["cost"]) -- main.canc.lua:493
solut_partial = solut_clone(solut_crnt) -- main.canc.lua:499
print("\9[+] Looking for the best Neighbor..") -- main.canc.lua:505
local iterILS = 0 -- main.canc.lua:506
while iterILS < Iils do -- main.canc.lua:507
RVND(solut_crnt, info) -- main.canc.lua:508
if solut_crnt["cost"] < solut_partial["cost"] - info["EPSILON"] then -- main.canc.lua:513
solut_partial["cost"] = solut_crnt["cost"] - info["EPSILON"] -- main.canc.lua:514
solut_partial["s"] = table["clone"](solut_crnt["s"]) -- main.canc.lua:515
iterILS = 0 -- main.canc.lua:516
end -- main.canc.lua:516
solut_crnt["s"] = perturb(solut_partial, info) -- main.canc.lua:519
subseq_load(solut_crnt, info) -- main.canc.lua:520
iterILS = iterILS + 1 -- main.canc.lua:522
end -- main.canc.lua:522
subseq_load(solut_partial, info) -- main.canc.lua:529
if solut_partial["cost"] < solut_best["cost"] then -- main.canc.lua:534
solut_best = solut_clone(solut_partial) -- main.canc.lua:535
end -- main.canc.lua:535
print("\9Current best solution cost: ", solut_best["cost"]) -- main.canc.lua:537
end -- main.canc.lua:537
print("COST: ", solut_best["cost"]) -- main.canc.lua:540
end -- main.canc.lua:540
protect = function(tbl) -- main.canc.lua:543
return setmetatable({}, { -- main.canc.lua:544
["__index"] = tbl, -- main.canc.lua:545
["__newindex"] = function(t, key, value) -- main.canc.lua:546
error("attempting to change constant " .. tostring(key) .. " to " .. tostring(value), 2) -- main.canc.lua:548
end -- main.canc.lua:548
}) -- main.canc.lua:548
end -- main.canc.lua:548
main = function() -- main.canc.lua:553
local info = { -- main.canc.lua:554
["c"] = {}, -- main.canc.lua:555
["T"] = 1, -- main.canc.lua:556
["C"] = 2, -- main.canc.lua:557
["W"] = 3, -- main.canc.lua:557
["EPSILON"] = 1e-15, -- main.canc.lua:558
["rnd"] = {}, -- main.canc.lua:559
["rnd_index"] = 1 -- main.canc.lua:560
} -- main.canc.lua:560
local a = 0 -- main.canc.lua:564
info["dimension"], a = readData(info["c"], info["rnd"]) -- main.canc.lua:565
print(info["rnd"][a]) -- main.canc.lua:566
math["randomseed"](os["time"]()) -- main.canc.lua:567
local Imax = 10 -- main.canc.lua:569
local Iils = math["min"](100, info["dimension"]) -- main.canc.lua:570
local R = { -- main.canc.lua:571
0.00, -- main.canc.lua:571
0.01, -- main.canc.lua:571
0.02, -- main.canc.lua:571
0.03, -- main.canc.lua:571
0.04, -- main.canc.lua:571
0.05, -- main.canc.lua:571
0.06, -- main.canc.lua:571
0.07, -- main.canc.lua:571
0.08, -- main.canc.lua:571
0.09, -- main.canc.lua:571
0.10, -- main.canc.lua:571
0.11, -- main.canc.lua:571
0.12, -- main.canc.lua:571
0.13, -- main.canc.lua:572
0.14, -- main.canc.lua:572
0.15, -- main.canc.lua:572
0.16, -- main.canc.lua:572
0.17, -- main.canc.lua:572
0.18, -- main.canc.lua:572
0.19, -- main.canc.lua:572
0.20, -- main.canc.lua:572
0.21, -- main.canc.lua:572
0.22, -- main.canc.lua:572
0.23, -- main.canc.lua:572
0.24, -- main.canc.lua:572
0.25 -- main.canc.lua:572
} -- main.canc.lua:572
local start = os["clock"]() -- main.canc.lua:576
GILS_RVND(Imax, Iils, R, info) -- main.canc.lua:577
print("TIME: ", os["clock"]() - start) -- main.canc.lua:579
end -- main.canc.lua:579
local profiler = require("profiler") -- main.canc.lua:582
profiler["start"]() -- main.canc.lua:583
main() -- main.canc.lua:585
profiler["stop"]() -- main.canc.lua:586
profiler["report"]("profiler.log") -- main.canc.lua:587
