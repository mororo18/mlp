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
solut_clone = function(src, dest) -- main.canc.lua:58
for i = 1, # src["s"] do -- main.canc.lua:62
dest["s"][i] = src["s"][i] -- main.canc.lua:63
end -- main.canc.lua:63
dest["cost"] = src["cost"] -- main.canc.lua:66
return dest -- main.canc.lua:68
end -- main.canc.lua:68
subseq_fill = function(seq, info) -- main.canc.lua:73
for i = 1, info["dimension"] + 1 do -- main.canc.lua:74
for j = i, info["dimension"] + 1 do -- main.canc.lua:75
for k = 1, 3 do -- main.canc.lua:76
seq[(3 * ((info["dimension"] + 1) * (i - 1) + (j - 1)) + k)] = 0.0 -- main.canc.lua:77
end -- main.canc.lua:77
end -- main.canc.lua:77
end -- main.canc.lua:77
end -- main.canc.lua:77
subseq_load = function(solut, info) -- main.canc.lua:83
local dimen = info["dimension"] -- main.canc.lua:84
local seq = solut["seq"] -- main.canc.lua:85
local T = info["T"] -- main.canc.lua:87
local W = info["W"] -- main.canc.lua:88
local C = info["C"] -- main.canc.lua:89
local c = info["c"] -- main.canc.lua:91
local s = solut["s"] -- main.canc.lua:92
local s_size = dimen + 1 -- main.canc.lua:93
for i = 1, s_size do -- main.canc.lua:97
local k = 1 - i - (i ~= 1 and 0 or 1) -- main.canc.lua:98
seq[(3 * ((s_size) * (i - 1) + (i - 1)) + T)] = 0.0 -- main.canc.lua:100
seq[(3 * ((s_size) * (i - 1) + (i - 1)) + C)] = 0.0 -- main.canc.lua:101
seq[(3 * ((s_size) * (i - 1) + (i - 1)) + W)] = (i ~= 1 and 1 or 0) -- main.canc.lua:102
local some_T = 0.0 -- main.canc.lua:104
local some_C = 0.0 -- main.canc.lua:105
local s_j_prev = s[i] -- main.canc.lua:107
for j = i + 1, s_size do -- main.canc.lua:109
local j_prev = j - 1 -- main.canc.lua:110
local s_j = s[j] -- main.canc.lua:112
local cost_T = c[s_j_prev][s_j] + some_T -- main.canc.lua:114
seq[(3 * ((s_size) * (i - 1) + (j - 1)) + T)] = cost_T -- main.canc.lua:115
local cost_C = cost_T + some_C -- main.canc.lua:117
seq[(3 * ((s_size) * (i - 1) + (j - 1)) + C)] = cost_C -- main.canc.lua:118
seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)] = j + k -- main.canc.lua:120
s_j_prev = s_j -- main.canc.lua:123
some_T = cost_T -- main.canc.lua:124
some_C = cost_C -- main.canc.lua:125
end -- main.canc.lua:125
end -- main.canc.lua:125
solut["cost"] = seq[(3 * ((dimen + 1) * (1 - 1) + (dimen + 1 - 1)) + C)] - info["EPSILON"] -- main.canc.lua:129
end -- main.canc.lua:129
sort = function(arr, r, info) -- main.canc.lua:133
for i = 1, # arr do -- main.canc.lua:134
for j = 1, # arr - i do -- main.canc.lua:135
if info["c"][r][arr[j]] > info["c"][r][arr[j + 1]] then -- main.canc.lua:136
local tmp = arr[j] -- main.canc.lua:137
arr[j] = arr[j + 1] -- main.canc.lua:138
arr[j + 1] = tmp -- main.canc.lua:139
end -- main.canc.lua:139
end -- main.canc.lua:139
end -- main.canc.lua:139
end -- main.canc.lua:139
construction = function(alpha, info) -- main.canc.lua:145
local s = { 1 } -- main.canc.lua:146
local cList = {} -- main.canc.lua:148
for i = 2, info["dimension"] do -- main.canc.lua:149
cList[# cList + 1] = i -- main.canc.lua:150
end -- main.canc.lua:150
local r = 1 -- main.canc.lua:155
while # cList > 0 do -- main.canc.lua:158
sort(cList, r, info) -- main.canc.lua:161
local range = math["floor"](# cList * alpha) + 1 -- main.canc.lua:164
local i = math["random"](1, range) -- main.canc.lua:165
i = info["rnd"][info["rnd_index"]] + 1 -- main.canc.lua:166
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:167
local c = table["remove"](cList, i) -- main.canc.lua:169
table["insert"](s, c) -- main.canc.lua:170
r = c -- main.canc.lua:172
end -- main.canc.lua:172
table["insert"](s, 1) -- main.canc.lua:174
return table["clone"](s) -- main.canc.lua:176
end -- main.canc.lua:176
swap = function(s, i, j) -- main.canc.lua:179
s[j], s[i] = s[i], s[j] -- main.canc.lua:180
end -- main.canc.lua:180
reverse = function(s, i, j) -- main.canc.lua:183
local l = j -- main.canc.lua:184
for k = i, math["floor"]((j + i) / 2) do -- main.canc.lua:185
swap(s, k, l) -- main.canc.lua:187
l = l - 1 -- main.canc.lua:188
end -- main.canc.lua:188
end -- main.canc.lua:188
reinsert = function(s, i, j, pos) -- main.canc.lua:192
if i < pos then -- main.canc.lua:193
for k = i, j do -- main.canc.lua:194
local e = s[i] -- main.canc.lua:195
table["insert"](s, pos, e) -- main.canc.lua:196
table["remove"](s, i) -- main.canc.lua:197
end -- main.canc.lua:197
else -- main.canc.lua:197
for k = i, j do -- main.canc.lua:200
local e = s[j] -- main.canc.lua:201
table["remove"](s, j) -- main.canc.lua:202
table["insert"](s, pos, e) -- main.canc.lua:203
end -- main.canc.lua:203
end -- main.canc.lua:203
end -- main.canc.lua:203
search_swap = function(solut, info) -- main.canc.lua:208
local cost_best = math["huge"] -- main.canc.lua:209
local I = - 1 -- main.canc.lua:210
local J = - 1 -- main.canc.lua:211
local dimen = info["dimension"] -- main.canc.lua:212
local cost_concat_1 = 0.0 -- main.canc.lua:214
local cost_concat_2 = 0.0 -- main.canc.lua:215
local cost_concat_3 = 0.0 -- main.canc.lua:216
local cost_concat_4 = 0.0 -- main.canc.lua:217
local cost_new = 0.0 -- main.canc.lua:218
local s = solut["s"] -- main.canc.lua:220
local seq = solut["seq"] -- main.canc.lua:221
local c = info["c"] -- main.canc.lua:222
local s_size = dimen + 1 -- main.canc.lua:224
local T = info["T"] -- main.canc.lua:226
local C = info["C"] -- main.canc.lua:227
local W = info["W"] -- main.canc.lua:228
local EP = info["EPSILON"] -- main.canc.lua:230
for i = 2, dimen - 1 do -- main.canc.lua:232
local i_prev = i - 1 -- main.canc.lua:233
local i_next = i + 1 -- main.canc.lua:234
cost_concat_1 = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + T)] + c[s[i_prev]][s[i_next]] -- main.canc.lua:236
cost_concat_2 = cost_concat_1 + seq[(3 * ((s_size) * (i - 1) + (i_next - 1)) + T)] + c[s[i]][s[i_next + 1]] -- main.canc.lua:237
cost_new = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + C)] + seq[(3 * ((s_size) * (i - 1) + (i_next - 1)) + W)] * (cost_concat_1) + c[s[i_next]][s[i]] + seq[(3 * ((s_size) * (i_next + 1 - 1) + (s_size - 1)) + W)] * (cost_concat_2) + seq[(3 * ((s_size) * (i_next + 1 - 1) + (s_size - 1)) + C)] -- main.canc.lua:241
if cost_new < cost_best then -- main.canc.lua:243
cost_best = cost_new - EP -- main.canc.lua:244
I = i -- main.canc.lua:245
J = i_next -- main.canc.lua:246
end -- main.canc.lua:246
for j = i_next + 1, dimen do -- main.canc.lua:249
local j_prev = j - 1 -- main.canc.lua:250
local j_next = j + 1 -- main.canc.lua:251
cost_concat_1 = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + T)] + c[s[i_prev]][s[j]] -- main.canc.lua:254
cost_concat_2 = cost_concat_1 + c[s[j]][s[i_next]] -- main.canc.lua:255
cost_concat_3 = cost_concat_2 + seq[(3 * ((s_size) * (i_next - 1) + (j_prev - 1)) + T)] + c[s[j_prev]][s[i]] -- main.canc.lua:256
cost_concat_4 = cost_concat_3 + c[s[i]][s[j_next]] -- main.canc.lua:257
cost_new = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + C)] + cost_concat_1 + seq[(3 * ((s_size) * (i_next - 1) + (j_prev - 1)) + W)] * cost_concat_2 + seq[(3 * ((s_size) * (i_next - 1) + (j_prev - 1)) + C)] + cost_concat_3 + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + W)] * cost_concat_4 + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + C)] -- main.canc.lua:264
if (cost_new < cost_best) then -- main.canc.lua:266
cost_best = cost_new - EP -- main.canc.lua:267
I = i -- main.canc.lua:268
J = j -- main.canc.lua:269
end -- main.canc.lua:269
end -- main.canc.lua:269
end -- main.canc.lua:269
if cost_best < solut["seq"][(3 * ((s_size) * (1 - 1) + (s_size - 1)) + info["C"])] - info["EPSILON"] then -- main.canc.lua:275
swap(solut["s"], I, J) -- main.canc.lua:278
subseq_load(solut, info) -- main.canc.lua:279
return true -- main.canc.lua:281
end -- main.canc.lua:281
return false -- main.canc.lua:284
end -- main.canc.lua:284
search_two_opt = function(solut, info) -- main.canc.lua:287
local cost_best = math["huge"] -- main.canc.lua:288
local I = - 1 -- main.canc.lua:289
local J = - 1 -- main.canc.lua:290
local dimen = info["dimension"] -- main.canc.lua:291
local cost_concat_1 = 0.0 -- main.canc.lua:293
local cost_concat_2 = 0.0 -- main.canc.lua:294
local cost_new = 0.0 -- main.canc.lua:295
local s = solut["s"] -- main.canc.lua:297
local seq = solut["seq"] -- main.canc.lua:298
local c = info["c"] -- main.canc.lua:299
local T = info["T"] -- main.canc.lua:301
local C = info["C"] -- main.canc.lua:302
local W = info["W"] -- main.canc.lua:303
local EP = info["EPSILON"] -- main.canc.lua:305
local s_size = dimen + 1 -- main.canc.lua:307
for i = 2, dimen - 1 do -- main.canc.lua:309
local i_prev = i - 1 -- main.canc.lua:310
local rev_seq_cost = seq[(3 * ((s_size) * (i - 1) + (i + 1 - 1)) + T)] -- main.canc.lua:311
for j = i + 2, dimen do -- main.canc.lua:312
local j_next = j + 1 -- main.canc.lua:313
rev_seq_cost = rev_seq_cost + c[s[j - 1]][s[j]] * (seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)] - 1.0) -- main.canc.lua:315
cost_concat_1 = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + T)] + c[s[j]][s[i_prev]] -- main.canc.lua:317
cost_concat_2 = cost_concat_1 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + T)] + c[s[j_next]][s[i]] -- main.canc.lua:318
cost_new = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + C)] + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)] * cost_concat_1 + rev_seq_cost + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + W)] * cost_concat_2 + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + C)] -- main.canc.lua:322
if cost_new < cost_best then -- main.canc.lua:324
cost_best = cost_new - EP -- main.canc.lua:325
I = i -- main.canc.lua:326
J = j -- main.canc.lua:327
end -- main.canc.lua:327
end -- main.canc.lua:327
end -- main.canc.lua:327
if cost_best < solut["seq"][(3 * ((s_size) * (1 - 1) + (s_size - 1)) + C)] - EP then -- main.canc.lua:332
reverse(solut["s"], I, J) -- main.canc.lua:333
subseq_load(solut, info) -- main.canc.lua:334
return true -- main.canc.lua:335
end -- main.canc.lua:335
return false -- main.canc.lua:338
end -- main.canc.lua:338
search_reinsertion = function(solut, info, opt) -- main.canc.lua:341
local cost_best = math["huge"] -- main.canc.lua:342
local I = - 1 -- main.canc.lua:343
local J = - 1 -- main.canc.lua:344
local POS = - 1 -- main.canc.lua:345
local dimen = info["dimension"] -- main.canc.lua:346
local cost_concat_1 = 0.0 -- main.canc.lua:348
local cost_concat_2 = 0.0 -- main.canc.lua:349
local cost_concat_3 = 0.0 -- main.canc.lua:350
local cost_new = 0.0 -- main.canc.lua:351
local s = solut["s"] -- main.canc.lua:353
local seq = solut["seq"] -- main.canc.lua:354
local c = info["c"] -- main.canc.lua:355
local s_size = dimen + 1 -- main.canc.lua:356
local T = info["T"] -- main.canc.lua:358
local C = info["C"] -- main.canc.lua:359
local W = info["W"] -- main.canc.lua:360
local EP = info["EPSILON"] -- main.canc.lua:362
for i = 2, dimen - opt + 1 do -- main.canc.lua:364
local j = opt + i - 1 -- main.canc.lua:365
local i_prev = i - 1 -- main.canc.lua:366
local j_next = j + 1 -- main.canc.lua:367
for k = 1, i_prev - 1 do -- main.canc.lua:369
local k_next = k + 1 -- main.canc.lua:370
cost_concat_1 = seq[(3 * ((s_size) * (1 - 1) + (k - 1)) + T)] + c[s[k]][s[i]] -- main.canc.lua:372
cost_concat_2 = cost_concat_1 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + T)] + c[s[j]][s[k_next]] -- main.canc.lua:373
cost_concat_3 = cost_concat_2 + seq[(3 * ((s_size) * (k_next - 1) + (i_prev - 1)) + T)] + c[s[i_prev]][s[j_next]] -- main.canc.lua:374
cost_new = seq[(3 * ((s_size) * (1 - 1) + (k - 1)) + C)] + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)] * cost_concat_1 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + C)] + seq[(3 * ((s_size) * (k_next - 1) + (i_prev - 1)) + W)] * cost_concat_2 + seq[(3 * ((s_size) * (k_next - 1) + (i_prev - 1)) + C)] + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + W)] * cost_concat_3 + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + C)] -- main.canc.lua:379
if cost_new < cost_best then -- main.canc.lua:381
cost_best = cost_new - EP -- main.canc.lua:382
I = i -- main.canc.lua:383
J = j -- main.canc.lua:384
POS = k -- main.canc.lua:385
end -- main.canc.lua:385
end -- main.canc.lua:385
for k = i + opt, dimen do -- main.canc.lua:390
local k_next = k + 1 -- main.canc.lua:391
cost_concat_1 = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + T)] + c[s[i_prev]][s[j_next]] -- main.canc.lua:393
cost_concat_2 = cost_concat_1 + seq[(3 * ((s_size) * (j_next - 1) + (k - 1)) + T)] + c[s[k]][s[i]] -- main.canc.lua:394
cost_concat_3 = cost_concat_2 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + T)] + c[s[j]][s[k_next]] -- main.canc.lua:395
cost_new = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + C)] + seq[(3 * ((s_size) * (j_next - 1) + (k - 1)) + W)] * cost_concat_1 + seq[(3 * ((s_size) * (j_next - 1) + (k - 1)) + C)] + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)] * cost_concat_2 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + C)] + seq[(3 * ((s_size) * (k_next - 1) + (s_size - 1)) + W)] * cost_concat_3 + seq[(3 * ((s_size) * (k_next - 1) + (s_size - 1)) + C)] -- main.canc.lua:400
if cost_new < cost_best then -- main.canc.lua:402
cost_best = cost_new - EP -- main.canc.lua:403
I = i -- main.canc.lua:404
J = j -- main.canc.lua:405
POS = k -- main.canc.lua:406
end -- main.canc.lua:406
end -- main.canc.lua:406
end -- main.canc.lua:406
if cost_best < solut["cost"] then -- main.canc.lua:413
reinsert(solut["s"], I, J, POS + 1) -- main.canc.lua:417
subseq_load(solut, info) -- main.canc.lua:419
if cost_best ~= solut["cost"] then -- main.canc.lua:422
print("ERROR") -- main.canc.lua:423
os["exit"](1) -- main.canc.lua:424
end -- main.canc.lua:424
return true -- main.canc.lua:426
end -- main.canc.lua:426
return false -- main.canc.lua:429
end -- main.canc.lua:429
RVND = function(solut, info) -- main.canc.lua:432
local SWAP = 0 -- main.canc.lua:433
local REINSERTION = 1 -- main.canc.lua:434
local OR_OPT_2 = 2 -- main.canc.lua:435
local OR_OPT_3 = 3 -- main.canc.lua:436
local TWO_OPT = 4 -- main.canc.lua:437
local neighbd_list = { -- main.canc.lua:439
SWAP, -- main.canc.lua:439
TWO_OPT, -- main.canc.lua:439
REINSERTION, -- main.canc.lua:439
OR_OPT_2, -- main.canc.lua:439
OR_OPT_3 -- main.canc.lua:439
} -- main.canc.lua:439
while # neighbd_list > 0 do -- main.canc.lua:447
local index = math["random"](1, # neighbd_list) -- main.canc.lua:448
index = info["rnd"][info["rnd_index"]] + 1 -- main.canc.lua:450
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:451
local neighbd = neighbd_list[index] -- main.canc.lua:453
local improve = false -- main.canc.lua:455
if neighbd == SWAP then -- main.canc.lua:457
improve = search_swap(solut, info) -- main.canc.lua:458
elseif neighbd == REINSERTION then -- main.canc.lua:460
improve = search_reinsertion(solut, info, REINSERTION) -- main.canc.lua:461
elseif neighbd == OR_OPT_2 then -- main.canc.lua:463
improve = search_reinsertion(solut, info, OR_OPT_2) -- main.canc.lua:464
elseif neighbd == OR_OPT_3 then -- main.canc.lua:466
improve = search_reinsertion(solut, info, OR_OPT_3) -- main.canc.lua:467
elseif neighbd == TWO_OPT then -- main.canc.lua:469
improve = search_two_opt(solut, info) -- main.canc.lua:470
end -- main.canc.lua:470
if improve == true then -- main.canc.lua:475
neighbd_list = { -- main.canc.lua:476
SWAP, -- main.canc.lua:476
TWO_OPT, -- main.canc.lua:476
REINSERTION, -- main.canc.lua:476
OR_OPT_2, -- main.canc.lua:476
OR_OPT_3 -- main.canc.lua:476
} -- main.canc.lua:476
else -- main.canc.lua:476
table["remove"](neighbd_list, index) -- main.canc.lua:479
end -- main.canc.lua:479
end -- main.canc.lua:479
end -- main.canc.lua:479
perturb = function(solut, info) -- main.canc.lua:486
local s = table["clone"](solut["s"]) -- main.canc.lua:487
local A_start = 1 -- main.canc.lua:489
local A_end = 1 -- main.canc.lua:490
local B_start = 1 -- main.canc.lua:491
local B_end = 1 -- main.canc.lua:492
local size_max = math["floor"](# s / 10) -- main.canc.lua:494
size_max = size_max >= 2 and size_max or 2 -- main.canc.lua:495
local size_min = 2 -- main.canc.lua:496
while (A_start <= B_start and B_start <= A_end) or (B_start <= A_start and A_start <= B_end) do -- main.canc.lua:498
A_start = math["random"](2, # s - 1 - size_max) -- main.canc.lua:499
A_end = A_start + math["random"](size_min, size_max) -- main.canc.lua:500
B_start = math["random"](2, # s - 1 - size_max) -- main.canc.lua:502
B_end = B_start + math["random"](size_min, size_max) -- main.canc.lua:503
A_start = info["rnd"][info["rnd_index"]] + 1 -- main.canc.lua:506
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:507
A_end = A_start + info["rnd"][info["rnd_index"]] -- main.canc.lua:508
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:509
B_start = info["rnd"][info["rnd_index"]] + 1 -- main.canc.lua:511
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:512
B_end = B_start + info["rnd"][info["rnd_index"]] -- main.canc.lua:513
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:514
end -- main.canc.lua:514
if A_start < B_start then -- main.canc.lua:517
reinsert(s, B_start, B_end - 1, A_end) -- main.canc.lua:518
reinsert(s, A_start, A_end - 1, B_end) -- main.canc.lua:519
else -- main.canc.lua:519
reinsert(s, A_start, A_end - 1, B_end) -- main.canc.lua:521
reinsert(s, B_start, B_end - 1, A_end) -- main.canc.lua:522
end -- main.canc.lua:522
return s -- main.canc.lua:525
end -- main.canc.lua:525
GILS_RVND = function(Imax, Iils, R, info) -- main.canc.lua:528
local solut_partial = { -- main.canc.lua:530
["s"] = {}, -- main.canc.lua:531
["seq"] = {} -- main.canc.lua:532
} -- main.canc.lua:532
local solut_crnt = { -- main.canc.lua:535
["s"] = {}, -- main.canc.lua:536
["seq"] = {} -- main.canc.lua:537
} -- main.canc.lua:537
local solut_best = { -- main.canc.lua:540
["s"] = {}, -- main.canc.lua:541
["seq"] = {} -- main.canc.lua:542
} -- main.canc.lua:542
subseq_fill(solut_partial["seq"], info) -- main.canc.lua:545
solut_partial["cost"] = 0 -- main.canc.lua:546
subseq_fill(solut_crnt["seq"], info) -- main.canc.lua:548
solut_crnt["cost"] = 0 -- main.canc.lua:549
subseq_fill(solut_best["seq"], info) -- main.canc.lua:551
solut_best["cost"] = 0 -- main.canc.lua:552
solut_best["cost"] = math["huge"] -- main.canc.lua:557
for i = 1, Imax do -- main.canc.lua:559
local Rsz = # R -- main.canc.lua:560
local alpha = R[math["random"](1, Rsz)] -- main.canc.lua:561
alpha = R[info["rnd"][info["rnd_index"]] + 1] -- main.canc.lua:562
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:564
print("[+] Local Search", i) -- main.canc.lua:567
print("\9[+] Constructing Inital Solution..") -- main.canc.lua:568
solut_crnt["s"] = construction(alpha, info) -- main.canc.lua:569
subseq_load(solut_crnt, info) -- main.canc.lua:570
s_print(solut_crnt) -- main.canc.lua:571
print("\9Construction cost  ", solut_crnt["cost"]) -- main.canc.lua:572
solut_clone(solut_crnt, solut_partial) -- main.canc.lua:578
print("\9[+] Looking for the best Neighbor..") -- main.canc.lua:585
local iterILS = 0 -- main.canc.lua:586
while iterILS < Iils do -- main.canc.lua:587
RVND(solut_crnt, info) -- main.canc.lua:588
if solut_crnt["cost"] < solut_partial["cost"] - info["EPSILON"] then -- main.canc.lua:593
solut_partial["cost"] = solut_crnt["cost"] - info["EPSILON"] -- main.canc.lua:594
solut_partial["s"] = table["clone"](solut_crnt["s"]) -- main.canc.lua:595
iterILS = 0 -- main.canc.lua:596
end -- main.canc.lua:596
solut_crnt["s"] = perturb(solut_partial, info) -- main.canc.lua:599
subseq_load(solut_crnt, info) -- main.canc.lua:600
iterILS = iterILS + 1 -- main.canc.lua:602
end -- main.canc.lua:602
subseq_load(solut_partial, info) -- main.canc.lua:609
if solut_partial["cost"] < solut_best["cost"] then -- main.canc.lua:614
solut_clone(solut_partial, solut_best) -- main.canc.lua:615
end -- main.canc.lua:615
print("\9Current best solution cost: ", solut_best["cost"]) -- main.canc.lua:618
end -- main.canc.lua:618
print("COST: ", solut_best["cost"]) -- main.canc.lua:621
end -- main.canc.lua:621
protect = function(tbl) -- main.canc.lua:624
return setmetatable({}, { -- main.canc.lua:625
["__index"] = tbl, -- main.canc.lua:626
["__newindex"] = function(t, key, value) -- main.canc.lua:627
error("attempting to change constant " .. tostring(key) .. " to " .. tostring(value), 2) -- main.canc.lua:629
end -- main.canc.lua:629
}) -- main.canc.lua:629
end -- main.canc.lua:629
main = function() -- main.canc.lua:634
local info = { -- main.canc.lua:635
["c"] = {}, -- main.canc.lua:636
["T"] = 1, -- main.canc.lua:637
["C"] = 2, -- main.canc.lua:638
["W"] = 3, -- main.canc.lua:638
["EPSILON"] = 1e-15, -- main.canc.lua:639
["rnd"] = {}, -- main.canc.lua:640
["rnd_index"] = 1 -- main.canc.lua:641
} -- main.canc.lua:641
local a = 0 -- main.canc.lua:645
info["dimension"], a = readData(info["c"], info["rnd"]) -- main.canc.lua:646
print(info["rnd"][a]) -- main.canc.lua:647
math["randomseed"](os["time"]()) -- main.canc.lua:648
local Imax = 10 -- main.canc.lua:650
local Iils = math["min"](100, info["dimension"]) -- main.canc.lua:651
local R = { -- main.canc.lua:652
0.00, -- main.canc.lua:652
0.01, -- main.canc.lua:652
0.02, -- main.canc.lua:652
0.03, -- main.canc.lua:652
0.04, -- main.canc.lua:652
0.05, -- main.canc.lua:652
0.06, -- main.canc.lua:652
0.07, -- main.canc.lua:652
0.08, -- main.canc.lua:652
0.09, -- main.canc.lua:652
0.10, -- main.canc.lua:652
0.11, -- main.canc.lua:652
0.12, -- main.canc.lua:652
0.13, -- main.canc.lua:653
0.14, -- main.canc.lua:653
0.15, -- main.canc.lua:653
0.16, -- main.canc.lua:653
0.17, -- main.canc.lua:653
0.18, -- main.canc.lua:653
0.19, -- main.canc.lua:653
0.20, -- main.canc.lua:653
0.21, -- main.canc.lua:653
0.22, -- main.canc.lua:653
0.23, -- main.canc.lua:653
0.24, -- main.canc.lua:653
0.25 -- main.canc.lua:653
} -- main.canc.lua:653
local start = os["clock"]() -- main.canc.lua:657
GILS_RVND(Imax, Iils, R, info) -- main.canc.lua:658
print("TIME: ", os["clock"]() - start) -- main.canc.lua:660
end -- main.canc.lua:660
local profiler = require("profiler") -- main.canc.lua:665
profiler["start"]() -- main.canc.lua:666
main() -- main.canc.lua:668
profiler["stop"]() -- main.canc.lua:669
profiler["report"]("profiler.log") -- main.canc.lua:670
