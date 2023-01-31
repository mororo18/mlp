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
local cost_T = c[s_j_prev][s_j] + some_T -- main.canc.lua:115
seq[(3 * ((s_size) * (i - 1) + (j - 1)) + T)] = cost_T -- main.canc.lua:116
local cost_C = cost_T + some_C -- main.canc.lua:119
seq[(3 * ((s_size) * (i - 1) + (j - 1)) + C)] = cost_C -- main.canc.lua:120
seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)] = j + k -- main.canc.lua:123
s_j_prev = s_j -- main.canc.lua:126
some_T = cost_T -- main.canc.lua:127
some_C = cost_C -- main.canc.lua:128
end -- main.canc.lua:128
end -- main.canc.lua:128
solut["cost"] = seq[(3 * ((dimen + 1) * (1 - 1) + (dimen + 1 - 1)) + C)] - info["EPSILON"] -- main.canc.lua:132
end -- main.canc.lua:132
sort = function(arr, r, info) -- main.canc.lua:136
for i = 1, # arr do -- main.canc.lua:137
for j = 1, # arr - i do -- main.canc.lua:138
if info["c"][r][arr[j]] > info["c"][r][arr[j + 1]] then -- main.canc.lua:139
local tmp = arr[j] -- main.canc.lua:140
arr[j] = arr[j + 1] -- main.canc.lua:141
arr[j + 1] = tmp -- main.canc.lua:142
end -- main.canc.lua:142
end -- main.canc.lua:142
end -- main.canc.lua:142
end -- main.canc.lua:142
construction = function(alpha, info) -- main.canc.lua:148
local s = { 1 } -- main.canc.lua:149
local cList = {} -- main.canc.lua:151
for i = 2, info["dimension"] do -- main.canc.lua:152
cList[# cList + 1] = i -- main.canc.lua:153
end -- main.canc.lua:153
local r = 1 -- main.canc.lua:158
while # cList > 0 do -- main.canc.lua:161
sort(cList, r, info) -- main.canc.lua:164
local range = math["floor"](# cList * alpha) + 1 -- main.canc.lua:167
local i = math["random"](1, range) -- main.canc.lua:168
i = info["rnd"][info["rnd_index"]] + 1 -- main.canc.lua:169
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:170
local c = table["remove"](cList, i) -- main.canc.lua:172
table["insert"](s, c) -- main.canc.lua:173
r = c -- main.canc.lua:175
end -- main.canc.lua:175
table["insert"](s, 1) -- main.canc.lua:177
return table["clone"](s) -- main.canc.lua:179
end -- main.canc.lua:179
swap = function(s, i, j) -- main.canc.lua:182
s[j], s[i] = s[i], s[j] -- main.canc.lua:183
end -- main.canc.lua:183
reverse = function(s, i, j) -- main.canc.lua:186
local l = j -- main.canc.lua:187
for k = i, math["floor"]((j + i) / 2) do -- main.canc.lua:188
swap(s, k, l) -- main.canc.lua:190
l = l - 1 -- main.canc.lua:191
end -- main.canc.lua:191
end -- main.canc.lua:191
reinsert = function(s, i, j, pos) -- main.canc.lua:195
if i < pos then -- main.canc.lua:196
for k = i, j do -- main.canc.lua:197
local e = s[i] -- main.canc.lua:198
table["insert"](s, pos, e) -- main.canc.lua:199
table["remove"](s, i) -- main.canc.lua:200
end -- main.canc.lua:200
else -- main.canc.lua:200
for k = i, j do -- main.canc.lua:203
local e = s[j] -- main.canc.lua:204
table["remove"](s, j) -- main.canc.lua:205
table["insert"](s, pos, e) -- main.canc.lua:206
end -- main.canc.lua:206
end -- main.canc.lua:206
end -- main.canc.lua:206
search_swap = function(solut, info) -- main.canc.lua:211
local cost_best = math["huge"] -- main.canc.lua:212
local I = - 1 -- main.canc.lua:213
local J = - 1 -- main.canc.lua:214
local dimen = info["dimension"] -- main.canc.lua:215
local cost_concat_1 = 0.0 -- main.canc.lua:217
local cost_concat_2 = 0.0 -- main.canc.lua:218
local cost_concat_3 = 0.0 -- main.canc.lua:219
local cost_concat_4 = 0.0 -- main.canc.lua:220
local cost_new = 0.0 -- main.canc.lua:221
local s = solut["s"] -- main.canc.lua:223
local seq = solut["seq"] -- main.canc.lua:224
local c = info["c"] -- main.canc.lua:225
local s_size = dimen + 1 -- main.canc.lua:227
local T = info["T"] -- main.canc.lua:229
local C = info["C"] -- main.canc.lua:230
local W = info["W"] -- main.canc.lua:231
local EP = info["EPSILON"] -- main.canc.lua:233
for i = 2, dimen - 1 do -- main.canc.lua:235
local i_prev = i - 1 -- main.canc.lua:236
local i_next = i + 1 -- main.canc.lua:237
local s_i = s[i] -- main.canc.lua:239
local s_i_next = s[i_next] -- main.canc.lua:240
cost_concat_1 = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + T)] + c[s[i_prev]][s_i_next] -- main.canc.lua:242
cost_concat_2 = cost_concat_1 + seq[(3 * ((s_size) * (i - 1) + (i_next - 1)) + T)] + c[s_i][s[i_next + 1]] -- main.canc.lua:243
cost_new = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + C)] + seq[(3 * ((s_size) * (i - 1) + (i_next - 1)) + W)] * (cost_concat_1) + c[s_i_next][s_i] + seq[(3 * ((s_size) * (i_next + 1 - 1) + (s_size - 1)) + W)] * (cost_concat_2) + seq[(3 * ((s_size) * (i_next + 1 - 1) + (s_size - 1)) + C)] -- main.canc.lua:247
if cost_new < cost_best then -- main.canc.lua:249
cost_best = cost_new - EP -- main.canc.lua:250
I = i -- main.canc.lua:251
J = i_next -- main.canc.lua:252
end -- main.canc.lua:252
for j = i_next + 1, dimen do -- main.canc.lua:255
local j_prev = j - 1 -- main.canc.lua:256
local j_next = j + 1 -- main.canc.lua:257
local s_j = s[j] -- main.canc.lua:258
cost_concat_1 = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + T)] + c[s[i_prev]][s_j] -- main.canc.lua:261
cost_concat_2 = cost_concat_1 + c[s_j][s_i_next] -- main.canc.lua:262
cost_concat_3 = cost_concat_2 + seq[(3 * ((s_size) * (i_next - 1) + (j_prev - 1)) + T)] + c[s[j_prev]][s_i] -- main.canc.lua:263
cost_concat_4 = cost_concat_3 + c[s_i][s[j_next]] -- main.canc.lua:264
cost_new = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + C)] + cost_concat_1 + seq[(3 * ((s_size) * (i_next - 1) + (j_prev - 1)) + W)] * cost_concat_2 + seq[(3 * ((s_size) * (i_next - 1) + (j_prev - 1)) + C)] + cost_concat_3 + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + W)] * cost_concat_4 + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + C)] -- main.canc.lua:271
if (cost_new < cost_best) then -- main.canc.lua:273
cost_best = cost_new - EP -- main.canc.lua:274
I = i -- main.canc.lua:275
J = j -- main.canc.lua:276
end -- main.canc.lua:276
end -- main.canc.lua:276
end -- main.canc.lua:276
if cost_best < solut["seq"][(3 * ((s_size) * (1 - 1) + (s_size - 1)) + info["C"])] - info["EPSILON"] then -- main.canc.lua:282
swap(solut["s"], I, J) -- main.canc.lua:285
subseq_load(solut, info) -- main.canc.lua:286
return true -- main.canc.lua:288
end -- main.canc.lua:288
return false -- main.canc.lua:291
end -- main.canc.lua:291
search_two_opt = function(solut, info) -- main.canc.lua:294
local cost_best = math["huge"] -- main.canc.lua:295
local I = - 1 -- main.canc.lua:296
local J = - 1 -- main.canc.lua:297
local dimen = info["dimension"] -- main.canc.lua:298
local cost_concat_1 = 0.0 -- main.canc.lua:300
local cost_concat_2 = 0.0 -- main.canc.lua:301
local cost_new = 0.0 -- main.canc.lua:302
local s = solut["s"] -- main.canc.lua:304
local seq = solut["seq"] -- main.canc.lua:305
local c = info["c"] -- main.canc.lua:306
local T = info["T"] -- main.canc.lua:308
local C = info["C"] -- main.canc.lua:309
local W = info["W"] -- main.canc.lua:310
local EP = info["EPSILON"] -- main.canc.lua:312
local s_size = dimen + 1 -- main.canc.lua:314
for i = 2, dimen - 1 do -- main.canc.lua:316
local i_prev = i - 1 -- main.canc.lua:317
local rev_seq_cost = seq[(3 * ((s_size) * (i - 1) + (i + 1 - 1)) + T)] -- main.canc.lua:318
local s_j_prev = s[i + 1] -- main.canc.lua:320
for j = i + 2, dimen do -- main.canc.lua:322
local j_next = j + 1 -- main.canc.lua:323
local s_j = s[j] -- main.canc.lua:324
rev_seq_cost = rev_seq_cost + c[s_j_prev][s_j] * (seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)] - 1.0) -- main.canc.lua:326
cost_concat_1 = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + T)] + c[s_j][s[i_prev]] -- main.canc.lua:328
cost_concat_2 = cost_concat_1 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + T)] + c[s[j_next]][s[i]] -- main.canc.lua:329
cost_new = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + C)] + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)] * cost_concat_1 + rev_seq_cost + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + W)] * cost_concat_2 + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + C)] -- main.canc.lua:333
if cost_new < cost_best then -- main.canc.lua:335
cost_best = cost_new - EP -- main.canc.lua:336
I = i -- main.canc.lua:337
J = j -- main.canc.lua:338
end -- main.canc.lua:338
s_j_prev = s_j -- main.canc.lua:341
end -- main.canc.lua:341
end -- main.canc.lua:341
if cost_best < solut["seq"][(3 * ((s_size) * (1 - 1) + (s_size - 1)) + C)] - EP then -- main.canc.lua:345
reverse(solut["s"], I, J) -- main.canc.lua:346
subseq_load(solut, info) -- main.canc.lua:347
return true -- main.canc.lua:348
end -- main.canc.lua:348
return false -- main.canc.lua:351
end -- main.canc.lua:351
search_reinsertion = function(solut, info, opt) -- main.canc.lua:354
local cost_best = math["huge"] -- main.canc.lua:355
local I = - 1 -- main.canc.lua:356
local J = - 1 -- main.canc.lua:357
local POS = - 1 -- main.canc.lua:358
local dimen = info["dimension"] -- main.canc.lua:359
local cost_concat_1 = 0.0 -- main.canc.lua:361
local cost_concat_2 = 0.0 -- main.canc.lua:362
local cost_concat_3 = 0.0 -- main.canc.lua:363
local cost_new = 0.0 -- main.canc.lua:364
local s = solut["s"] -- main.canc.lua:366
local seq = solut["seq"] -- main.canc.lua:367
local c = info["c"] -- main.canc.lua:368
local s_size = dimen + 1 -- main.canc.lua:369
local T = info["T"] -- main.canc.lua:371
local C = info["C"] -- main.canc.lua:372
local W = info["W"] -- main.canc.lua:373
local EP = info["EPSILON"] -- main.canc.lua:375
for i = 2, dimen - opt + 1 do -- main.canc.lua:377
local j = opt + i - 1 -- main.canc.lua:378
local i_prev = i - 1 -- main.canc.lua:379
local j_next = j + 1 -- main.canc.lua:380
local s_i = s[i] -- main.canc.lua:382
local s_j = s[j] -- main.canc.lua:383
local s_i_prev = s[i_prev] -- main.canc.lua:385
local s_j_next = s[j_next] -- main.canc.lua:386
for k = 1, i_prev - 1 do -- main.canc.lua:388
local k_next = k + 1 -- main.canc.lua:389
cost_concat_1 = seq[(3 * ((s_size) * (1 - 1) + (k - 1)) + T)] + c[s[k]][s_i] -- main.canc.lua:391
cost_concat_2 = cost_concat_1 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + T)] + c[s_j][s[k_next]] -- main.canc.lua:392
cost_concat_3 = cost_concat_2 + seq[(3 * ((s_size) * (k_next - 1) + (i_prev - 1)) + T)] + c[s_i_prev][s_j_next] -- main.canc.lua:393
cost_new = seq[(3 * ((s_size) * (1 - 1) + (k - 1)) + C)] + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)] * cost_concat_1 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + C)] + seq[(3 * ((s_size) * (k_next - 1) + (i_prev - 1)) + W)] * cost_concat_2 + seq[(3 * ((s_size) * (k_next - 1) + (i_prev - 1)) + C)] + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + W)] * cost_concat_3 + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + C)] -- main.canc.lua:398
if cost_new < cost_best then -- main.canc.lua:400
cost_best = cost_new - EP -- main.canc.lua:401
I = i -- main.canc.lua:402
J = j -- main.canc.lua:403
POS = k -- main.canc.lua:404
end -- main.canc.lua:404
end -- main.canc.lua:404
for k = i + opt, dimen do -- main.canc.lua:409
local k_next = k + 1 -- main.canc.lua:410
cost_concat_1 = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + T)] + c[s_i_prev][s_j_next] -- main.canc.lua:412
cost_concat_2 = cost_concat_1 + seq[(3 * ((s_size) * (j_next - 1) + (k - 1)) + T)] + c[s[k]][s_i] -- main.canc.lua:413
cost_concat_3 = cost_concat_2 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + T)] + c[s_j][s[k_next]] -- main.canc.lua:414
cost_new = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + C)] + seq[(3 * ((s_size) * (j_next - 1) + (k - 1)) + W)] * cost_concat_1 + seq[(3 * ((s_size) * (j_next - 1) + (k - 1)) + C)] + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)] * cost_concat_2 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + C)] + seq[(3 * ((s_size) * (k_next - 1) + (s_size - 1)) + W)] * cost_concat_3 + seq[(3 * ((s_size) * (k_next - 1) + (s_size - 1)) + C)] -- main.canc.lua:419
if cost_new < cost_best then -- main.canc.lua:421
cost_best = cost_new - EP -- main.canc.lua:422
I = i -- main.canc.lua:423
J = j -- main.canc.lua:424
POS = k -- main.canc.lua:425
end -- main.canc.lua:425
end -- main.canc.lua:425
end -- main.canc.lua:425
if cost_best < solut["cost"] then -- main.canc.lua:432
reinsert(solut["s"], I, J, POS + 1) -- main.canc.lua:436
subseq_load(solut, info) -- main.canc.lua:438
if cost_best ~= solut["cost"] then -- main.canc.lua:441
print("ERROR") -- main.canc.lua:442
os["exit"](1) -- main.canc.lua:443
end -- main.canc.lua:443
return true -- main.canc.lua:445
end -- main.canc.lua:445
return false -- main.canc.lua:448
end -- main.canc.lua:448
RVND = function(solut, info) -- main.canc.lua:451
local SWAP = 0 -- main.canc.lua:452
local REINSERTION = 1 -- main.canc.lua:453
local OR_OPT_2 = 2 -- main.canc.lua:454
local OR_OPT_3 = 3 -- main.canc.lua:455
local TWO_OPT = 4 -- main.canc.lua:456
local neighbd_list = { -- main.canc.lua:458
SWAP, -- main.canc.lua:458
TWO_OPT, -- main.canc.lua:458
REINSERTION, -- main.canc.lua:458
OR_OPT_2, -- main.canc.lua:458
OR_OPT_3 -- main.canc.lua:458
} -- main.canc.lua:458
while # neighbd_list > 0 do -- main.canc.lua:466
local index = math["random"](1, # neighbd_list) -- main.canc.lua:467
index = info["rnd"][info["rnd_index"]] + 1 -- main.canc.lua:469
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:470
local neighbd = neighbd_list[index] -- main.canc.lua:472
local improve = false -- main.canc.lua:474
if neighbd == SWAP then -- main.canc.lua:476
improve = search_swap(solut, info) -- main.canc.lua:477
elseif neighbd == REINSERTION then -- main.canc.lua:479
improve = search_reinsertion(solut, info, REINSERTION) -- main.canc.lua:480
elseif neighbd == OR_OPT_2 then -- main.canc.lua:482
improve = search_reinsertion(solut, info, OR_OPT_2) -- main.canc.lua:483
elseif neighbd == OR_OPT_3 then -- main.canc.lua:485
improve = search_reinsertion(solut, info, OR_OPT_3) -- main.canc.lua:486
elseif neighbd == TWO_OPT then -- main.canc.lua:488
improve = search_two_opt(solut, info) -- main.canc.lua:489
end -- main.canc.lua:489
if improve == true then -- main.canc.lua:494
neighbd_list = { -- main.canc.lua:495
SWAP, -- main.canc.lua:495
TWO_OPT, -- main.canc.lua:495
REINSERTION, -- main.canc.lua:495
OR_OPT_2, -- main.canc.lua:495
OR_OPT_3 -- main.canc.lua:495
} -- main.canc.lua:495
else -- main.canc.lua:495
table["remove"](neighbd_list, index) -- main.canc.lua:498
end -- main.canc.lua:498
end -- main.canc.lua:498
end -- main.canc.lua:498
perturb = function(solut, info) -- main.canc.lua:505
local s = table["clone"](solut["s"]) -- main.canc.lua:506
local A_start = 1 -- main.canc.lua:508
local A_end = 1 -- main.canc.lua:509
local B_start = 1 -- main.canc.lua:510
local B_end = 1 -- main.canc.lua:511
local size_max = math["floor"](# s / 10) -- main.canc.lua:513
size_max = size_max >= 2 and size_max or 2 -- main.canc.lua:514
local size_min = 2 -- main.canc.lua:515
while (A_start <= B_start and B_start <= A_end) or (B_start <= A_start and A_start <= B_end) do -- main.canc.lua:517
A_start = math["random"](2, # s - 1 - size_max) -- main.canc.lua:518
A_end = A_start + math["random"](size_min, size_max) -- main.canc.lua:519
B_start = math["random"](2, # s - 1 - size_max) -- main.canc.lua:521
B_end = B_start + math["random"](size_min, size_max) -- main.canc.lua:522
A_start = info["rnd"][info["rnd_index"]] + 1 -- main.canc.lua:525
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:526
A_end = A_start + info["rnd"][info["rnd_index"]] -- main.canc.lua:527
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:528
B_start = info["rnd"][info["rnd_index"]] + 1 -- main.canc.lua:530
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:531
B_end = B_start + info["rnd"][info["rnd_index"]] -- main.canc.lua:532
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:533
end -- main.canc.lua:533
if A_start < B_start then -- main.canc.lua:536
reinsert(s, B_start, B_end - 1, A_end) -- main.canc.lua:537
reinsert(s, A_start, A_end - 1, B_end) -- main.canc.lua:538
else -- main.canc.lua:538
reinsert(s, A_start, A_end - 1, B_end) -- main.canc.lua:540
reinsert(s, B_start, B_end - 1, A_end) -- main.canc.lua:541
end -- main.canc.lua:541
return s -- main.canc.lua:544
end -- main.canc.lua:544
GILS_RVND = function(Imax, Iils, R, info) -- main.canc.lua:547
local solut_partial = { -- main.canc.lua:549
["s"] = {}, -- main.canc.lua:550
["seq"] = {} -- main.canc.lua:551
} -- main.canc.lua:551
local solut_crnt = { -- main.canc.lua:554
["s"] = {}, -- main.canc.lua:555
["seq"] = {} -- main.canc.lua:556
} -- main.canc.lua:556
local solut_best = { -- main.canc.lua:559
["s"] = {}, -- main.canc.lua:560
["seq"] = {} -- main.canc.lua:561
} -- main.canc.lua:561
subseq_fill(solut_partial["seq"], info) -- main.canc.lua:564
solut_partial["cost"] = 0 -- main.canc.lua:565
subseq_fill(solut_crnt["seq"], info) -- main.canc.lua:567
solut_crnt["cost"] = 0 -- main.canc.lua:568
subseq_fill(solut_best["seq"], info) -- main.canc.lua:570
solut_best["cost"] = 0 -- main.canc.lua:571
solut_best["cost"] = math["huge"] -- main.canc.lua:576
for i = 1, Imax do -- main.canc.lua:578
local Rsz = # R -- main.canc.lua:579
local alpha = R[math["random"](1, Rsz)] -- main.canc.lua:580
alpha = R[info["rnd"][info["rnd_index"]] + 1] -- main.canc.lua:581
info["rnd_index"] = info["rnd_index"] + 1 -- main.canc.lua:583
print("[+] Local Search", i) -- main.canc.lua:586
print("\9[+] Constructing Inital Solution..") -- main.canc.lua:587
solut_crnt["s"] = construction(alpha, info) -- main.canc.lua:588
subseq_load(solut_crnt, info) -- main.canc.lua:589
s_print(solut_crnt) -- main.canc.lua:590
print("\9Construction cost  ", solut_crnt["cost"]) -- main.canc.lua:591
solut_clone(solut_crnt, solut_partial) -- main.canc.lua:597
print("\9[+] Looking for the best Neighbor..") -- main.canc.lua:604
local iterILS = 0 -- main.canc.lua:605
while iterILS < Iils do -- main.canc.lua:606
RVND(solut_crnt, info) -- main.canc.lua:607
if solut_crnt["cost"] < solut_partial["cost"] - info["EPSILON"] then -- main.canc.lua:612
solut_partial["cost"] = solut_crnt["cost"] - info["EPSILON"] -- main.canc.lua:613
solut_partial["s"] = table["clone"](solut_crnt["s"]) -- main.canc.lua:614
iterILS = 0 -- main.canc.lua:615
end -- main.canc.lua:615
solut_crnt["s"] = perturb(solut_partial, info) -- main.canc.lua:618
subseq_load(solut_crnt, info) -- main.canc.lua:619
iterILS = iterILS + 1 -- main.canc.lua:621
end -- main.canc.lua:621
subseq_load(solut_partial, info) -- main.canc.lua:628
if solut_partial["cost"] < solut_best["cost"] then -- main.canc.lua:633
solut_clone(solut_partial, solut_best) -- main.canc.lua:634
end -- main.canc.lua:634
print("\9Current best solution cost: ", solut_best["cost"]) -- main.canc.lua:637
end -- main.canc.lua:637
print("COST: ", solut_best["cost"]) -- main.canc.lua:640
end -- main.canc.lua:640
protect = function(tbl) -- main.canc.lua:643
return setmetatable({}, { -- main.canc.lua:644
["__index"] = tbl, -- main.canc.lua:645
["__newindex"] = function(t, key, value) -- main.canc.lua:646
error("attempting to change constant " .. tostring(key) .. " to " .. tostring(value), 2) -- main.canc.lua:648
end -- main.canc.lua:648
}) -- main.canc.lua:648
end -- main.canc.lua:648
main = function() -- main.canc.lua:653
local info = { -- main.canc.lua:654
["c"] = {}, -- main.canc.lua:655
["T"] = 1, -- main.canc.lua:656
["C"] = 2, -- main.canc.lua:657
["W"] = 3, -- main.canc.lua:657
["EPSILON"] = 1e-15, -- main.canc.lua:658
["rnd"] = {}, -- main.canc.lua:659
["rnd_index"] = 1 -- main.canc.lua:660
} -- main.canc.lua:660
local a = 0 -- main.canc.lua:664
info["dimension"], a = readData(info["c"], info["rnd"]) -- main.canc.lua:665
print(info["rnd"][a]) -- main.canc.lua:666
math["randomseed"](os["time"]()) -- main.canc.lua:667
local Imax = 10 -- main.canc.lua:669
local Iils = math["min"](100, info["dimension"]) -- main.canc.lua:670
local R = { -- main.canc.lua:671
0.00, -- main.canc.lua:671
0.01, -- main.canc.lua:671
0.02, -- main.canc.lua:671
0.03, -- main.canc.lua:671
0.04, -- main.canc.lua:671
0.05, -- main.canc.lua:671
0.06, -- main.canc.lua:671
0.07, -- main.canc.lua:671
0.08, -- main.canc.lua:671
0.09, -- main.canc.lua:671
0.10, -- main.canc.lua:671
0.11, -- main.canc.lua:671
0.12, -- main.canc.lua:671
0.13, -- main.canc.lua:672
0.14, -- main.canc.lua:672
0.15, -- main.canc.lua:672
0.16, -- main.canc.lua:672
0.17, -- main.canc.lua:672
0.18, -- main.canc.lua:672
0.19, -- main.canc.lua:672
0.20, -- main.canc.lua:672
0.21, -- main.canc.lua:672
0.22, -- main.canc.lua:672
0.23, -- main.canc.lua:672
0.24, -- main.canc.lua:672
0.25 -- main.canc.lua:672
} -- main.canc.lua:672
local start = os["clock"]() -- main.canc.lua:676
GILS_RVND(Imax, Iils, R, info) -- main.canc.lua:677
print("TIME: ", os["clock"]() - start) -- main.canc.lua:679
end -- main.canc.lua:679
local profiler = require("profiler") -- main.canc.lua:684
profiler["start"]() -- main.canc.lua:685
main() -- main.canc.lua:687
profiler["stop"]() -- main.canc.lua:688
profiler["report"]("profiler.log") -- main.canc.lua:689
