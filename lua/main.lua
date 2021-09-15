dofile("Data.lua")

function table.clone(org)
      return {table.unpack(org)}
  end

function subseq_fill(seq, info)
    for i=1,info.dimension+1 do
        seq[i] = {}
        for j=1,info.dimension+1 do
            seq[i][j] = {0, 0, 0}
        end
    end
end

function subseq_load(solut, info)
    for i=1,info.dimension+1 do
        local k = 1 - i - (i ~= 1 and 0 or 1)

        solut.seq[i][i][info.T] = 0.0
        solut.seq[i][i][info.C] = 0.0
        solut.seq[i][i][info.W] = (i ~= 1 and 1 or 0)
        for j=i+1,info.dimension+1 do
            local j_prev = j-1

            solut.seq[i][j][info.T] = info.c[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev][info.T]
            solut.seq[i][j][info.C] = solut.seq[i][j][info.T] + solut.seq[i][j_prev][info.C]
            solut.seq[i][j][info.W] = j + k
        end
    end
end

function construction(alpha, info) 
    math.randomseed(os.time()) 
    local s = {1}

    local cList = {}
    for i=2,info.dimension do
        cList[#cList + 1] = i
        io.write(cList[i-1], " ")
    end
    print()

    local r = 1

    while #cList > 0 do
    --while table.getn(cList) > 0 do
        table.sort(cList, function(i, j) return info.c[r][i] < info.c[r][j] end)
        io.write("Soted ")
        --[[
        for i=1,#cList do
            io.write(cList[i], " ")
        end
        print()
        ]]--

        local range = #cList * alpha + 1
        local i = math.random(1, range)

        print(#cList)
        local c = table.remove(cList, i)
        table.insert(s, c)
        r = c
    end
    table.insert(s, 1)

    for i=1,#s do
        io.write(s[i], " ")
    end
    print()

    return total.clone(s)
end

function RVND(solut, info)
end

function GILS_RVND(Imax, Iils, R, info)

    local solut_partial = {
        sl = {},
        seq = {}
    }

    local solut_crnt = {
        s = {},
        seq = {}
    }

    local solut_best = {
        sl = {},
        seq = {}
    }

    subseq_fill(solut.seq, info)

    for i=1,Imax do
        local Rsz = #R
        local alpha = R[math.random(1, Rsz)]

        solut.s = construction(alpha, info)
        solut_partial.sl = total.clone(solut.s)

        subseq_load(solut.s, info)

        local rnvd_cost_best = solut.seq[1][info.dimension+1][info.C]

        local iterILS = 0
        while iterILS < Iils do
            RVND(solut, info)
            local rnvd_cost_crnt = solut.seq[1][info.dimension+1][info.C]
            if rnvd_cost_crnt < rvnd_cost_best then
               rnvd_cost_best = rnvd_cost_crnt - info.EPSILON
               solution_best.sl = table.clone(solut.s)
               iterILS = 0
            end
            iterILS = iterILS + 1
        end

        subseq_load(solut_partial, info)
        local cost_partial = solut_partial.seq[1][info.dimension+1][info.C]

        if cost_partial < cost_best then
            cost_best  = cost_partial
            solut_best = table.clone(solut_partial)
        end
    end
end

function main() 
    local info = {
        c = {}, 
        T = 1,
        C = 2,
        W = 3,
        EPSILON = 1^-15
    }


    info.dimension = readData(info.c)

    print(info.dimension)
    for i=1,info.dimension do
        for j=1,info.dimension do
            io.write(info.c[i][j], " ")
        end
        print()
    end


    --[[
    for i=1,info.dimension+1 do
        for j=1,info.dimension+1 do
            for s=1,3 do
                io.write(solut.seq[i][j][s] , " ") --= {0, 0, 0}
            end
            io.write("| ")
        end
        print()
    end

    for i=1,info.dimension do
        solut.sl[#solut.sl + 1] = i
        io.write(solut.sl[i], " ")
    end
    solut.sl[#solut.sl + 1] = 1
    io.write(solut.sl[info.dimension+1], " ")
    print()
    ]]--

    Imax = 10
    Iils = math.min(100, info.dimension)
    R = = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
            0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25}

    GILS_RVND(Imax, Iils, R, info)

end

main()
