dofile("Data.lua")

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

function main() 
    local info = {
        c = {}, 
        T = 1,
        C = 2,
        W = 3
    }

    info.dimension = readData(info.c)

    print(info.dimension)
    for i=1,info.dimension do
        for j=1,info.dimension do
            io.write(info.c[i][j], " ")
        end
        print()
    end

    local solut = {
        s = {},
        seq = {}
    }

    subseq_fill(solut.seq, info)

    for i=1,info.dimension+1 do
        solut.seq[i] = {}
        for j=1,info.dimension+1 do
            solut.seq[i][j] = {0, 0, 0}
            for s=1,3 do
                io.write(solut.seq[i][j][s] , " ") --= {0, 0, 0}
            end
            io.write("| ")
        end
        print()
    end

    for i=1,info.dimension do
        solut.s[#solut.s + 1] = i
        io.write(solut.s[i], " ")
    end
    solut.s[#solut.s + 1] = 1
    io.write(solut.s[info.dimension+1], " ")
    print()



    subseq_load(solut, info)

    print(solut.seq[1][info.dimension][info.C])

end

main()

--[[
--comments
]]--

