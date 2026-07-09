-- Quantifies total allocation volume (not just live memory) during a
-- real solve of main.canc.lua, to check whether neighbd_list's
-- per-improvement table literal (main.canc.lua:509) is a meaningful GC
-- pressure source, as flagged (but not measured) in IMPLEMENTATION.md.
--
-- Technique: disable automatic GC (collectgarbage("stop")) before
-- running, so nothing gets freed -- collectgarbage("count") growth then
-- equals total bytes allocated (garbage + still-live), not just current
-- live heap. Must run from mlp-main/lua/ (same cwd main.canc.lua expects
-- for its relative dofile/io.open paths).

collectgarbage("collect")
collectgarbage("collect")
local before = collectgarbage("count")

collectgarbage("stop")
dofile("main.canc.lua") -- runs Data.lua load + full GILS_RVND solve, prints COST/TIME itself
local after = collectgarbage("count")

print(("\nTotal allocated during run (GC disabled): %.1f KB (%.1f MB)")
    :format(after - before, (after - before) / 1024))
