data = MCMCData(pnts)
@test_nothrow global result = runsampler(data; verbose = false) # test defaults
data = MCMCData(distM)
@test_nothrow global result = runsampler(data; verbose = false) # test defaults