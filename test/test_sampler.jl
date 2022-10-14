data = MCMCData(pnts)
@test_nothrow global result = runsampler(data; verbose = false) # test defaults
data = MCMCData(distM)
@test_nothrow global result = runsampler(data; verbose = false) # test defaults
@test_nothrow global result = runsampler(data, MCMCOptionsList(numMH = 0); verbose = false) # test pure Gibbs sampling