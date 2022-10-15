data = MCMCData(pnts)
@test_nothrow global result = runsampler(data; verbose = false) # test defaults
data = MCMCData(distM)
@test show(data) # test pretty printing
@test display(data) # test multiline pretty printing
@test_nothrow global result = runsampler(data; verbose = false) # test defaults
options = MCMCOptionsList(numMH = 0)
@test_nothrow show(options) # test pretty printing
@test_nothrow display(options) # test multiline pretty printing
@test_nothrow global result = runsampler(data, options; verbose = false) # test pure Gibbs sampling
@test_nothrow show(result) # test pretty printing
@test_nothrow display(result) # test multiline pretty printing