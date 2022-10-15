# Test point estimation
temp = sample(1:K, N)
@test infodist(temp, temp; normalised = true) ≈ 0 atol = 1e-9
@test_throws ArgumentError infodist(temp, [temp; 1]) # check that lengths are equal
@test infodist(temp, temp; normalised = false) ≈ 0 atol = 1e-9
@test binderloss(temp, temp; normalised = true) ≈ 0 atol = 1e-9
@test_throws ArgumentError binderloss(temp, [temp; 1]) # check that lengths are equal
@test binderloss(temp, temp; normalised = false) ≈ 0 atol = 1e-9
# test the methods that should work
@test_nothrow getpointestimate(result; method = "MAP")
@test_nothrow getpointestimate(result; method = "MLE")
@test_nothrow getpointestimate(result; loss = "binder", method = "MPEL")
@test_nothrow getpointestimate(result; loss = "omARI", method = "MPEL")
@test_nothrow getpointestimate(result; loss = "VI", method = "MPEL")
@test_nothrow getpointestimate(result; loss = "ID", method = "MPEL")
@test_nothrow getpointestimate(result; loss = "binder", method = "SALSO")
@test_nothrow getpointestimate(result; loss = "omARI", method = "SALSO")
@test_nothrow getpointestimate(result; loss = "VI", method = "SALSO")
@test_nothrow getpointestimate(result; loss = "ID", method = "SALSO")
@test_nothrow getpointestimate(result; loss = varinfo, method = "MPEL")
# test bad inputs
@test_throws ArgumentError getpointestimate(result; loss = (x, y) -> 1.0, method = "SALSO")
@test_throws ArgumentError getpointestimate(result; loss = "ID", method = "some other method")
@test_throws ArgumentError getpointestimate(result; loss = "some other loss function", method = "SALSO")
@test_throws ArgumentError getpointestimate(result; loss = "some other loss function", method = "MPEL")