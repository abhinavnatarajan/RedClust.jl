@test_nothrow fitprior(pnts, "k-means", false; useR = true, verbose = false) # kmeans with points, R
@test_nothrow fitprior(pnts, "k-medoids", false; useR = true, verbose = false) # kmedoids with points, R
@test_nothrow fitprior(distM, "k-medoids", true; useR = true, verbose = false) # kmedoids with distances, R
@test_nothrow fitprior(pnts, "k-means", false; useR = false, verbose = false) # kmeans with points, Julia
@test_nothrow fitprior(pnts, "k-medoids", false; useR = false, verbose = false) # kmedoids with points, Julia
@test_nothrow fitprior(distM, "k-medoids", true; useR = false, verbose = false) # kmedoids with distances, Julia
params = ()
@test_logs (:warn,) global params = fitprior(pnts, "k-means", false; useR = false, verbose = false, Kmin = 1, Kmax = 1) # check edge case Kmax = 2
@test params.K_initial == 1
@test_logs (:warn,) global params = fitprior(pnts, "k-means", false; useR = false, verbose = false, Kmin = N, Kmax = N) # check edge case Kmin = N-1
@test params.K_initial == N
@test_throws ArgumentError fitprior(distM, "hierarchical", true; useR = false, verbose = false) # check that only 2 methods allowed
@test_throws ArgumentError fitprior(pnts, "k-means", true; useR = false, verbose = false) # check that dist is disregarded when pnts is supplied
@test_throws ArgumentError fitprior(distM, "k-means", true; useR = false, verbose = false) # check that kmeans does not work when distances are given
@test_throws ArgumentError fitprior(pnts, "k-means", false; Kmin = 0, verbose = false) # check 1 ≤ Kmin
@test_throws ArgumentError fitprior(pnts, "k-means", false; Kmin = N, Kmax = 1, verbose = false) # check Kmin ≤ Kmax
@test_throws ArgumentError fitprior(pnts, "k-means", false; Kmax = length(pnts) + 1, verbose = false) # check Kmax ≤ N