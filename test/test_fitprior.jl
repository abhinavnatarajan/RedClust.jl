@test_nothrow fitprior(pnts, "k-means", false; verbose = false) # kmeans with points, Julia
@test_nothrow fitprior(pnts, "k-medoids", false; verbose = false) # kmedoids with points, Julia
@test_nothrow fitprior(distM, "k-medoids", true; verbose = false) # kmedoids with distances, Julia
@test_logs (:warn,) global params = fitprior(pnts, "k-means", false; verbose = false, Kmin = 1, Kmax = 1) # check edge case Kmax = 2
@test params.K_initial == 1
@test_logs (:warn,) global params = fitprior(pnts, "k-means", false; verbose = false, Kmin = N, Kmax = N) # check edge case Kmin = N-1
@test params.K_initial == N
@test_throws ArgumentError fitprior(distM, "hierarchical", true; verbose = false) # check that only 2 methods allowed
@test_throws ArgumentError fitprior(pnts, "k-means", true; verbose = false) # check that dist is disregarded when pnts is supplied
@test_throws ArgumentError fitprior(distM, "k-means", true; verbose = false) # check that kmeans does not work when distances are given
@test_throws ArgumentError fitprior(pnts, "k-means", false; Kmin = 0, verbose = false) # check 1 ≤ Kmin
@test_throws ArgumentError fitprior(pnts, "k-means", false; Kmin = N, Kmax = 1, verbose = false) # check Kmin ≤ Kmax
@test_throws ArgumentError fitprior(pnts, "k-means", false; Kmax = length(pnts) + 1, verbose = false) # check Kmax ≤ N
@test_nothrow show(params)
@test_nothrow display(params)

@test_nothrow fitprior2(pnts, "k-means", false; verbose = false) # kmeans with points, Julia
@test_nothrow fitprior2(pnts, "k-medoids", false; verbose = false) # kmedoids with points, Julia
@test_nothrow fitprior2(distM, "k-medoids", true; verbose = false) # kmedoids with distances, Julia
@test_logs (:warn,) global params = fitprior(pnts, "k-means", false; verbose = false, Kmin = 1, Kmax = 1) # check edge case Kmax = 2
@test params.K_initial == 1
@test_logs (:warn,) global params = fitprior(pnts, "k-means", false; verbose = false, Kmin = N, Kmax = N) # check edge case Kmin = N-1
@test params.K_initial == N