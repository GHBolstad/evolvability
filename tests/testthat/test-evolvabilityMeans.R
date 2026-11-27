test_that("an identity G has unit average parameters", {
  expect_equal(evolvabilityMeans(diag(3)),
               c(e_mean = 1, e_min = 1, e_max = 1, c_mean = 1, r_mean = 1,
                 a_mean = 1, i_mean = 0))
})

test_that("evolvabilityMeansMCMC agrees with evolvabilityMeans", {
  G1 <- crossprod(matrix(c(1, 2, 0, 1, 3, 1, 0, 1, 2), 3))
  G2 <- G1 + diag(3)
  res <- evolvabilityMeansMCMC(rbind(c(G1), c(G2)))
  expect_equal(unname(res$post.dist[1, ]),
               unname(evolvabilityMeans(G1)[colnames(res$post.dist)]))
  expect_equal(unname(res$post.dist[2, ]),
               unname(evolvabilityMeans(G2)[colnames(res$post.dist)]))
})
