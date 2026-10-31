test_that("Davies and Liu quantiles match a noncentral chi-square across units", {
    tree <- ape::read.tree(text = "(a:1,b:1,c:1,d:1);")
    y <- c(a = -1, b = 0.5, c = 2, d = 3)
    probs <- c(.25, .5, .75)
    exact <- .25 * stats::qchisq(
        probs, df = 3, ncp = sum((y - mean(y))^2)
    ) / 3
    for (method in c("davies", "liu")) {
        for (scale in c(1e-7, 1e-3, 1, 1e7)) {
            result <- reconstruct_empirical_moments(
                tree, y * scale, times = .5, sigma2 = scale^2,
                quantities = "variance", variance_probs = probs,
                quantile_method = method
            )
            actual <- unlist(result$data[c("q25","q50","q75")], use.names=FALSE)
            expect_equal(actual / scale^2, exact, tolerance = 2e-6)
        }
    }
})

test_that("rank-one quantiles retain precision for small trait units", {
    for (scale in c(1e-10, 1, 1e10)) {
        for (prob in c(.025, .5, .975)) {
            actual <- .empvar_quantile_rank1(
                1-prob, lambda = .7, eta = 1.2*scale, sigma2 = scale^2
            )
            reference <- .7 * stats::qchisq(prob, 1, ncp = 1.2^2)
            expect_equal(actual/scale^2, reference, tolerance = 1e-10)
        }
    }
})


test_that("direct Liu inversion retains the CompQuadForm approximation", {
    lambda <- c(.1,.8,2)
    eta <- c(1,-.5,.2)
    sigma2 <- .7
    for (prob in c(.025,.25,.5,.75,.975)) {
        q <- .empvar_quantile_liu(1-prob,lambda,eta,sigma2)
        expect_equal(
            CompQuadForm::liu(q/sigma2,lambda,delta=eta^2/sigma2),
            1-prob,tolerance=1e-9
        )
    }
})

test_that("quantile integration failures cannot silently become quantiles", {
    testthat::local_mocked_bindings(
        .quadratic_form_tail_davies = function(...) {
            list(probability=-1, convergence=FALSE, diagnostic="forced failure")
        },
        .package = "PastMoments"
    )
    reference <- stats::qchisq(.5, 2)
    expect_equal(.empvar_quantile_davies(.5, c(1,1), c(0,0), 1),
                 reference, tolerance = 2e-6)
    testthat::local_mocked_bindings(
        .quadratic_form_tail_imhof = function(...) {
            list(probability=NA_real_, convergence=FALSE, diagnostic="forced failure")
        },
        .package = "PastMoments"
    )
    expect_error(.empvar_quantile_davies(.5, c(1,1), c(0,0), 1),
                 "did not converge")
})

test_that("time units preserve covariance, spread, quantiles and disparity", {
    tree <- example_tree()
    trait <- example_trait()
    for (method in c("qapprox", "davies", "liu")) {
        if (method == "qapprox" && !requireNamespace("Qapprox", quietly=TRUE)) next
        reference <- reconstruct_empirical_moments(
            tree, trait, times=c(.5,1.5,2), sigma2=1,
            variance_probs=c(.25,.5,.75), quantile_method=method
        )
        for (scale in c(1e-13,1e13)) {
            scaled_tree <- tree
            scaled_tree$edge.length <- tree$edge.length * scale
            result <- reconstruct_empirical_moments(
                scaled_tree, trait, times=c(.5,1.5,2)*scale, sigma2=1/scale,
                variance_probs=c(.25,.5,.75), quantile_method=method
            )
            expect_equal(result$data[names(result$data)!="time"],
                         reference$data[names(reference$data)!="time"],
                         tolerance=2e-6)
            tip <- result$data[result$data$time==2*scale &
                               result$data$quantity=="empirical_variance", ]
            expect_equal(tip$q50, stats::var(trait), tolerance=1e-10)
        }
    }
})
